#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import torch
from torch import tensor as tt
import scanpy as sc
import anndata
from sklearn.decomposition import PCA
from noise_model import Noise_Model
from collections import namedtuple

####################################################################################################
# big functions apparing in main.py


def get_data_from_anndata(path, gene_list=None, cell_list=None):
    """
    This function is used to load the data from the anndata object.
    It also performs some filtering and normalization. In this case we get rid of a specific sample,
    filters too low mito content.
    """
    adatas = {"Tomaz": sc.read_h5ad(path)}
    data = anndata.concat(adatas, label="dataset", join="inner")
    sample_names = data.obs["orig.ident"]

    # removing the bad replicate
    remove = "scRNA_Seq_Tomaz_220531-01-M-ZT22-H-lib_1-129S"
    remove_mask = np.invert(sample_names == remove)

    # update data and sample list after the removal of the sample
    data = data[remove_mask, :]

    # filter mitochoncdrial genes
    genes_mito = [s.startswith("mt-") for s in data.var.index]
    genes_Mup = [s.startswith("Mup") for s in data.var.index]
    genes_keep = np.logical_not(np.logical_or(genes_mito, genes_Mup))
    # cut out cells with an unusual number of reads + too low mito content
    cells_keep = np.logical_and(
        data.obs["nCount_RNA"] < 15000, data.obs["percent.mt"] > 3
    )
    data = data[cells_keep, genes_keep]

    # adding fields to the data object
    data.layers["n_cg"] = data.X.toarray()
    data.obs["n_c"] = data.layers["n_cg"].sum(axis=1)
    data.layers["f_cg"] = data.layers["n_cg"] / data.obs["n_c"][:, None]

    sample_names = data.obs["orig.ident"]
    nn = sample_names.unique()

    sample_id = np.zeros(data.n_obs, dtype=np.int64)
    for i, s in enumerate(nn):
        sample_id[sample_names == s] = i
    dm = make_design_matrix(torch.tensor(sample_id, dtype=float))

    data_zonated = data[:, gene_list]
    return data, data_zonated, data.obs["n_c"].values, dm, sample_id, nn


def do_pca(data, pc=0):
    """
    in this fuynction we do pca on the data.
    In our datset the zonation source of variation is clrearly the zonation one,
    however in other datasets, could not be the same case.
    """

    pca = PCA(n_components=5, whiten=False)
    X_pca = data

    # normalize yourself
    X_pca = X_pca - X_pca.mean(axis=0)[None, :]
    X_pca = X_pca / np.sqrt(X_pca.var(axis=0))[None, :]
    PC = pca.fit_transform(X_pca)

    # normalization of the first PC
    x_unif = PC[:, pc]
    x_unif = x_unif - x_unif.mean()
    x_unif = x_unif / np.sqrt(x_unif.var())
    return x_unif


def fit_coeff(data, x_unif, genes):
    """
    This function is used to fit the coefficients a0 and a1 of the model
    using a negative binomial regression.
    """
    noise = "NB"
    D = np.stack((np.repeat(1, data.n_obs), x_unif), axis=1)
    coef_pau = np.zeros((len(genes), D.shape[1]))

    alpha = np.zeros(len(genes))
    logN = np.log(data.obs["n_c"].values)

    for gi, g in enumerate(genes):

        yy = data[:, [g]].layers["n_cg"].toarray().squeeze()
        model_n = Noise_Model(yy, D, logN, noise)

        iterations = 50
        mf = model_n.fit(iterations)
        if noise == "Poisson":
            coef_pau[gi, :] = mf
        else:
            coef_pau[gi, :] = mf[:-1]
            alpha[gi] = mf[-1]
    return coef_pau


def training(data, x_unif, coef_pau, n_c, dm, clamp, n_iter, batch_size, dev):
    """
    This function is used to train the model on the data.
    It's the only function where pytorch is used.
    """
    # set the device to use
    torch.set_default_device(dev)
    DATA = torch.tensor(data)
    NC, NG = DATA.shape
    NS = dm.shape[1]
    if batch_size == 0:
        batch_size = NC

    # preparing the starting values for the optimization
    a0_pau = coef_pau[:, 0]
    a1_pau = coef_pau[:, 1]
    scale = a1_pau[clamp]
    a1_scaled = a1_pau / scale
    x_scaled = x_unif * scale

    mask = torch.eye(NG, dtype=float)
    mask[clamp, clamp] = 0
    fix = tt(1.0).detach()
    log_n_UMI = torch.log(tt(n_c))
    mp = dict(
        log_n_UMI=log_n_UMI,
        clamp=clamp,
        dm=dm,
        fix=fix,
        mask=mask,
        cutoff=50,
    )
    MyTuple = namedtuple("param", mp)
    MP = MyTuple(**mp)

    # initalizing the parameters (leafs)
    disp = tt(np.log(0.3), requires_grad=True)
    x = tt(x_scaled, requires_grad=True, dtype=float)
    a1 = tt(a1_scaled, requires_grad=True, dtype=float)
    a0 = tt(a0_pau, dtype=torch.float32)
    a0 = a0.repeat(NS, 1)
    a0.requires_grad = True

    # training the model
    losses = []
    optimizer = torch.optim.Adam([x, a0, a1, disp], lr=0.001)
    batch_size = NC
    # Optimize the latent variables to minimize the loss
    for step in range(n_iter):
        optimizer.zero_grad()  # zero the gradients
        output = loss_clamp_batch(x, a0, a1, disp, batch_size, MP, DATA)
        output.backward()  # compute the gradients
        optimizer.step()  # update the variable
        losses.append(output.detach())

    x_final = x.clone().cpu().detach().numpy()
    disp_final = torch.exp(disp.clone()).cpu().detach().numpy()
    a0_final = a0.clone().cpu().detach().numpy()
    a1_final = a1.clone().cpu().detach().numpy()

    return x_final, a0_final, a1_final, disp_final, losses


def training_gene_selection(
    DATA, x_unif, coef_pau, n_c, dm, clamp, n_iter, batch_size, dev
):
    """
    This function is used to train the model on the data.
    """
    NC, NG = DATA.shape
    NS = dm.shape[1]
    if batch_size == 0:
        batch_size = NC

    # preparing the starting values for the optimization
    a0_pau = coef_pau[:, 0]
    a1_pau = coef_pau[:, 1]
    scale = a1_pau[clamp]
    a1_scaled = a1_pau / scale
    x_scaled = x_unif * scale

    mask = torch.eye(NG, device=dev, dtype=float)
    mask[clamp, clamp] = 0
    fix = tt(1.0, device=dev).detach()
    log_n_UMI = torch.log(tt(n_c, device=dev))
    mp = dict(
        log_n_UMI=log_n_UMI,
        clamp=clamp,
        dm=dm,
        fix=fix,
        mask=mask,
        cutoff=50,
    )
    MyTuple = namedtuple("param", mp)
    MP = MyTuple(**mp)

    # initalizing the parameters (leafs)
    disp = tt(np.log(0.3), requires_grad=True, device=dev)
    x = tt(x_scaled, requires_grad=True, dtype=float, device=dev)
    a1 = tt(a1_scaled, requires_grad=True, dtype=float, device=dev)
    a0 = tt(a0_pau, dtype=torch.float32, device=dev)
    a0 = a0.repeat(NS, 1)
    a0.requires_grad = True

    # training the model
    losses = []
    optimizer = torch.optim.Adam([x, a0, a1, disp], lr=0.001)
    batch_size = NC
    # Optimize the latent variables to minimize the loss
    for step in range(n_iter):
        optimizer.zero_grad()  # zero the gradients
        output = loss_clamp_batch(x, a0, a1, disp, batch_size, MP, DATA)
        output.backward()  # compute the gradients
        optimizer.step()  # update the variable
        losses.append(output.detach())

    x_final = x.clone().cpu().detach().numpy()
    disp_final = torch.exp(disp.clone()).cpu().detach().numpy()
    a0_final = a0.clone().cpu().detach().numpy()
    a1_final = a1.clone().cpu().detach().numpy()

    return x_final, a0_final, a1_final, disp_final, losses


def save_parameters(x, a0, a1, sample_names_uniq, name, data):
    """
    This function is used to save the parameters to a .txt file.
    """
    df_a0 = pd.DataFrame(a0, index=sample_names_uniq, columns=data.var.index)
    df_a1 = pd.Series(a1, index=data.var.index)
    df_x = pd.Series(x, index=data.obs.index)

    # Write the DataFrame to a .txt file
    df_a0.to_csv("coeff_values/" + name + "_a0.txt", sep=",")
    df_a1.to_csv("coeff_values/" + name + "_a1.txt", sep=",")
    df_x.to_csv("coeff_values/" + name + "_x.txt", sep=",")


def shift_samples_per_mouse(x, a0, a1, sample_id, central, data):
    """
    This function is used to shift the samples per mouse.
    This is done because diferent sample (mice) have different
    are x values, and different mice look completely shifted to one another.
    As this is not desirable, we alligned our sumple such as the are aligned
    on the central (defualt) side.
    """
    xs = []
    # computing the individual shifts
    for s in np.unique(sample_id):
        idx = sample_id == s
        central_score = data[idx, central].layers["f_cg"].sum(1)
        xx = x[idx]
        rc = central_score.argsort().argsort() / idx.sum()
        col = np.logical_and(rc > 0.9, rc < 0.99)
        offset = np.median(xx[col])
        xs.append(offset)

    x_shifted, a0_shifted, a1_shifted = shift_parameters(x, a0, a1, xs, sample_id)

    return x_shifted, a0_shifted, a1_shifted, xs


def shift_parameters(x, a0, a1, shift, sample_id):
    """
    This function transforms all the model's parameters
    accordingly to the shift that were found beforehand.
    In particular it shifts the x, and trasform the intercept a0,
    accordingly to the shift.
    """

    x_scaled = np.zeros(len(x), dtype=float)
    a0_scaled = np.zeros(a0.shape, dtype=float)
    a1_scaled = a1
    for s in np.unique(sample_id):
        idx = sample_id == s
        xx = x[idx]
        x_range = 1.0
        x_scaled[idx] = (xx - shift[s]) / x_range
        a0_scaled[s, :] = a0[s, :] + a1[:] * shift[s]
        a1_scaled[:] = a1[:] * x_range
    return x_scaled, a0_scaled, a1_scaled


####################################################################
# smaller functions, called many times


def loss_clamp_batch(x, a0, a1, disp, batch_size, mp, DATA):
    """
    function that takes as input the leaf parameters and returns the loss.
    Beware disp will need to be exp() in the main function
    a0 is sample specific, disp and a1 only gene specific.
    If you want to use all datapoints, set batch_size = DATA.shape[0]
    The 'clamp' gene slope coefficient a1 is set to the fix value (1).

    """
    NC = DATA.shape[0]
    # killing the gradient
    a1_ = torch.matmul(mp.mask, a1)
    a1_[mp.clamp] = mp.fix

    idx = torch.randperm(DATA.shape[0])[:batch_size]
    y = x[idx, None] * a1_[None, :] + mp.log_n_UMI[idx, None]
    y += torch.matmul(mp.dm[idx, :], a0)
    alpha = torch.exp(disp)

    y = mp.cutoff * torch.tanh(y / mp.cutoff)
    lmbda = torch.exp(y)

    r = 1 / alpha
    p = alpha * lmbda / (1 + alpha * lmbda)
    NB = torch.distributions.NegativeBinomial(
        total_count=r, probs=p, validate_args=None
    )
    return -NB.log_prob(DATA[idx, :]).sum() * (
        NC / batch_size
    )  # correct the likelihood rescaling


def loss_simple(x, a0, a1, disp, mp, DATA):

    NC = DATA.shape[0]
    # killing the gradient
    a1_ = torch.matmul(mp.mask, a1)
    a1_[mp.clamp] = mp.fix

    # idx = torch.randperm(DATA.shape[0])[:batch_size]
    y = x[:, None] * a1_[None, :] + a0[None, :] + mp.log_n_UMI[:, None]
    alpha = torch.exp(disp)

    y = mp.cutoff * torch.tanh(y / mp.cutoff)
    lmbda = torch.exp(y)

    r = 1 / alpha
    p = alpha * lmbda / (1 + alpha * lmbda)
    NB = torch.distributions.NegativeBinomial(
        total_count=r, probs=p, validate_args=None
    )
    return -NB.log_prob(DATA[:, :]).sum()


def loss_simple_batch(x, a0, a1, disp, batch_size, mp, DATA):

    NC = DATA.shape[0]
    # killing the gradient
    a1_ = torch.matmul(mp.mask, a1)
    a1_[mp.clamp] = mp.fix

    idx = torch.randperm(DATA.shape[0])[:batch_size]
    y = x[idx, None] * a1_[None, :] + a0[None, :] + mp.log_n_UMI[idx, None]
    alpha = torch.exp(disp)

    y = mp.cutoff * torch.tanh(y / mp.cutoff)
    lmbda = torch.exp(y)

    r = 1 / alpha
    p = alpha * lmbda / (1 + alpha * lmbda)
    NB = torch.distributions.NegativeBinomial(
        total_count=r, probs=p, validate_args=None
    )
    return -NB.log_prob(DATA[idx, :]).sum() * (NC / batch_size)


def loss_gene_selection(x, a0, a1, disp, batch_size, mp, DATA):
    """
    This loss function does not optimize x, but only a0 and a1.
    It is used copute the likelihood of the data, and fits the
    parameters than will later be used for the optimization of x.
    """
    NC = DATA.shape[0]

    idx = torch.randperm(DATA.shape[0])[:batch_size]
    y = (
        mp.log_n_UMI[idx, None]
        + torch.matmul(mp.dm[idx, :], a0)
        + torch.matmul(mp.dm[idx, :], a1) * x[idx, None]
    )
    alpha = disp.exp()

    y = mp.cutoff * torch.tanh(y / mp.cutoff)
    lmbda = torch.exp(y)

    r = 1 / alpha
    p = alpha * lmbda / (1 + alpha * lmbda)
    NB = torch.distributions.NegativeBinomial(
        total_count=r, probs=p, validate_args=None
    )
    return -NB.log_prob(DATA[idx, :]).sum() * (
        NC / batch_size
    )  # correct the likelihood rescaling


def make_design_matrix(cell_identifiers):
    """A function to create a generic design matrix from cell identifiers

    For example:
    input, cell_identifiers = [1,1,1,2,2,2]
    output, design_matrix: [[1,1,1,0,0,0][0,0,0,1,1,1]]

    Arguments
    ---------
    cell_identifiers: torch tensor
        List of cells identifiers, should be of type torch.int

    Returns
    -------
    design_matrix: torch tensor
        Design matrix of shape (num_cells, num_unique_identifiers)

    """
    design_matrix = torch.hstack(
        [
            (cell_identifiers == v).type(torch.float).reshape(len(cell_identifiers), 1)
            for i, v in enumerate(torch.unique(cell_identifiers))
        ]
    )
    return design_matrix


def scale_parameters(x, a0, a1):
    x_range = x.clone().max() - x.clone().min()
    x_min = x.clone().min()
    x_scaled = (x.clone() - x_min) / x_range
    a0_scaled = a0.clone() + a1.clone() * x_min
    a1_scaled = a1.clone() * x_range
    return x_scaled, a0_scaled, a1_scaled


def scale_parameters2(x, a0, a1, x_min, x_max):
    x_range = x_max - x_min
    x_scaled = (x - x_min) / x_range
    a0_scaled = a0 + a1 * x_min
    a1_scaled = a1 * x_range
    return x_scaled, a0_scaled, a1_scaled


def gene_index(data, gene):
    return np.where(gene == data.var.index)[0][0]


def import_data_txt(location):
    df = pd.read_csv(location, delimiter="\t")
    z = df.values
    return z


# functions use for plotting results

# def f_norm(x, a0_pyro, a1_pyro, DM):
#     y = x[:, None] * a1_pyro[None, :]
#     y += np.matmul(DM[:, :], a0_pyro)
#     return np.exp(y)


# def f_single(x, a0_pyro, a1_pyro, DM):
#     y = x[:, None] * np.matmul(DM[:, :], a1_pyro)
#     y += np.matmul(DM[:, :], a0_pyro)
#     return np.exp(y)


# def f_norm_a0_mean(x, a0_pyro, a1_pyro):
#     y = x[:, None] * a1_pyro[None, :] + a0_pyro.mean(axis=0)
#     return np.exp(y)


# def split_vector_into_subsets(vector, num_subsets):
#     avg = len(vector) // num_subsets
#     remainder = len(vector) % num_subsets
#     subsets = []
#     i = 0
#     for _ in range(num_subsets):
#         subset_size = avg + (1 if remainder > 0 else 0)
#         subset = vector[i : i + subset_size]
#         subsets.append(subset)
#         i += subset_size
#         remainder -= 1
#     return subsets
