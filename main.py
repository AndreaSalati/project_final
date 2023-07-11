import config as cfg
from utils import *
from torch_losses import *
from noise_model import Noise_Model

import numpy as np
import pandas as pd

# import torch
# from torch import tensor as tt

import matplotlib.pyplot as plt

data_0, data, n_c, dm, sample_id, sample_names = get_data_from_anndata(
    cfg.path, cfg.genes
)
x_unif = do_pca(data.layers["f_cg"])

coef_pau = fit_coeff(data, x_unif, cfg.genes)

# DATA = torch.tensor(data.layers["n_cg"], device=cfg.dev)

NC, NG = data.shape
NS = dm.shape[1]
print(NC, NG)
# here we decide the gene to clamp in order to get rid of one of the likelihood symmetries
clamp = gene_index(data, cfg.clamp_gene)

x_, a0_, a1_, disp, losses = training(
    data.layers["n_cg"],
    x_unif,
    coef_pau,
    n_c,
    dm,
    clamp,
    cfg.n_iter,
    cfg.batch_size,
    cfg.dev,
)


plt.scatter(x_, x_unif[:], c=sample_id, s=1)
plt.show()

x_shifted, a0_shifted, a1_shifted, xs = shift_samples_per_mouse(
    x_, a0_, a1_, sample_id, cfg.central, data
)

# save_parameters(
#     x_shifted,
#     a0_shifted,
#     a1_shifted,
#     sample_names,
#     cfg.name,
#     data,
# )
