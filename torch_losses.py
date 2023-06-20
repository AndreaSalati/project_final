# import numpy as np
import torch
from torch import tensor as tt


def loss_clamp(x, a0, a1, disp, mp, DATA):

    NC = DATA.shape[0]
    # killing the gradient for the clamp entry
    a1_ = torch.matmul(mp.mask, a1)
    a1_[mp.clamp] = mp.fix

    y = x[:, None] * a1_[None, :] + mp.log_n_UMI[:, None]
    y += torch.matmul(mp.dm[:, :], a0)
    alpha = disp.abs()

    y = mp.cutoff * torch.tanh(y / mp.cutoff)

    lmbda = torch.exp(y)
    r = 1 / alpha
    p = alpha * lmbda / (1 + alpha * lmbda)
    NB = torch.distributions.NegativeBinomial(
        total_count=r, probs=p, validate_args=None
    )
    return -NB.log_prob(DATA).sum()


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
