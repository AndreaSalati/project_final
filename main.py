import config as cfg
from utils import *
from torch_losses import *
from noise_model import Noise_Model

import numpy as np
import pandas as pd
import torch
from torch import tensor as tt

import matplotlib.pyplot as plt

data_0, data, n_c, dm = get_data_from_anndata(cfg.path, cfg.genes)

x_unif = do_pca(data.layers["f_cg"])

coef_pau = fit_coeff(data, x_unif, cfg.genes)

DATA = torch.tensor(data.layers["n_cg"], device=cfg.dev)
NC, NG = DATA.shape
NS = dm.shape[1]

# here we decide the gene to clamp in order to get rid of one of the likelihood symmetries
clamp = gene_index(data, "Cyp2e1")

x_final, disp_final, a0_final, a1_final, losses = training(
    DATA, x_unif, coef_pau, n_c, dm, clamp, cfg.n_iter, cfg.dev
)


plt.scatter(x_final, x_unif[:])
plt.show()
