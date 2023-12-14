path = "data/130723_data_hepatocytes.h5ad"

# small set of genes
central = [
    "Oat",
    "Cyp2e1",
    "Lect2",
    "Cyp2c37",
    "Gulo",
    "Cyp2a5",
    "Glul",
    "Aldh1a1",
    "Cyp1a2",
    "Slc22a1",
    "Slc1a2",
]
portal = ["Pck1", "Aldh1b1", "Ctsc", "Sds", "Hal", "Hsd17b13", "Cyp2f2"]


# full set of genes
# central = [
#     "Rnase4",
#     "Glul",
#     "Oat",
#     "Cyp2c37",
#     "Airn",
#     "Lect2",
#     "Slc22a1",
#     "Cyp2e1",
#     "Cyp7a1",
#     "Gm30117",
#     "Aldh1a1",
#     "Cyp2c50",
#     "Lhpp",
#     "Cyp1a2",
#     "Slc1a2",
#     "Cyp2a5",
#     "Pon1",
#     "Slco1b2",
#     "Cyp2c54",
#     "Mgst1",
#     "Gulo",
#     "Lgr5",
# ]
# portal = [
#     "Ctsc",
#     "Aldh1b1",
#     "Pigr",
#     "Serpina1d",
#     "Apoc2",
#     "Hal",
#     "Gc",
#     "Hsd17b13",
#     "Itih4",
#     "Apoc3",
#     "Hpx",
#     "Hsd17b6",
#     "Apoa4",
#     "Sds",
#     "Serpina1b",
#     "Serpinc1",
#     "Alb",
#     "Cyp2f2",
#     "Pck1",
# ]

genes = central + portal

dev = "cpu"
clamp_gene = "Cyp2e1"  # gene to clamp

# training parameters
batch_size = 0  # batch size, put zero for full batch
n_iter = 10000  # number of iterations

# name of the output files
save = False
name = "big_set"
