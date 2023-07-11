path = "data/270323_data_hepatocytes.h5ad"

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
genes = central + portal


dev = "cpu"
clamp_gene = "Cyp2e1"  # gene to clamp


# training parameters
batch_size = 128  # batch size, put zero for full batch
n_iter = 3  # number of iterations

# name of the output files
name = "first_run"
