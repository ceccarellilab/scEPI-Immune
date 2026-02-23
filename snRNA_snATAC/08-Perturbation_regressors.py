import mudata
import os
import scanpy as sc
import anndata
import matplotlib
import matplotlib.pyplot as plt
import adjustText
import numpy as np
import pandas as pd
import pickle

# Function to save an object
def save_object(obj, filename):
    with open(filename, 'wb') as f:
        pickle.dump(obj, f)
        
# Function to load an object
def load_object(filename):
    with open(filename, 'rb') as f:
        return pickle.load(f)

from scenicplus.simulation import (
    train_gene_expression_models,
    simulate_perturbation,
    plot_perturbation_effect_in_embedding
)

# Directories
data_dir = "/home3/ciervo/scMULTIOME/Analisi/scRNA_scATAC/SCENICplus/"
os.mkdir(data_dir+'Perturbation')

print('loading data')
scplus_mdata = mudata.read(os.path.join(data_dir, "SCENICplus_output_2/scplusmdata.h5mu"))

print('eRegulon_gene_AUC')
eRegulon_gene_AUC = anndata.concat(
    [scplus_mdata["direct_gene_based_AUC"], scplus_mdata["extended_gene_based_AUC"]],
    axis = 1,
)
eRegulon_gene_AUC.obs = scplus_mdata.obs
eRegulon_gene_AUC

sc.pp.pca(eRegulon_gene_AUC)
color_dict_line = {
    'MP_1' : "#D55E00", 
    'MP_2' : "#E69F00", 
    'MP_3' : "#F0E442", 
    'MP_4' : '#8BC34A', 
    'MP_5' : "#009E73", 
    'MP_6' : "#0072B2", 
    'MP_7' : "#6A3D9A"
    }

def plot_mm_line_pca(ax):
    texts = []
    # Plot PCA
    ax.scatter(
        eRegulon_gene_AUC.obsm["X_pca"][:, 0],
        eRegulon_gene_AUC.obsm["X_pca"][:, 1],
        color = [color_dict_line[line] for line in eRegulon_gene_AUC.obs["scRNA_counts:Metaprogram_assignment"]]
    )
    # Plot labels
    for line in set(eRegulon_gene_AUC.obs["scRNA_counts:Metaprogram_assignment"]):
        line_bc_idc = np.arange(len(eRegulon_gene_AUC.obs_names))[eRegulon_gene_AUC.obs["scRNA_counts:Metaprogram_assignment"] == line]
        avg_x, avg_y = eRegulon_gene_AUC.obsm["X_pca"][line_bc_idc, 0:2].mean(0)
        texts.append(
            ax.text(
                avg_x,
                avg_y,
                line,
                fontweight = "bold"
            )
        )
    adjustText.adjust_text(texts)

print('plot pca')
fig, ax = plt.subplots()
plot_mm_line_pca(ax)
fig.savefig(data_dir+'Perturbation/pca_mps.pdf')   # save the figure to file
plt.close(fig)   
save_object(eRegulon_gene_AUC, data_dir+'Perturbation/eRegulon_gene_AUC.pkl')

gene_tf_direct_extended = pd.concat(
    [
        scplus_mdata.uns["direct_e_regulon_metadata"][["Gene", "TF"]].drop_duplicates(),
        scplus_mdata.uns["extended_e_regulon_metadata"][["Gene", "TF"]].drop_duplicates()
    ]
).drop_duplicates()
gene_to_TF = gene_tf_direct_extended.groupby("Gene")["TF"].apply(lambda tfs: list(tfs)).to_dict()

print('start regressors')
regressors = train_gene_expression_models(
    df_EXP = scplus_mdata["scRNA_counts"].to_df(),
    gene_to_TF = gene_to_TF
)

save_object(regressors, data_dir+'Perturbation/regressors.pkl')
