import warnings
warnings.simplefilter(action='ignore')
import pycisTopic
pycisTopic.__version__
import pandas as pd
import scanpy as sc
import scipy
import pickle
import numpy as np
from scipy import sparse
import os

# Project directory
projDir = '/home3/ciervo/scMULTIOME/Analisi/'

# Output directory
outDir = projDir + 'scRNA_scATAC/SCENICplus/output/'
import os
if not os.path.exists(outDir):
    os.makedirs(outDir)

# Temp dir
tmpDir = projDir + 'scRNA_scATAC/SCENICplus/tmp/'
if not os.path.exists(tmpDir):
    os.makedirs(tmpDir)

# __Creating a cisTopic object__
# Create cisTopic object
#from pycisTopic.cistopic_class import *
path_to_blacklist='/home3/ciervo/scMULTIOME/Analisi/scRNA_scATAC/SCENICplus/pycisTopic/blacklist/hg38-blacklist.v2.bed'

# load count matrix
count_matrix= scipy.io.mmread('/home3/ciervo/scMULTIOME/Analisi/scRNA_scATAC/SCENICplus/malignant_counts.mtx').todense()
print(count_matrix)

# load cell annotation
cell_annotation=pd.read_csv(projDir+'scRNA_scATAC/SCENICplus/cell_annotation.tsv', sep='\t')
# load peak annotation
peak_annotation=pd.read_csv(projDir+'scRNA_scATAC/SCENICplus/peak_annotation.tsv', sep='\t')

cell_annotation = cell_annotation['x'].to_list()
peak_annotation = peak_annotation['x'].to_list()

count_matrix = np.asarray(count_matrix)
count_df = pd.DataFrame(count_matrix, index=peak_annotation, columns=cell_annotation)
count_matrix_sparse = sparse.csr_matrix(count_df.values)

from pycisTopic.cistopic_class import *
cistopic_obj = create_cistopic_object(fragment_matrix=count_matrix_sparse, path_to_blacklist=path_to_blacklist, cell_names=cell_annotation, region_names=peak_annotation)

# Adding cell information
cell_data=pd.read_csv('/home3/ciervo/scMULTIOME/Analisi/scRNA_scATAC/SCENICplus/metadata.tsv', sep='\t')
cistopic_obj.add_cell_data(cell_data)

#print(cistopic_obj)

pickle.dump(cistopic_obj, 
            open(os.path.join(outDir, 'cisTopic_malignant.pkl'), 'wb'))

cistopic_obj = pickle.load(open(os.path.join(outDir, 'cisTopic_malignant.pkl'), 'rb'))

path_to_mallet_binary='/home3/ciervo/scMULTIOME/Analisi/scRNA_scATAC/SCENICplus/pycisTopic/Mallet-202108/bin/mallet'
os.environ['MALLET_MEMORY'] = '200G'

from pycisTopic.lda_models import run_cgs_models_mallet

# Run models
models=run_cgs_models_mallet(
  mallet_path = path_to_mallet_binary,
  cistopic_obj = cistopic_obj,
  n_topics=[2,5,10,15,20,25,30,35,40,45,50],
  n_cpu=100,
  n_iter=500,
  random_state=555,
  alpha=50,
  alpha_by_topic=True,
  eta=0.1,
  eta_by_topic=False,
  tmp_path= tmpDir, 
  save_path = None)

# Save
import pickle 
with open(outDir+'Mallet_models.pkl', 'wb') as f:
  pickle.dump(models, f)
  
os.mkdir(outDir+'models')
from pycisTopic.lda_models import evaluate_models
model=evaluate_models(models,
                     select_model=None, 
                     return_model=True, 
                     metrics=['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
                     plot_metrics=False,
                     save= outDir + 'models/model_selection.pdf')

# Add model to cisTopicObject
cistopic_obj.add_LDA_model(model)

# Save
with open(outDir + 'cisTopic_malignant.pkl', 'wb') as f:
  pickle.dump(cistopic_obj, f)

# Inferring candidate enhancer regions
from pycisTopic.topic_binarization import *
region_bin_topics_otsu = binarize_topics(cistopic_obj, method='otsu')
region_bin_topics_top3k = binarize_topics(cistopic_obj, method='ntop', ntop = 3000)

from pycisTopic.diff_features import *
imputed_acc_obj = impute_accessibility(cistopic_obj, selected_cells=None, selected_regions=None, scale_factor=10**6)
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)
variable_regions = find_highly_variable_features(normalized_imputed_acc_obj, plot = False)
markers_dict = find_diff_features(cistopic_obj, imputed_acc_obj, variable = 'Metaprogram_assignment', var_features=variable_regions, split_pattern = '-')

os.mkdir(outDir+'candidate_enhancers')

import pickle
pickle.dump(region_bin_topics_otsu, open(outDir + '/candidate_enhancers/region_bin_topics_otsu.pkl', 'wb'))
pickle.dump(region_bin_topics_top3k, open(outDir+'/candidate_enhancers/region_bin_topics_top3k.pkl', 'wb'))
pickle.dump(markers_dict, open(outDir+'/candidate_enhancers/markers_dict.pkl', 'wb'))

import os
import pandas as pd
import pickle

# Create output directories
os.makedirs(outDir+'region_sets/Topics_otsu', exist_ok=True)
os.makedirs(outDir+'region_sets/Topics_top_3k', exist_ok=True)
os.makedirs(outDir+'region_sets/DARs_cell_type', exist_ok=True)

# Function to extract BED fields from row names like "chr6:14000921-14002088"
def index_to_bed(df_index):
    chrom = df_index.str.extract(r'(chr[\w]+):')[0]
    start = df_index.str.extract(r':(\d+)-')[0].astype(int)
    end = df_index.str.extract(r'-(\d+)$')[0].astype(int)
    return pd.DataFrame({'chrom': chrom, 'start': start, 'end': end})

# Function to save BED file from DataFrame with genomic index
def save_bed_from_dict(data_dict, out_dir):
    for key, df in data_dict.items():
        bed = index_to_bed(df.index.to_series())
        filename = key.replace(' ', '_').replace('/', '_vs_') + '.bed'
        bed.to_csv(os.path.join(out_dir, filename), sep='\t', header=False, index=False)

# Load and process each file
with open(outDir+'candidate_enhancers/region_bin_topics_otsu.pkl', 'rb') as f:
    topics = pickle.load(f)
    save_bed_from_dict(topics, outDir+'region_sets/Topics_otsu')

with open(outDir+'candidate_enhancers/region_bin_topics_top3k.pkl', 'rb') as f:
    top3k = pickle.load(f)
    save_bed_from_dict(top3k, outDir+'region_sets/Topics_top_3k')

with open(outDir+'candidate_enhancers/markers_dict.pkl', 'rb') as f:
    dars = pickle.load(f)
    save_bed_from_dict(dars, outDir+'region_sets/DARs_cell_type')


### Check names
with open(outDir + 'cisTopic_malignant.pkl', 'rb') as f:
  cistopic_obj = pickle.load(f)
