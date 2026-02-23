import os
import pickle
import gzip
import argparse
import mudata
import warnings
from scenicplus.simulation import simulate_perturbation
def load_object(filename):
    with open(filename, 'rb') as f:
        return pickle.load(f)

def save_object(obj, filename):
    with gzip.open(filename + '.gz', 'wb') as f:
        pickle.dump(obj, f)
        
warnings.simplefilter(action='ignore', category=FutureWarning)

# Transcription factor
parser = argparse.ArgumentParser()
parser.add_argument('--tf', required=True)
args = parser.parse_args()
tf_name = args.tf

# Path
data_dir = "/home3/ciervo/scMULTIOME/Analisi/scRNA_scATAC/SCENICplus/"

output_dir = os.path.join(data_dir, 'Perturbation/MP_4/')
os.makedirs(output_dir, exist_ok=True)

# Loading data
print(f"Loading data: {tf_name}")
scplus_mdata = mudata.read(os.path.join(data_dir, "SCENICplus_output_2/scplusmdata.h5mu"))
expression_df = scplus_mdata["scRNA_counts"].to_df()
regressors = load_object(os.path.join(data_dir, 'Perturbation/regressors.pkl'))
del scplus_mdata

# Perturbation
if tf_name not in regressors:
    print(f"[SKIP] {tf_name} not in regressors")
else:
    try:
        print(f"Starting TF perturbation: {tf_name}")
        result = simulate_perturbation(
            df_EXP=expression_df,
            perturbation={tf_name: 0},
            keep_intermediate=True,
            n_iter=5,
            regressors=regressors
        )
        save_path = os.path.join(output_dir, f"{tf_name}_perturbation.pkl")
        save_object(result, save_path)
        print(f"Saved: {save_path}")
    except Exception as e:
        print(f"[ERROR] {tf_name}: {e}")
