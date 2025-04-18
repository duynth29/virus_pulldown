from utils import get_pdb_AFDB
import sys, os
import pandas as pd
human_receptor_csv = sys.argv[1]
output_dir = sys.argv[2]

human_receptor_df = pd.read_csv(human_receptor_csv)

#create the output directory
os.makedirs(output_dir, exist_ok=True)

for _, row in human_receptor_df.iterrows():
    get_pdb_AFDB(row['uniprot_id'], row['name'], output_dir)




