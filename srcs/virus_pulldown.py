# Parse the segmentation results

from utils import parse_segmentation_output, get_seq_from_pdb
import pandas as pd
import sys
import os
OUTDIR = sys.argv[1]
VIRUS_PDB = sys.argv[2]
# For chainsaw
df_chainsaw = parse_segmentation_output(
    segmentation_file=f'{OUTDIR}/virus_pdb_segmentation/chopping_chainsaw.txt',
)

# For merizo
df_merizo = parse_segmentation_output(
    segmentation_file=f'{OUTDIR}/virus_pdb_segmentation/chopping_merizo.txt'
)

#merge the two dataframes and remove duplicates
df = pd.concat([df_chainsaw, df_merizo], ignore_index=True)
df = df.drop_duplicates(subset='domain_id')

#get the sequence of the virus protein
virus_seq = get_seq_from_pdb(pdb_file=f'{VIRUS_PDB}')

#prepare the input for colabfold
with open(f'{OUTDIR}/virus_domains_colabfold.csv', 'w') as f:
    f.write('id,sequence\n')
    for _, row in df.iterrows():
        # 1-based indexing
        f.write(f"{row['domain_id']},{virus_seq[int(row['start'])-1:int(row['end'])]}\n")
    #write a full length sequence
    #get the base name of the virus pdb
    virus_pdb_name = os.path.basename(VIRUS_PDB)
    virus_pdb_name = virus_pdb_name.split('.')[0]
    f.write(f"{virus_pdb_name}_FL,{virus_seq}\n")
