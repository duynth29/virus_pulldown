# Parse the segmentation results

from utils import parse_segmentation_output, get_seq_from_pdb
import pandas as pd
import sys
import os
OUTDIR = sys.argv[1]
HUMAN_RECEPTOR_PDBS_DIR = sys.argv[2]
# For chainsaw
df_chainsaw = parse_segmentation_output(
    segmentation_file=f'{OUTDIR}/human_receptor_segmentation/chopping_chainsaw.txt',
)

# For merizo
df_merizo = parse_segmentation_output(
    segmentation_file=f'{OUTDIR}/human_receptor_segmentation/chopping_merizo.txt'
)

#merge the two dataframes and remove duplicates
df = pd.concat([df_chainsaw, df_merizo], ignore_index=True)
df = df.drop_duplicates(subset='domain_id')

#group df by protein
df_by_protein = df.groupby('protein')
outdir = f'{OUTDIR}/human_receptor_domains_colabfold'
os.makedirs(outdir, exist_ok=True)
os.makedirs(f'{outdir}/inputs', exist_ok=True)
for protein, df_protein in df_by_protein:
    #get the sequence of the protein
    protein_seq = get_seq_from_pdb(pdb_file=f'{HUMAN_RECEPTOR_PDBS_DIR}/{protein}.pdb')
    #prepare the input for colabfold
    with open(f'{outdir}/inputs/{protein}_domains_colabfold.csv', 'w') as f:
        f.write('id,sequence\n')
        for _, row in df_protein.iterrows():
            # 1-based indexing for domain sequence extraction
            try:
                f.write(f"{row['domain_id']},{protein_seq[int(row['start'])-1:int(row['end'])]}\n")
            except IndexError as e:
                print(f"Error processing domain {row['domain_id']}: {e}")
        #write a full length sequence
        #get the base name of the virus pdb
        protein_name = os.path.basename(protein)
        protein_name = protein_name.split('.')[0]
        f.write(f"{protein_name}_FL,{protein_seq}\n")
