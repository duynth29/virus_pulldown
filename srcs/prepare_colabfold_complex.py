import os, sys
import pandas as pd
import glob

#read the input files
virus_input_file = sys.argv[1]
human_receptor_input_dir = sys.argv[2]
output_dir = sys.argv[3]
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

#read the virus input file
virus_df = pd.read_csv(virus_input_file)

human_receptor_files = glob.glob(os.path.join(human_receptor_input_dir, '*.csv'))


with open(os.path.join(output_dir, 'complex_input.csv'), 'w') as f:
    f.write('id,sequence\n')
    for file in human_receptor_files:
        df = pd.read_csv(file)
        for _, row in df.iterrows():
            for _, virus_row in virus_df.iterrows():
                f.write(f"{virus_row['id']}_{row['id']},{virus_row['sequence']}:{row['sequence']}\n")












