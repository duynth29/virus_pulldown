import os, glob, json
import requests
import pandas as pd
from Bio import PDB
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB import PDBParser, Superimposer
# from Bio.PDB.Polypeptide import three_to_one 
from Bio.Align import PairwiseAligner
def get_pdb_AFDB(id, name, out_dir):
    """Fetch the PDB file from the AFDB and save it to the specified directory.
    
    Args:
        id (str): Uniprot ID
        out_dir (str): The directory to save the PDB file
    """
    url = f"https://alphafold.ebi.ac.uk/files/AF-{id}-F1-model_v4.pdb"
    response = requests.get(url)
    
    if response.status_code == 200:
        # Ensure the output directory exists
        os.makedirs(out_dir, exist_ok=True)
        
        # Construct the full file path
        file_path = os.path.join(out_dir, f"{name}.pdb")
        
        # Write the PDB file to the specified directory
        with open(file_path, "w") as f:
            f.write(response.text)
        print(f"Successfully fetched and saved PDB file for {name} to {file_path}")
    else:
        print(f"Failed to fetch PDB file for {name}. HTTP Status Code: {response.status_code}")



def get_seq_from_pdb(pdb_file: str, chain: str = None) -> None:
    """
    Convert a PDB file to a sequence (one-letter code).
    """
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_file)
    # Use the first model and the first chain
    model = structure[0] #only use the first model
    if chain is not None:
        chain = model[chain]
    else:
        chain = list(model.get_chains())[0] #only use the first chain if no chain is specified
    ppb = PPBuilder()
    seq = ""
    for pp in ppb.build_peptides(chain):
        seq += str(pp.get_sequence())
    return seq


def parse_segmentation_output(segmentation_file: str, output_dir: str = None) -> pd.DataFrame:
    """
    Parse the segmentation output from the chainsaw/merizo tool and return a pandas DataFrame.
    Example line:
    AF-A0A6C0LIE9-F1-model_v4	53ea6ddca5bc224a1e351b4739d7087e	304	1	3-302	0.0241
    AF-A0A1V6M2Y0-F1-model_v4	17c1e5cdd97ce14fb41a72ac2f0a15a8	1272	4	101-134,201-314_351-419,425-1125,1178-1269	0.0443
    
    Returns a DataFrame with the following columns:
    - domain_id: a unique identifier for each domain
    - start: start position of the domain
    - end: end position of the domain
    """
    raw_df = pd.read_csv(segmentation_file, sep="\t", header=None, usecols=[0, 4], names=["protein_id", "list_of_domains"])
    
    domain_records = []
    for _, row in raw_df.iterrows():
        domains = row["list_of_domains"].replace("_", "-").split(",")
        for domain in domains:
            start = domain.split("-")[0]
            end = domain.split("-")[-1]
            domain_id = f"{row['protein_id']}_{start}_{end}"
            domain_records.append({"protein": row['protein_id'], "domain_id": domain_id, "start": start, "end": end})
    
    domain_df = pd.DataFrame(domain_records)
    base_name = os.path.basename(segmentation_file).split(".")[0]
    
    # Determine output directory
    output_dir = output_dir or os.path.dirname(segmentation_file)
    domain_df.to_csv(os.path.join(output_dir, f"{base_name}_domains.csv"), index=False)
    
    return domain_df


def align_and_extract_domain_RMSD(full_pdb_path: str, domain_pdb_path: str) -> float:
    """
    Aligns the domain structure with a segment in the full protein structure and returns the minimal RMSD value.

    Assumptions:
     - First perform a pairwise sequence alignment (local) between the domain sequence and the full protein chain's sequence.
     - Then, using the alignment, extract the corresponding CA atoms that form a contiguous segment.

    Parameters:
      full_pdb_path (str): Path to the full protein PDB file.
      domain_pdb_path (str): Path to the domain PDB file.

    Returns:
      float: Minimal RMSD achieved when aligning the domain structure with a segment in the full protein.

    Raises:
      ValueError: If no CA atoms are found in the domain or if no matching segment is found in the full protein.
    """


    # Parse PDB files
    parser = PDBParser(QUIET=True)
    full_structure = parser.get_structure("full", full_pdb_path)
    domain_structure = parser.get_structure("domain", domain_pdb_path)

    # Use the first chain of the domain structure for consistency
    domain_model = domain_structure[0]
    domain_chain = list(domain_model.get_chains())[0]
    domain_atoms = [res['CA'] for res in domain_chain if 'CA' in res]
    if not domain_atoms:
        raise ValueError("No CA atoms found in the domain structure.")

    # Create one-letter domain sequence using three_to_one conversion so each residue is one character.
    domain_sequence = ''.join(three_to_one(res.get_resname(), warn=False) for res in domain_chain if 'CA' in res)
    num_domain_atoms = len(domain_atoms)
    best_rmsd = None

    # Align with the first model of the full protein
    full_model = full_structure[0]

    # Iterate over each chain in the full protein
    for chain in full_model.get_chains():
        chain_atoms = [res for res in chain if 'CA' in res]
        # Build one-letter code sequence for the chain
        chain_sequence = ''.join(three_to_one(res.get_resname(), warn=False) for res in chain_atoms)
        if len(chain_atoms) < num_domain_atoms:
            continue  # Skip chains with insufficient residues

        # Perform local alignment between domain_sequence and chain_sequence
        aligner = PairwiseAligner()
        aligner.mode = 'local'
        alignments = aligner.align(domain_sequence, chain_sequence)
        if not alignments:
            continue
        best_alignment = max(alignments, key=lambda x: x[2])
        aligned_domain, aligned_chain, score, _, _ = best_alignment

        # Skip alignment if the chain part contains gaps corresponding to domain residues.
        # This ensures we extract a contiguous segment from chain_atoms.
        if any(c == '-' for d, c in zip(aligned_domain, aligned_chain) if d != '-'):
            continue

        # Map the aligned region to indices in chain_atoms.
        candidate_start = None
        candidate_end = None
        chain_index = 0
        for i in range(len(aligned_chain)):
            if aligned_chain[i] != '-':
                if candidate_start is None and aligned_domain[i] != '-':
                    candidate_start = chain_index
                if aligned_domain[i] != '-':
                    candidate_end = chain_index + 1
                chain_index += 1
        if candidate_start is None or candidate_end is None:
            continue
        candidate_atoms = chain_atoms[candidate_start:candidate_end]
        if len(candidate_atoms) != num_domain_atoms:
            continue

        sup = Superimposer()
        try:
            sup.set_atoms(candidate_atoms, domain_atoms)
        except Exception:
            continue  # Skip alignment errors
        current_rmsd = sup.rms
        if best_rmsd is None or current_rmsd < best_rmsd:
            best_rmsd = current_rmsd

    if best_rmsd is None:
        raise ValueError("No matching domain segment found in the full protein structure.")

    return best_rmsd

def color_domain_in_chimera(csv_file: str, output_dir: str, chain: str = "*") -> None:
    """
    Color the domains in the PDB file based on the segmentation results.

    Args:
        csv_file (str): Path to the CSV file containing domain information.
        output_dir (str): Directory to save the output PML file.
    """
    df = pd.read_csv(csv_file)
    base_name = os.path.basename(csv_file).split(".")[0]
    colors = [
        "#FFB3BA", "#FFDFBA", "#FFFFBA", "#BAFFC9", "#BAE1FF",
        "#FFABAB", "#FFC3A0", "#FF677D", "#D4A5A5", "#392F5A"
    ]
    df = df.groupby("protein")
    for protein, df_protein in df:
        with open(f"{output_dir}/{protein}_{base_name}_chimera_color.cxc", "w") as f:
            for index, row in df_protein.iterrows():
                color_index = index % len(colors)
                color = colors[color_index]
                cmd = f"color #*/{chain}:{row['start']}-{row['end']} {color}"  # Define the color in Chimera
                f.write(cmd + "\n")



def get_scores_colabfold(colabfold_dir):
    """Get scores from colabfold output.
    
    Args:
        colabfold_dir (str): Path to colabfold output.
        
    Returns:
        pd.DataFrame: DataFrame containing columns:
            model (str): Model identifier.
            rank (int): Rank of the model.
            iptm (float): ipTM score (if it exists).
            ptm (float): pTM score.
            
    Example file output:
        OC43_1230_1276_TMPS2_377_491_scores_rank_001_alphafold2_multimer_v3_model_3_seed_000.json
    """
    json_files = glob.glob(os.path.join(colabfold_dir, "*scores*.json"))
    scores_list = []
    
    for json_file in json_files:
        with open(json_file, "r") as f:
            data = json.load(f)
        
        # If iptm exists, get it
        if "iptm" in data:
            iptm = data.get("iptm")
        else:
            iptm = None

        #get ptm
        ptm = data.get("ptm")
       
        
        # Get the file name without its extension for further parsing
        file_name = os.path.splitext(os.path.basename(json_file))[0]
        parts = file_name.split("_")
        
        # Extract rank and model using the keywords in the filename
        try:
            rank = parts[parts.index("rank") + 1]
        except (ValueError, IndexError):
            rank = None
        
        try:
            model = parts[parts.index("model") + 1]
        except (ValueError, IndexError):
            model = None
        
        complex_name = file_name.split("_scores_")[0]
        scores_list.append({
            "model": model,
            "rank": rank,
            "complex_name": complex_name,
            "iptm": iptm,
            "ptm": ptm,
        })
    
    return pd.DataFrame(scores_list)

def convert_colabfold_input_to_AF3_input(csv_file_path, job_limit=2, output_dir="output/AF3_input", overwrite=False):
    """
    Convert the ColabFold input CSV to AlphaFold3 JSON input format.

    Input CSV format:
      - Columns: "id", "sequence"
      - The "sequence" column should be in the format "SEQA:SEQB"

    Output:
      A list of dictionaries in the following format:
      [
          {
              "name": "Complex_1",
              "modelSeeds": [],
              "sequences": [
                  {
                      "proteinChain": {
                          "sequence": "SEQA",
                          "count": 1,
                          "useStructureTemplate": True
                      }
                  },
                  {
                      "proteinChain": {
                          "sequence": "SEQB",
                          "count": 1,
                          "useStructureTemplate": True
                      }
                  }
              ],
              "dialect": "alphafoldserver",
              "version": 1
          },
          {
              "name": "Complex_2",
              "modelSeeds": [],
              "sequences": [
                  {
                      "proteinChain": {
                          "sequence": "SEQA",
                          "count": 1,
                          "useStructureTemplate": True
                      }
                  },
                  {
                      "proteinChain": {
                          "sequence": "SEQB",
                          "count": 1,
                          "useStructureTemplate": True
                      }
                  }
              ],
              "dialect": "alphafoldserver",
              "version": 1
          }
      ]

    Parameters:
      csv_file_path (str): Path to the ColabFold input CSV file.
      job_limit (int): Maximum number of jobs (rows) to process (default: 2).

    Returns:
      List[dict]: A list of dictionaries in JSON format ready for AF3 input.
    """

    if not os.path.exists(output_dir):  # Create output directory if it doesn't exist
        os.makedirs(output_dir)
    # Check if the json file already exists
    if glob.glob(os.path.join(output_dir, "AF3_input_job_*.json")):
        print("AF3 input json file already exists")
        if overwrite:
            print("Deleting the existing json file")
            for file in glob.glob(os.path.join(output_dir, "AF3_input_job_*.json")):
                os.remove(file)
        else:
            print("Skipping the json file")
            return None
    df = pd.read_csv(csv_file_path)  # Read the input CSV file
    n_jobs = (df.shape[0] + job_limit - 1) // job_limit  # Calculate number of jobs

    for i in range(n_jobs):  # Iterate over jobs
        df_subset = df.iloc[i * job_limit: (i + 1) * job_limit]  # Get subset of DataFrame
        complexes = []  # Initialize list for complex entries
        for idx, row in df_subset.iterrows():  # Iterate over rows in subset
            identifier = row['id']  # Extract identifier
            sequence_field = row['sequence']  # Extract sequence field
            if ':' not in sequence_field:  # Validate sequence format
                raise ValueError(
                    f"Invalid sequence format in row {idx}: '{sequence_field}'. Expected format 'SEQA:SEQB'."
                )
            parts = sequence_field.split(':')  # Split sequences
            if len(parts) < 2:  # Ensure two parts are present
                raise ValueError(
                    f"Invalid sequence format in row {idx}: '{sequence_field}'. Expected two parts separated by ':'."
                )
            seq_a = parts[0].strip()  # Clean up sequence A
            seq_b = parts[1].strip()  # Clean up sequence B

            complex_entry = {  # Create complex entry dictionary
                "name": f"{identifier}",
                "modelSeeds": [],
                "sequences": [
                    {
                        "proteinChain": {
                            "sequence": seq_a,
                            "count": 1,
                            "useStructureTemplate": True
                        }
                    },
                    {
                        "proteinChain": {
                            "sequence": seq_b,
                            "count": 1,
                            "useStructureTemplate": True
                        }
                    }
                ],
                "dialect": "alphafoldserver",
                "version": 1
            }
            complexes.append(complex_entry)  # Append complex entry to list
        output_file_path = os.path.join(output_dir, f"AF3_input_job_{i}.json")  # Define output file path
        with open(output_file_path, "w") as f:  # Write complexes to JSON file
            json.dump(complexes, f, indent=4)  # Save JSON with indentation
    print(f"AF3 input json file saved to {output_dir}")