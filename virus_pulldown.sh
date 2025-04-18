#! /bin/bash
#SBATCH --output=logs/%j_virus_pulldown.out
#SBATCH --error=logs/%j_virus_pulldown.err
#SBATCH --partition=cyaa #change partition to your own GPU partition
#SBATCH --gres=gpu:1
#SBATCH --mem=40G
#SBATCH --cpus-per-task=8


#Load python
module load Python/3.10.7

#Get the input arguments
virus_pdb=$1
human_receptor_csv=$2
output=$3
process_virus=$4
#Copy the input files to the working directory
mkdir -p $output

# Check that the environment exists and activate it 
VENV_DIR="ted_consensus_1.0/ted_consensus"
source "$VENV_DIR/bin/activate"

########################################################
##### PROCESS THE VIRUS PROTEIN ########################
########################################################
if [ "$process_virus" = true ]; then
    mkdir -p "$output/virus/virus_pdb"
    cp "$virus_pdb" "$output/virus/virus_pdb/"
    # Get the base name of the virus pdb
    virus_pdb_name=$(basename "$virus_pdb")
    virus_pdb="$output/virus/virus_pdb/$virus_pdb_name"

    VIRUS_PDB_DIR="$output/virus/virus_pdb"

    # Run segmentation on virus protein
    bash ted_consensus_1.0/run_segmentation.sh -i "$VIRUS_PDB_DIR" -o "$output/virus/virus_pdb_segmentation" 

    # Prepare the input for colabfold
    python srcs/virus_pulldown.py "$output/virus" "$virus_pdb"

    file="$output/virus/virus_domains_colabfold.csv"
    out_dir="$output/virus/virus_domains_colabfold"

    # module load MMseqs2/15-6f452 Kalign/2.04 cuda/11.6 cudnn/11.x-v8.7.0.84 ColabFold/1.5.3
    # module load blast/2.2.26 dssp/2.2.1 psipred/4.02 gcc/9.2.0 openmpi/4.0.5 hhsuite/3.3.0

    # start_time=$(date +%s)
    # colabfold_batch "$file" "$out_dir" \
    #     --use-gpu-relax \
    #     --amber \
    #     --templates \
    #     --num-relax 1 \
    #     --num-models 5 \
    #     --num-recycle 6
    # end_time=$(date +%s)
    # duration=$((end_time - start_time))
    # echo "Completed processing of $virus_pdb_name in $duration seconds at $(date)"
fi

##########################################################
##### PROCESS THE HUMAN RECEPTOR PROTEINS ################
##########################################################
mkdir -p "$output/human_receptor"
cp "$human_receptor_csv" "$output/human_receptor/"

# Get the base name of the human receptor csv
human_receptor_csv_name=$(basename "$human_receptor_csv")
human_receptor_csv="$output/human_receptor/$human_receptor_csv_name"

# Get the human receptor pdbs by fetching the AFDB entries
python srcs/human_receptor.py "$human_receptor_csv" "$output/human_receptor/human_receptor_pdbs"

# Run segmentation on human receptor protein
bash ted_consensus_1.0/run_segmentation.sh -i "$output/human_receptor/human_receptor_pdbs" -o "$output/human_receptor/human_receptor_segmentation"

# Prepare the input for colabfold
python srcs/human_segmentation.py "$output/human_receptor" "$output/human_receptor/human_receptor_pdbs"

files="$output/human_receptor/human_receptor_domains_colabfold/inputs/*.csv"
out_dir="$output/human_receptor/human_receptor_domains_colabfold/outputs"
if [ ! -d "$out_dir" ]; then
    mkdir -p "$out_dir"
fi


# for file in $files; do
#     module load MMseqs2/15-6f452 Kalign/2.04 cuda/11.6 cudnn/11.x-v8.7.0.84 ColabFold/1.5.3
#     module load blast/2.2.26 dssp/2.2.1 psipred/4.02 gcc/9.2.0 openmpi/4.0.5 hhsuite/3.3.0
#     start_time=$(date +%s)
#     colabfold_batch "$file" "$out_dir" \
#         --use-gpu-relax \
#         --amber \
#         --templates \
#         --num-relax 1 \
#         --num-models 5 \
#         --num-recycle 6
#     end_time=$(date +%s)
#     duration=$((end_time - start_time))
#     echo "Completed processing of $file in $duration seconds at $(date)"
# done


###Prepare the input for colabfold complex

# Get the virus domains
virus_input_file="$output/virus/virus_domains_colabfold.csv"

# Get the human receptor domains
human_receptor_input_dir="$output/human_receptor/human_receptor_domains_colabfold/inputs"

# Prepare the input for colabfold complex
python srcs/prepare_colabfold_complex.py "$virus_input_file" "$human_receptor_input_dir" "$output/colabfold_complex/inputs"
input_file="$output/colabfold_complex/inputs/complex_input.csv"
output_dir="$output/colabfold_complex/outputs"

module load MMseqs2/15-6f452 Kalign/2.04 cuda/11.6 cudnn/11.x-v8.7.0.84 ColabFold/1.5.3
module load blast/2.2.26 dssp/2.2.1 psipred/4.02 gcc/9.2.0 openmpi/4.0.5 hhsuite/3.3.0

# start_time=$(date +%s)
# colabfold_batch "$input_file" "$output_dir" \
#     --use-gpu-relax \
#     --amber \
#     --templates \
#     --num-relax 1 \
#     --num-models 5 \
#     --num-recycle 6
# end_time=$(date +%s)
# duration=$((end_time - start_time))
# echo "Completed processing of $input_file in $duration seconds at $(date)"
