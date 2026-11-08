#!/usr/bin/env bash
#SBATCH --job-name=af2_multimer_P04637(bait)
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --array=0-9
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1
#SBATCH --mem=48G
#SBATCH --time=5:00:00
#SBATCH --partition=gpulong
set -euo pipefail

#Settings
FASTA_DIR="Analyze_co-IP_hits/fasta_output"
OUTPUT_DIR="Analyze_co-IP_hits/af2_multimer_P04637(bait)"
MODEL_PRESET="multimer"
GPU_RELAX="true"

#Try not to change anything here
mkdir -p "$OUTPUT_DIR" logs

mapfile -d '' FASTAS < <(find "$FASTA_DIR" -maxdepth 1 -type f \( -iname '*.fasta' -o -iname '*.fa' \) -print0 | sort -z)

if [[ ${SLURM_ARRAY_TASK_ID} -ge ${#FASTAS[@]} ]]; then
  echo "Array index ${SLURM_ARRAY_TASK_ID} out of range (found ${#FASTAS[@]} files)."
  exit 1
fi

#reverses the order of files drawn for computation (last file comes first..)
#rev_idx=$(( ${#FASTAS[@]} - 1 - SLURM_ARRAY_TASK_ID ))
#FASTA_PATH="${FASTAS[$rev_idx]}"

FASTA_PATH="${FASTAS[$SLURM_ARRAY_TASK_ID]}" # uncomment this line and delete 2 lines above to get the default order
BASENAME="$(basename "$FASTA_PATH")"
SAMPLE="${BASENAME%.*}"

#Per-job output directory
JOB_OUT="$OUTPUT_DIR/$SAMPLE"
mkdir -p "$JOB_OUT"

echo "Task ${SLURM_ARRAY_TASK_ID}: FASTA=${FASTA_PATH}"
echo "Output dir: ${JOB_OUT}"

#Environment (AlphaFold recommendations)
export TF_FORCE_UNIFIED_MEMORY=1
export XLA_PYTHON_CLIENT_MEM_FRACTION=4.0

# Dates: keep max_template_date at or before today (here: 2025-10-29)
if [[ $MODEL_PRESET == "multimer" ]]; then
  singularity run --nv --nvccli --writable-tmpfs --bind /home,/fsimb,/media,/mnt,/tmp \
    /mnt/storage/alphafold/v232/alphafold_2.3.2.sif \
      --use_gpu_relax="$GPU_RELAX" \
      --model_preset="$MODEL_PRESET" \
      --max_template_date=2025-10-29 \
      --data_dir=/mnt/storage/alphafold/v232 \
      --bfd_database_path=/mnt/storage/alphafold/v232/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
      --mgnify_database_path=/mnt/storage/alphafold/v232/mgnify/mgy_clusters_2022_05.fa \
      --obsolete_pdbs_path=/mnt/storage/alphafold/v232/pdb_mmcif/obsolete.dat \
      --pdb_seqres_database_path=/mnt/storage/alphafold/v232/pdb_seqres/pdb_seqres.txt \
      --template_mmcif_dir=/mnt/storage/alphafold/v232/pdb_mmcif/mmcif_files \
      --uniprot_database_path=/mnt/storage/alphafold/v232/uniprot/uniprot.fasta \
      --uniref30_database_path=/mnt/storage/alphafold/v232/uniref30/UniRef30_2021_03 \
      --uniref90_database_path=/mnt/storage/alphafold/v232/uniref90/uniref90.fasta \
      --output_dir="$JOB_OUT" \
      --fasta_paths="$FASTA_PATH"
else
  singularity run --nv --nvccli --writable-tmpfs --bind /home,/fsimb,/media,/mnt,/tmp \
    /mnt/storage/alphafold/v232/alphafold_2.3.2.sif \
      --use_gpu_relax="$GPU_RELAX" \
      --model_preset="$MODEL_PRESET" \
      --max_template_date=2020-05-14 \
      --data_dir=/mnt/storage/alphafold/v232 \
      --bfd_database_path=/mnt/storage/alphafold/v232/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
      --mgnify_database_path=/mnt/storage/alphafold/v232/mgnify/mgy_clusters_2022_05.fa \
      --obsolete_pdbs_path=/mnt/storage/alphafold/v232/pdb_mmcif/obsolete.dat \
      --pdb70_database_path=/mnt/storage/alphafold/v232/pdb70/pdb70 \
      --template_mmcif_dir=/mnt/storage/alphafold/v232/pdb_mmcif/mmcif_files \
      --uniref30_database_path=/mnt/storage/alphafold/v232/uniref30/UniRef30_2021_03 \
      --uniref90_database_path=/mnt/storage/alphafold/v232/uniref90/uniref90.fasta \
      --output_dir="$JOB_OUT" \
      --fasta_paths="$FASTA_PATH"
fi
