#load packages
library(readr)
library(dplyr)
#for fetching fastas from uniprot
library(httr2) 
#as posed on uniprot API doc
library(jsonlite) 

#bait: protein on which Co-IP was performed
#prey: proteins which are bound directly or indirectly to bait

#Settings
  input_file <- "C:/Users/morit/OneDrive - JGU/Desktop/Arbeitsgruppe_Beli/R_Beli_Lab/Analyze_co-IP_hits_project/Test_input_file.csv"
  prey_protein_ids_coloumn_name <- "partner_uniprot"
  bait <- "P04637"
#do not add "/" here
  base_output_dir <- "Analyze_co-IP_hits" 
#create folder if it does not exist
  if (!dir.exists(base_output_dir)) {
    dir.create(base_output_dir, recursive = TRUE)
  }
  fasta_output_dir <- file.path(paste0(base_output_dir,"/fasta_output"))
#create folder if it does not exist
  if (!dir.exists(fasta_output_dir)) {
    dir.create(fasta_output_dir, recursive = TRUE)
  }
#create folder for AF2 multimer job
  job_name <- paste0("af2_multimer_",bait,"(bait)")
  job_dir  <- file.path(base_output_dir, job_name)
  if (!dir.exists(job_dir)) {
    dir.create(job_dir, recursive = TRUE)
  }
  
  
#read .csv input
  MS_data <- read_csv(input_file)
#get list of prey protein_ids
  prey_list <- unique(MS_data[[prey_protein_ids_coloumn_name]])

#get bait sequence from uniprot
#build the url
  url_bait <- paste0("https://rest.uniprot.org/uniprotkb/",bait,".fasta")
  params_bait <- list(fields = "sequence", version = "last")
#build the request
  req_bait <- request(url_bait)
#perform the request
  resp_bait <- req_perform(req_bait)
#check (if error) status
  if (resp_status(resp_bait) != 200) {
    stop(sprintf(
      "Error %d: %s",
      resp_status(resp_bait),
      resp_body_string(resp_bait)
    ))
  }
  
#get prey sequences from uniprot
  for (prey in prey_list) {
#build the url
    url <- paste0("https://rest.uniprot.org/uniprotkb/",prey, ".fasta")
    params <- list(fields = "sequence", version = "last")
#build the request
    req <- request(url)
#perform the request
    resp <- req_perform(req)
#check (if error) status
    if (resp_status(resp) != 200) {
      stop(sprintf(
        "Error %d: %s",
        resp_status(resp),
        resp_body_string(resp)
      ))
      
     
      
    }
#output: write one fasta file per prey
    fasta_txt <- resp_body_string(resp)
    bait_block <- resp_body_string(resp_bait)
    fasta_txt <- paste0(bait_block,fasta_txt)
    out_file <- file.path(fasta_output_dir, paste0(bait,"_",prey,".fasta"))
    writeLines(fasta_txt, out_file)  
  }
  
#create .sbatch file in job_dir ## "#SBATCH --array=0-16    
  sbatch_file <- file.path(job_dir, paste0(job_name, ".sbatch"))
  fasta_dir <- fasta_output_dir
  output_dir <- job_dir
  gpu_which <- "#SBATCH --partition=gpulong"
#adjusting array
    fasta_files <- list.files(fasta_dir)
    n_jobs <- length(fasta_files)
    adjusted_array <- paste0("#SBATCH --array=0-",n_jobs - 1)
  sbatch_lines <- c(
    "#!/usr/bin/env bash",
    paste0("#SBATCH --job-name=", job_name),
    "#SBATCH --output=logs/%x_%A_%a.out",
    "#SBATCH --error=logs/%x_%A_%a.err",
    adjusted_array,
    "#SBATCH --ntasks=1",
    "#SBATCH --cpus-per-task=8",
    "#SBATCH --gres=gpu:1",
    "#SBATCH --mem=48G",
    "#SBATCH --time=5:00:00",
    gpu_which,
    "set -euo pipefail",
    "",
    "######################## User settings ########################",
    paste0("FASTA_DIR=\"", fasta_dir, "\""),
    paste0("OUTPUT_DIR=\"", output_dir, "\""),
    "MODEL_PRESET=\"multimer\"",
    "GPU_RELAX=\"true\"",
    "###############################################################",
    "",
    "###################### Do not change unless necessary ######################",
    "mkdir -p \"$OUTPUT_DIR\" logs",
    "",
    "mapfile -d '' FASTAS < <(find \"$FASTA_DIR\" -maxdepth 1 -type f \\( -iname '*.fasta' -o -iname '*.fa' \\) -print0 | sort -z)",
    "",
    "if [[ ${SLURM_ARRAY_TASK_ID} -ge ${#FASTAS[@]} ]]; then",
    "  echo \"Array index ${SLURM_ARRAY_TASK_ID} out of range (found ${#FASTAS[@]} files).\"",
    "  exit 1",
    "fi",
    "",
    "# reverses the order of files drawn for computation (last file comes first..)",
    "#rev_idx=$(( ${#FASTAS[@]} - 1 - SLURM_ARRAY_TASK_ID ))",
    "#FASTA_PATH=\"${FASTAS[$rev_idx]}\"",
    "",
    "FASTA_PATH=\"${FASTAS[$SLURM_ARRAY_TASK_ID]}\" # uncomment this line and delete 2 lines above to get the default order",
    "BASENAME=\"$(basename \"$FASTA_PATH\")\"",
    "SAMPLE=\"${BASENAME%.*}\"",
    "",
    "# Per-job output directory",
    "JOB_OUT=\"$OUTPUT_DIR/$SAMPLE\"",
    "mkdir -p \"$JOB_OUT\"",
    "",
    "echo \"Task ${SLURM_ARRAY_TASK_ID}: FASTA=${FASTA_PATH}\"",
    "echo \"Output dir: ${JOB_OUT}\"",
    "",
    "# Environment (AlphaFold recommendations)",
    "export TF_FORCE_UNIFIED_MEMORY=1",
    "export XLA_PYTHON_CLIENT_MEM_FRACTION=4.0",
    "",
    "# Dates: keep max_template_date at or before today (here: 2025-10-29)",
    "if [[ $MODEL_PRESET == \"multimer\" ]]; then",
    "  singularity run --nv --nvccli --writable-tmpfs --bind /home,/fsimb,/media,/mnt,/tmp \\",
    "    /mnt/storage/alphafold/v232/alphafold_2.3.2.sif \\",
    "      --use_gpu_relax=\"$GPU_RELAX\" \\",
    "      --model_preset=\"$MODEL_PRESET\" \\",
    "      --max_template_date=2025-10-29 \\",
    "      --data_dir=/mnt/storage/alphafold/v232 \\",
    "      --bfd_database_path=/mnt/storage/alphafold/v232/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \\",
    "      --mgnify_database_path=/mnt/storage/alphafold/v232/mgnify/mgy_clusters_2022_05.fa \\",
    "      --obsolete_pdbs_path=/mnt/storage/alphafold/v232/pdb_mmcif/obsolete.dat \\",
    "      --pdb_seqres_database_path=/mnt/storage/alphafold/v232/pdb_seqres/pdb_seqres.txt \\",
    "      --template_mmcif_dir=/mnt/storage/alphafold/v232/pdb_mmcif/mmcif_files \\",
    "      --uniprot_database_path=/mnt/storage/alphafold/v232/uniprot/uniprot.fasta \\",
    "      --uniref30_database_path=/mnt/storage/alphafold/v232/uniref30/UniRef30_2021_03 \\",
    "      --uniref90_database_path=/mnt/storage/alphafold/v232/uniref90/uniref90.fasta \\",
    "      --output_dir=\"$JOB_OUT\" \\",
    "      --fasta_paths=\"$FASTA_PATH\"",
    "else",
    "  singularity run --nv --nvccli --writable-tmpfs --bind /home,/fsimb,/media,/mnt,/tmp \\",
    "    /mnt/storage/alphafold/v232/alphafold_2.3.2.sif \\",
    "      --use_gpu_relax=\"$GPU_RELAX\" \\",
    "      --model_preset=\"$MODEL_PRESET\" \\",
    "      --max_template_date=2020-05-14 \\",
    "      --data_dir=/mnt/storage/alphafold/v232 \\",
    "      --bfd_database_path=/mnt/storage/alphafold/v232/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \\",
    "      --mgnify_database_path=/mnt/storage/alphafold/v232/mgnify/mgy_clusters_2022_05.fa \\",
    "      --obsolete_pdbs_path=/mnt/storage/alphafold/v232/pdb_mmcif/obsolete.dat \\",
    "      --pdb70_database_path=/mnt/storage/alphafold/v232/pdb70/pdb70 \\",
    "      --template_mmcif_dir=/mnt/storage/alphafold/v232/pdb_mmcif/mmcif_files \\",
    "      --uniref30_database_path=/mnt/storage/alphafold/v232/uniref30/UniRef30_2021_03 \\",
    "      --uniref90_database_path=/mnt/storage/alphafold/v232/uniref90/uniref90.fasta \\",
    "      --output_dir=\"$JOB_OUT\" \\",
    "      --fasta_paths=\"$FASTA_PATH\"",
    "fi"
  )
  
  writeLines(sbatch_lines, sbatch_file)
  
#run the sbatch file (only works, if Rscript is executed on cluster (hpcgpu))
  submit_sbatch <- try(system2("sbatch"), sbatch_file, stdout = TRUE, stderr = TRUE)
  print("Job submission result: \n","submit_result")