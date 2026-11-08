# alphafold-multimer-coip-screening
R scripts and a Slurm pipeline to prepare AlphaFold Multimer inputs (fragmentation approach) and rank Co-IP hits by ipTM.


This repository contains R scripts for preparing AlphaFold Multimer input files (FASTA) using a fragmentation-based approach prior to prediction. It also includes scripts for ranking prediction results based on ipTM values, as well as a prototype pipeline designed to perform simple AlphaFold Multimer predictions from a CSV file containing Co-IP data. The fragmentation approach is following the strategy described by Lee et al. (2024).

 

The following files are contained:

 

Manual:

Creates_fastas_from_sequences_for_AFM.R (Fragmentation approach yielding fasta sequences based on 2 proteins with distinct domain and non-domain regions)

RankingAlphaFoldServerOutput.R (Ranks the output of Alphafoldserver.com (AF3, Dec. 2025))

Ranks_IMB_HPC_GPU_AFM_Output_.R (Ranks AlphaFold Multimer output from an institutes HPC)

Script_summary_NPLOC4ganz.R (Summarizes AlphaFold Multimer with prior filtering by keywords)

Script_summary_out.R (Summarizes AlphaFold Multimer output)

 

Pipeline:

Analyze_co-IP_hits_prepare_fastas.R (Pipeline)

Analyze_co-IP_hits_with_2_args.R (Pipeline runnable from command line)

Test_input_file.csv (example input file for Co-IP Pipeline)

af2_multimer_P04637(bait).sbatch (example slurm submission script)

 

Pipeline Workflow:

Login to your HPC containing AlphaFold Multimer
Create a folder containing Analyze_co-IP_hits_prepare_fastas.R, Analyze_co-IP_hits_with_2_args.R and a csv file containig Co-IP hits results
in home directory: module load R
start the preparation script: Rscript /path/to/Analyze_co-IP_hits_prepare_fastas.R path/to/csv_Co-IP_hits bait_protein_ID output_dir
start sbatch job: sbatch af2_multimer_bait_protein_ID.sbatch (make sure to change the name according to your bait protein)
The pipeline prepares fastas, assigning the bait protein to all prey proteins 1:1. When the sbatch job is started, AlphaFold Multimere predictions run and thereafter a Rscript ranks the output by their ipTM score.
