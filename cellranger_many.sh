#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --time=0-12:00:00
#SBATCH --output=01-Cellranger_Count_log.slurm.%a.stdout
#SBATCH --array=1-10
#SBATCH --mail-user=EMAIL
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="CRC"

# declare the name of the array
declare -A DATA

# each DATA below corresponds to the first part of the FASTQ file name for example "21-01877_S1_L001_R2_001.fastq.gz"
# cellranger requires Illumina format FASTQ file names
# paired-end sequencing

DATA[1]='21-01749'
DATA[2]='21-01750'
DATA[3]='21-01751'
DATA[4]='21-01756'
DATA[5]='21-01757'
DATA[6]='21-01758'
DATA[7]='21-01872'
DATA[8]='21-01873'
DATA[9]='21-01874'
DATA[10]='21-01875'
DATA[1]='21-01876'
DATA[2]='21-01877'

# add commands as needed
cellranger count --id=${DATA[$SLURM_ARRAY_TASK_ID]} --fastqs=./ --sample=${DATA[$SLURM_ARRAY_TASK_ID]} --transcriptome=/PATH/refdata-gex-mm10-2020-A --no-bam
