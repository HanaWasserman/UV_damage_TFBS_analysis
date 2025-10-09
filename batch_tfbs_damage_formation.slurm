#!/bin/bash
#SBATCH --job-name=prep_tfbs
#SBATCH --output=logs/output_%A_%a.out
#SBATCH --error=logs/error_%A_%a.err
#SBATCH --mem=8G
#SBATCH --array=0-196
#SBATCH --cpus-per-task=16
#SBATCH --partition=compsci

source /home/users/hiw4/miniconda3/etc/profile.d/conda.sh
conda activate tfbs_cpd

set -e
set -o pipefail

CONFIG_FILE="$1"

source "$CONFIG_FILE"

REQUIRED_VARS=("MAIN_DIR" "MOTIF_CLUSTER_FILE")
for var in "${REQUIRED_VARS[@]}"; do
  if [ -z "${!var}" ]; then
    echo "Error: Variable '$var' not set in config file."
    exit 1
  fi
done


line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${MAIN_DIR}/data/${MOTIF_CLUSTER_FILE})
IFS=',' read -r TF_CLUSTER TF_LEN TF <<< "$line"

echo "Analyzing damage formation at $TF sites..."
bash run_analyze_TFBS_damage_formation.bs "$CONFIG_FILE" "$TF_LEN" "$TF"

echo "Finished for $TF."
