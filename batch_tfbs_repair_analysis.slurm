#!/bin/bash
#SBATCH --job-name=analyze_repair
#SBATCH --output=results/logs/output_%A_%a.out
#SBATCH --error=results/logs/error_%A_%a.err
#SBATCH --mem=8G
#SBATCH --array=0-37
#SBATCH --cpus-per-task=16
#SBATCH --partition=compsci

source /home/users/hiw4/miniconda3/etc/profile.d/conda.sh
conda activate tfbs_cpd

set -e
set -o pipefail

source /home/users/hiw4/miniconda3/etc/profile.d/conda.sh
conda activate tfbs_cpd

CONFIG_FILE="$1"

source "$CONFIG_FILE"

REQUIRED_VARS=("MAIN_DIR" "MOTIF_CLUSTER_FILE" "TIMECOURSE_EXP1" "TIMECOURSE_EXP2" \
"DHS_FILE" "BOOTSTRAP_SAMPS" "DAMAGE_REPAIR_WINDOW")
for var in "${REQUIRED_VARS[@]}"; do
  if [ -z "${!var}" ]; then
    echo "Error: Variable '$var' not set in config file."
    exit 1
  fi
done

export PYTHONPATH="${MAIN_DIR}:${PYTHONPATH}"

line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${MAIN_DIR}/${MOTIF_CLUSTER_FILE})
IFS=',' read -r TF_CLUSTER TF_LEN TF <<< "$line"

echo "Analyzing damage repair at $TF sites..."
python3 scripts/repair/analyze_TFBS_repair_sims.py \
  "$MAIN_DIR" \
  "$TIMECOURSE_EXP1" \
  "$TIMECOURSE_EXP2" \
  "${DHS_FILE%.bed}" \
  "$BOOTSTRAP_SAMPS" \
  "$DAMAGE_REPAIR_WINDOW" \
  "$TF_LEN" \
  "$TF"

echo "Finished for $TF."