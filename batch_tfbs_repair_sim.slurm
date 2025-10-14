!/bin/bash
#SBATCH --job-name=sim_repair
#SBATCH --output=results/logs/sim_repair_%A_%a.out
#SBATCH --error=results/logs/sim_repair_%A_%a.err
#SBATCH --array=0-60
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --partition=compsci

source /home/users/hiw4/miniconda3/etc/profile.d/conda.sh
conda activate tfbs_cpd

set -e
set -o pipefail

CONFIG_FILE="$1"
TF="$2"

source "$CONFIG_FILE"

REQUIRED_VARS=("MAIN_DIR" "BOOTSTRAP_SAMPS" "DHS_FILE" "TIMECOURSE_EXP1" "TIMECOURSE_EXP2" \
"TF_WINDOW_SIZE" "DAMAGE_REPAIR_WINDOW")
for var in "${REQUIRED_VARS[@]}"; do
    if [ -z "${!var}" ]; then
        echo "Error: Variable '$var' not set in config file."
        exit 1
    fi
done

export PYTHONPATH="${MAIN_DIR}:${PYTHONPATH}"

echo "Beginning repair simulation for $TF on same strand..."
python3 scripts/repair/simulate_TFBS_repair.py "$MAIN_DIR" "$TF" "$TIMECOURSE_EXP1" "$TIMECOURSE_EXP2" \
    "same" \
    "${DHS_FILE%.bed}" \
    "$TF_WINDOW_SIZE" \
    "$BOOTSTRAP_SAMPS" \
    "$DAMAGE_REPAIR_WINDOW"

echo "Beginning repair simulation for $TF on opp strand..."
python3 scripts/repair/simulate_TFBS_repair.py "$MAIN_DIR" "$TF" "$TIMECOURSE_EXP1" "$TIMECOURSE_EXP2" \
    "opp" \
    "${DHS_FILE%.bed}" \
    "$TF_WINDOW_SIZE" \
    "$BOOTSTRAP_SAMPS" \
    "$DAMAGE_REPAIR_WINDOW"

echo "Finished repair simulations for $TF at $TIMECOURSE_EXP2..."
