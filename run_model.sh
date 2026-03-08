#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=64GB
#SBATCH --job-name=cabac_sig_gt1
#SBATCH --output=/scratch/zr2245/VVCSoftware_VTM-master/exp/bbc1/result/cabac_sig_gt1_%j.out
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --account=torch_pr_38_tandon_advanced
#SBATCH --constraint="l40s"

module purge

cd /scratch/zr2245/VVCSoftware_VTM-master

mkdir -p /scratch/zr2245/VVCSoftware_VTM-master/exp/bbc1/result

PYTHON_BIN=/scratch/zr2245/anaconda3/envs/gs/bin/python
$PYTHON_BIN -c "import sys, torch; print(sys.executable); print(torch.__version__); print(torch.cuda.is_available())"

echo "Running on host: $(hostname)"
echo "CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES"
echo "Using python: $PYTHON_BIN"
echo "Start time: $(date)"

$PYTHON_BIN /scratch/zr2245/VVCSoftware_VTM-master/train_cabac_sig_gt1.py

echo "End time: $(date)"