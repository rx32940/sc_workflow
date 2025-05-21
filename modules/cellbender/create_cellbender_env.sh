#BSUB -P cellbender_gpu
#BSUB -q gpu
#BSUB -gpu "num=1/host"
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=30GB]"
#BSUB -J cellbender_gpu
#BSUB -o cellbender_gpu.%J.out
#BSUB -e cellbender_gpu.%J.out

# cellbender need to be installed in a machine with a GPU to install the appropriate version of PyTorch with CUDA support
conda create -p /research/groups/northcgrp/home/common/Rachel/envs/cellbender python=3.7
source $HOME/miniconda3/bin/activate /research/groups/northcgrp/home/common/Rachel/envs/cellbender

pip install cellbender

conda $HOME/miniconda3/bin/deactivate
