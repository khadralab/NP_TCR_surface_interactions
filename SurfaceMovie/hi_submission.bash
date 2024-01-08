#!/bin/bash
#SBATCH --account=def-akhadra

#SBATCH --array=13-15
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G

#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=louis.richez@mail.mcgill.ca

module load nixpkgs/16.09
module load matlab/2018b
# Remove -singleCompThread below if you are using parallel commands:
srun -c $SLURM_CPUS_PER_TASK matlab -nosplash -nodesktop -nodisplay -r  "DoseResponseFunction($SLURM_ARRAY_TASK_ID)" > stdout 2>&1
