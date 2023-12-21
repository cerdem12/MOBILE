#!/bin/bash
#SBATCH -A hpc2n2023-130
#SBATCH -J 1203PRADERG
#SBATCH -n 16
#SBATCH -t 12:00:00
#SBATCH --error=matlab_%J.err
#SBATCH --output=matlab_%J.out

module add MATLAB/2023a.Update4
module add GCCcore/11.3.0
module add X11/20220504
module add binutils/2.38

matlab -nojvm -nodisplay -r "PRADERG"
