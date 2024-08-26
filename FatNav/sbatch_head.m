function sbatch_head(fid,name,mem,time)

fprintf(fid,'#!/bin/bash\n');

fprintf(fid,'
#SBATCH --job-name=recon
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=400000
#SBATCH --partition=bigmem
#SBATCH --qos bigmem_access
