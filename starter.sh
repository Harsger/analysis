#!/usr/bin/env bash
#
## run root benchmark
#
## Slurm parameters
#SBATCH --job-name=SM2inCRF
#SBATCH --partition=lsschaile
#SBATCH --export=NONE  # do not export current environemnt to job(recommended)
#SBATCH --get-user-env
#SBATCH --mem=3000mb
#SBATCH --time=01:00:00
#SBATCH -o /home/m/Maximilian.Herrmann/logslurm/job%j.txt   # out file name (in submission dir)

# module load marabou/6.10.02

# source /software/opt/xenial/x86_64/root/6.10.02/bin/thisroot.sh

source /etc/profile.d/modules.sh

module load root

command=""

for input; do
    command="$command ${input}"
done

echo " @ $(hostname) : $command"

date +%Y-%m-%dT%H:%M:%S
$command
date +%Y-%m-%dT%H:%M:%S

