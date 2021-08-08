#!/bin/bash
#SBATCH --time=0-1  # days-hours
#SBATCH --job-name=cellrngr  # Job name
#SBATCH --nodes=1
#SBATCH --ntasks=4 # Number of cores
#SBATCH --mem=10000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production # Partition to submit to
#SBATCH --reservation=workshop
#SBATCH --account=workshop
#SBATCH --output=counts-cellrngr.out # File to which STDOUT will be written
#SBATCH --error=counts-cellrngr.err # File to which STDERR will be written
#SBATCH --mail-type=END # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=youremail@whatever.edu # Email to which notifications will be sent

## Record the start time
start=`date +%s`

## Record the host being run on
echo "Hostname: $(eval hostname)"


THREADS=${SLURM_NTASKS}
MEM=$(expr ${SLURM_MEM_PER_NODE} / 1024)

echo "Allocated threads: " $THREADS
echo "Allocated memory: " $MEM

## Where cellranger executable is located
## a) by loading a module
#module load cellranger/6.0.1

## b) or, by placing the location of the executables on the path (edit to your location)
export PATH=/share/workshop/intro_scrnaseq/software/cellranger-6.0.2/bin:$PATH

## c) or if they are already on the path, do nothing

## Set the parameters for the run
basedir="/share/workshop/intro_scrnaseq"
transcriptome=${basedir}/software/refdata-gex-GRCh38-2020-A
fastqs="/share/workshop/intro_scrnaseq/${USER}/scrnaseq_example/00-RawData"


## provide the script the row # of the sample to be run
sample=`sed "$1q;d" samples.txt`

## https://support.10xgenomics.com/single-cell-gene-expression/software/overview/welcome
## Create the call
call="cellranger count \
  --id=${sample} \
  --sample=${sample} \
  --transcriptome=${transcriptome} \
  --fastqs=${fastqs} \
  --localcores=${THREADS} \
  --localmem=${MEM}"

## Some other parameters that may be usefull/needed
## --expect-cells=NUM, number of cells expected
## --include-introns         Include intronic reads in count
## --nosecondary, skip the unnecessary secondary analysis
## --r2-length=NUM, if your R2 qualities are really poor
## --chemistry=CHEM, should it fail chemistry detection

## Echo the call
echo $call
## Evaluate the call
#eval $call

## Record the start time, and output runtime
end=`date +%s`
runtime=$((end-start))
echo $runtime
