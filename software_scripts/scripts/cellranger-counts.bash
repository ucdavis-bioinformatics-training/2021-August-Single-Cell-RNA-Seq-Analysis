#!/bin/bash

## Record the start time
start=`date +%s`

## Record the host being run on
echo "Hostname: $(eval hostname)"

THREADS=4
MEM=10

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

for sample in `cat samples.txt`
do
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
done

## Record the start time, and output runtime
end=`date +%s`
runtime=$((end-start))
echo $runtime
