#!/bin/bash -l
#SBATCH -J array_job
#SBATCH -o array_job_out_%A_%a.txt
#SBATCH -e array_job_err_%A_%a.txt
#SBATCH --array=1-96
#SBATCH -p med
#SBATCH --cpus-per-task=2
###### number of nodes
###### number of processors
####SBATCH -n 1

#merges BWA bam files

#load software for this job
module load Bowtie2/2.0.2-goolf-1.4.10;module load samtools;module load bwa/0.7.5a

#various relevant directories containing alignments
newps=/home/nreid/popgen/alignments/bwamem/F11_12_13_14/

#output directories
discdir=/home/nreid/popgen/alignments/bwamem/discordant
splitdir=/home/nreid/popgen/alignments/bwamem/split
outdir=/home/nreid/popgen/alignments/bwamem/ALL

#root name for alignments
root=$(ls $newps | sed 's/\..*//' | sort | uniq | sed -n "$SLURM_ARRAY_TASK_ID p" )
echo $root

FULL=I=$newps/$root.all.bam
DISC=I=$newps/$root.disc.bam
SPLIT=I=$newps/$root.split.bam

cd $newps
for file in $(ls | grep $root.F | grep .rg.bam); do FULL=$FULL\ I=$newps/$file; done
for file in $(ls | grep $root.F | grep .disc.bam); do DISC=$DISC\ I=$newps/$file; done
for file in $(ls | grep $root.F | grep .split.bam); do SPLIT=$SPLIT\ I=$newps/$file; done


java -jar /home/nreid/bin/picard-tools-1.96/MergeSamFiles.jar $FULL O=$outdir/$root.all.bam

java -jar /home/nreid/bin/picard-tools-1.96/MergeSamFiles.jar $DISC O=$discdir/$root.disc.bam

java -jar /home/nreid/bin/picard-tools-1.96/MergeSamFiles.jar $SPLIT O=$splitdir/$root.split.bam

