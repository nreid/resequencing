#!/bin/bash -l
#SBATCH -J array_job
#SBATCH -o array_job_out_%A_%a.txt
#SBATCH -e array_job_err_%A_%a.txt
#SBATCH --array=1-1000
#SBATCH -p hi
#SBATCH --cpus-per-task=2
#SBATCH --mem=15000

#load software for this job
module load freebayes/0.9.14-15-gc6f49c0


#script generates jobs for scaffolds 9871-10179
	#another script generates jobs for scaffolds 0-9870

scafnum=$(expr $SLURM_ARRAY_TASK_ID + 4999)
scaf=Scaffold$scafnum
endpos=$(expr $(grep -P "$scaf\t" /home/nreid/popgen/kfish3/killifish20130322asm.fa.fai | cut -f 2) - 1)

region=$scaf:1..$endpos
echo $region

bamdir=/home/nreid/popgen/alignments/bowtie/ALL
outdir=/home/nreid/popgen/variants/bowfree/ALL1
listfile=/home/nreid/popgen/alignments/bowtie/ALL/meta/all.list
popsfile=/home/nreid/popgen/alignments/bowtie/ALL/meta/populations.txt
hicov=/home/nreid/popgen/alignments/bowtie/ALL/coverage/hicov.merge.bed

outfile=$scaf.vcf

runfb=/home/nreid/bin/freebayes/bin/freebayes 
runbamt=/home/nreid/bin/bedtools/bin//bedtools2
runbt=/home/nreid/bin/bamtools/bin/bamtools
refgen=/home/nreid/popgen/kfish3/killifish20130322asm.fa

cd $bamdir

$runbt merge -list $listfile -region $region | \
$runbt filter -in stdin -mapQuality ">30" -isProperPair true | \
$runbamt intersect -v -a stdin -b $hicov | \
$runfb -f $refgen --populations $popsfile --stdin  >$outdir/$outfile

echo $region done
