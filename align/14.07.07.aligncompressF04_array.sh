#!/bin/bash -l
#SBATCH -J array_job
#SBATCH -o array_job_out_%A_%a.txt
#SBATCH -e array_job_err_%A_%a.txt
#SBATCH --array=1-96
#SBATCH -p med
#SBATCH --cpus-per-task=24
###### number of nodes
###### number of processors
####SBATCH -n 1

#load software for this job
module load Bowtie2/2.0.2-goolf-1.4.10;module load samtools;module load bwa/0.7.5a

#tell bowtie where the indexed genome is
export BOWTIE2_INDEXES=/home/nreid/popgen/kfish3

#commandfile contains col 1 = line number, col2 = sample, col3 = illumina barcode
cmd=/home/nreid/popgen/scripts/alignscripts3/F04cmdfile.txt

#various relevant directories and files
sourcedir=/home/nreid/popgen/140402_D00255_0089_AC3UTVACXX/Unaligned/Project_DefaultProject/
bwadir=/home/nreid/popgen/alignments/bwamem/F04/
bowtiedir=/home/nreid/popgen/alignments/bowtie/F04
bwagenind=/home/nreid/popgen/kfish3/killifish20130322asm.fa

#parts of fastq file names
root=$(sed -n "/^$SLURM_ARRAY_TASK_ID\t/p" $cmd | cut -f 2)
ind=$(sed -n "/^$SLURM_ARRAY_TASK_ID\t/p" $cmd | cut -f 3)

cd $sourcedir/Sample_$root

#decompress gzipped fastq files in current directory
gzip -d *gz

##this loop aligns each lane within a sample
for lane in {2..8}
	do
	#these lines set filenames
	fq1=$(echo $root $ind L00$lane R1 001.fastq | sed 's/ /_/g')
	fq2=$(echo $root $ind L00$lane R2 001.fastq | sed 's/ /_/g')
	rg=$(echo \@RG\\tID:$root.F04.L00$lane\\tPL:Illumina\\tPU:x\\tLB:BI_NBH_1\\tSM:$root)
	outroot=$(echo $root.F04.L00$lane.sort.rg)
		
	echo $sourcedir/Sample_$root
	echo $fq1
	echo $fq2
	echo $rg
	echo $outroot
	
	#bwa command line
	cmdline=bwa\ mem\ $bwagenind\ -t\ 24\ -R\ $rg\ $fq1\ $fq2
	echo $cmdline
		
	#execute bwa command line, pipe to samblaster to mark duplicates and create files containing discordant and split alignments, then to samtools to sort output. 
	$cmdline | /home/nreid/bin/samblaster/samblaster -e -d $bwadir/$outroot.disc.sam -s $bwadir/$outroot.split.sam | samtools view -S -h -u - | samtools sort - $bwadir/$outroot
	
	samtools view -S -h -u $bwadir/$outroot.disc.sam | samtools sort - $bwadir/$outroot.disc
	rm $bwadir/$outroot.disc.sam
	samtools view -S -h -u $bwadir/$outroot.split.sam | samtools sort - $bwadir/$outroot.split
	rm $bwadir/$outroot.split.sam
	
	#bowtie command line	
	cmdline2=bowtie2\ --very-sensitive-local\ -q\ -p\ 24\ -x\ killifish20130322asm\ --rg-id\ $root.F04.L00$lane\ --rg\ PL:Illumina\ --rg\ PU:x\ --rg\ LB:BI_NBH_1\ --rg\ SM:$root\ -1\ $fq1\ -2\ $fq2
	echo $cmdline2
	
	#execute bowtie2 command line as above, but do not output split or discordant alignments.
	$cmdline2 | /home/nreid/bin/samblaster/samblaster | samtools view -S -h -u - | samtools sort - $bowtiedir/$outroot


done

gzip *fastq

