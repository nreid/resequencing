#!/bin/bash -l
#SBATCH -J array_job
#SBATCH -o array_job_out_%A_%a.txt
#SBATCH -e array_job_err_%A_%a.txt
#SBATCH --array=3-96
#SBATCH -p med
#SBATCH --cpus-per-task=24
###### number of nodes
###### number of processors
####SBATCH -n 1

###align, mdup, compress UNH data. 

#load software for this job
module load Bowtie2/2.0.2-goolf-1.4.10;module load samtools;module load bwa/0.7.5a

#tell bowtie where the indexed genome is
export BOWTIE2_INDEXES=/home/nreid/popgen/kfish3

#various relevant directories and files
flowcell5=/home/nreid/popgen/UNH/140612/GSF642_1
flowcell6=/home/nreid/popgen/UNH/140612/GSF642_3
flowcell7=/home/nreid/popgen/UNH/140614/GSF642_1
flowcell8=/home/nreid/popgen/UNH/140614/GSF642_3
flowcell9=/home/nreid/popgen/UNH/140617/GSF642_2
flowcell10=/home/nreid/popgen/UNH/140620/GSF642_2

flowcelllist=$flowcell5\ $flowcell6\ $flowcell7\ $flowcell8\ $flowcell9\ $flowcell10

bwadir=/home/nreid/popgen/alignments/bwamem/P020304/
bowtiedir=/home/nreid/popgen/alignments/bowtie/P020304/
bwagenind=/home/nreid/popgen/kfish3/killifish20130322asm.fa

#slurm loops 1:96

#in each array job:

	for flowcell in {1..6}
		do		
		echo $flowcell
		flowdir=$(echo $flowcelllist | cut -f $flowcell -d ' ')
		cd $flowdir

		lib=$(echo $(ls | sed 's/Sample_//' | sed -r 's/-[0-9]+//' | sort | uniq ) | sed 's/[ \t]/_/g')
		sample=$(ls | sed -n "$SLURM_ARRAY_TASK_ID p" | sed 's/Sample_//')
		sampledir=$(ls | sed -n "$SLURM_ARRAY_TASK_ID p")
		cd $sampledir
		
		echo $sample
		
		gzip -d *gz
		
		for lane in $(ls | grep -o 'L00.' | uniq) 
			do fq1=$(ls | grep $lane| sed -n '1p')
			fq2=$(ls | grep $lane| sed -n '2p')			
			fcell=$(expr 4 + $flowcell)
			outroot=$(echo $sample.F0$fcell.$lane.sort.rg)	
			rg=$(echo \@RG\\tID:$sample.F0$fcell.$lane\\tPL:Illumina\\tPU:x\\tLB:$lib\\tSM:$sample)
			echo $fq1 $fq2 $outroot 
			echo $rg
			
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
			cmdline2=bowtie2\ --very-sensitive-local\ -q\ -p\ 24\ -x\ killifish20130322asm\ --rg-id\ $sample.F0$fcell.$lane\ --rg\ PL:Illumina\ --rg\ PU:x\ --rg\ LB:$lib\ --rg\ SM:$sample\ -1\ $fq1\ -2\ $fq2
			echo $cmdline2

			#execute bowtie2 command line as above, but do not output split or discordant alignments.
			$cmdline2 | /home/nreid/bin/samblaster/samblaster | samtools view -S -h -u - | samtools sort - $bowtiedir/$outroot
			
			done
			
		gzip *fastq
		
		done