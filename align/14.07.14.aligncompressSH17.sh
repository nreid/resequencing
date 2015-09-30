#!/bin/bash

#load software for this job
module load Bowtie2/2.0.2-goolf-1.4.10;module load samtools;module load bwa/0.7.5a

#tell bowtie where the indexed genome is
export BOWTIE2_INDEXES=/home/nreid/popgen/kfish3

bwadir=/home/nreid/popgen/alignments/bwamem/P020304/
bowtiedir=/home/nreid/popgen/alignments/bowtie/P020304/
bwagenind=/home/nreid/popgen/kfish3/killifish20130322asm.fa

cd /home/nreid/popgen/UNH/140620/GSF642_2/Sample_SH-17
	
gzip -d *gz

	outroot1=$(echo SH-17.F010.L001.sort.rg)	
	outroot2=$(echo SH-17.F010.L002.sort.rg)	
	rg1=$(echo \@RG\\tID:SH-17.F010.L001\\tPL:Illumina\\tPU:x\\tLB:BP_F_NYC_SH\\tSM:SH-17)
	rg2=$(echo \@RG\\tID:SH-17.F010.L002\\tPL:Illumina\\tPU:x\\tLB:BP_F_NYC_SH\\tSM:SH-17)
	
	#execute bwa command line, pipe to samblaster to mark duplicates and create files containing discordant and split alignments, then to samtools to sort output. 
	bwa mem $bwagenind -t 24 -R $rg1 '<cat SH-17_AGGCTAAC_L001_R1_001.fastq SH-17_AGGCTAAC_L001_R1_002.fastq' '<cat SH-17_AGGCTAAC_L001_R2_001.fastq SH-17_AGGCTAAC_L001_R2_002.fastq'| /home/nreid/bin/samblaster/samblaster -e -d $bwadir/$outroot1.disc.sam -s $bwadir/$outroot1.split.sam | samtools view -S -h -u - | samtools sort - $bwadir/$outroot1
	bwa mem $bwagenind -t 24 -R $rg2 '<cat SH-17_AGGCTAAC_L002_R1_001.fastq SH-17_AGGCTAAC_L002_R1_002.fastq' '<cat SH-17_AGGCTAAC_L002_R2_001.fastq SH-17_AGGCTAAC_L002_R2_002.fastq'| /home/nreid/bin/samblaster/samblaster -e -d $bwadir/$outroot2.disc.sam -s $bwadir/$outroot2.split.sam | samtools view -S -h -u - | samtools sort - $bwadir/$outroot2

	samtools view -S -h -u $bwadir/$outroot1.disc.sam | samtools sort - $bwadir/$outroot1.disc
	rm $bwadir/$outroot1.disc.sam
	samtools view -S -h -u $bwadir/$outroot1.split.sam | samtools sort - $bwadir/$outroot1.split
	rm $bwadir/$outroot1.split.sam

	samtools view -S -h -u $bwadir/$outroot2.disc.sam | samtools sort - $bwadir/$outroot2.disc
	rm $bwadir/$outroot2.disc.sam
	samtools view -S -h -u $bwadir/$outroot2.split.sam | samtools sort - $bwadir/$outroot2.split
	rm $bwadir/$outroot2.split.sam

	#execute bowtie2 command line as above, but do not output split or discordant alignments.
	bowtie2 --very-sensitive-local -q -p 24 -x killifish20130322asm --rg-id SH-17.F010.L001 --rg PL:Illumina --rg PU:x --rg LB:BP_F_NYC_SH --rg SM:SH-17 -1 SH-17_AGGCTAAC_L001_R1_001.fastq,SH-17_AGGCTAAC_L001_R1_002.fastq -2 SH-17_AGGCTAAC_L001_R2_001.fastq,SH-17_AGGCTAAC_L001_R2_002.fastq | /home/nreid/bin/samblaster/samblaster | samtools view -S -h -u - | samtools sort - $bowtiedir/$outroot1
	bowtie2 --very-sensitive-local -q -p 24 -x killifish20130322asm --rg-id SH-17.F010.L002 --rg PL:Illumina --rg PU:x --rg LB:BP_F_NYC_SH --rg SM:SH-17 -1 SH-17_AGGCTAAC_L002_R1_001.fastq,SH-17_AGGCTAAC_L002_R1_002.fastq -2 SH-17_AGGCTAAC_L002_R2_001.fastq,SH-17_AGGCTAAC_L002_R2_002.fastq | /home/nreid/bin/samblaster/samblaster | samtools view -S -h -u - | samtools sort - $bowtiedir/$outroot2
	
	done
	
gzip *fastq
