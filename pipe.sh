#!/bin/sh

# Midterm assignment
# Bianca Mocanu

#===================================================================================================================================================
# Usage:
# 1. replace <username> with your BBC account username in this script
# 2. rename your input ctrl sample "Input.fastq"
# 3. Run the script by typing "bash <script path>"


#This script will analyze ChIP seq datasets and generate QC reports for the original and preprocessed data and BAM, BED and BEDGRAPH files for hg19

#===================================================================================================================================================
# Requirements:

# 1.fastq files in the /home/<username>/data prefix
# 2.indexed genome or genome file to be indexed with "bowtie-build <infile> <outfile_handle>"
# 3.fastqc module (check for latest version: fastqc/0.11.5/, retrieved on Sep. 28, 2017)
# 4.fastx_tools (check for availability on BBC in /share/apps/ - add to $PATH if needed!
# 5.bowtie2 (check for latest version: bowtie2/2.3.1/, retrieved on Sep. 28, 2017)
# 6.samtools (check for latest version: samtools/1.3.1/, retrieved on Oct. 19th, 2017)
# 7.bedtools (check for latest version: BedTools/2.26.0/, retrieved on Oct. 19th, 2017)

#===================================================================================================================================================
# Required modules load here:

module load bowtie2/2.3.1/
module load fastqc/0.11.5/
module load samtools/1.3.1/
module load BedTools/2.26.0/

#===================================================================================================================================================
# Global variables 

inPATH="/tempdata3/MCB5430/midterm/midterm/fastq/" # uncomment this for the real (very large) files
# inPATH="/home/bim16102/midterm/" #used this on 1 mil reads files to test the script
hg19index="/tempdata3/MCB5430/genomes/hg19/bowtieIndex/hg19"
hg19chromInfo="/tempdata3/MCB5430/genomes/hg19/hg19_chromInfo.txt"
gencode="/tempdata3/MCB5430/annotations/hs/bed/hg19_gencode_ENSG_geneID.bed"
chr12="/tempdata3/MCB5430/genomes/hg19/fasta/chr12.fasta"
hg19="/tempdata3/MCB5430/genomes/hg19/fasta/hg19.fasta"
TSSbackground="/tempdata3/MCB5430/midterm/hg19_unique_TSSonly_bkgrnd.txt"
outPATH="/home/bim16102/midterm/processed_data/"
jaspar_meme="/tempdata3/MCB5430/TF_db/jaspar.meme"
adapter="GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTGAAA"
fastqfiles=$(find ${inPATH} -maxdepth 1 -type f)
#===================================================================================================================================================

#if [ -s $outPATH ]
#	then
#		cd ${outPATH}
#		mkdir ${outPATH}logfiles
#		touch ${outPATH}logfiles/log.txt
#		echo "$outPATH directory already exists" | tee -a ${outPATH}logfiles/log.txt
#
#	else
#		mkdir ${outPATH}
#		cd ${outPATH}
#		mkdir ${outPATH}logfiles
#		touch ${outPATH}logfiles/log.txt
#		echo "New directory created: ${outPATH}" | tee -a ${outPATH}logfiles/log.txt
#fi |


#echo -e "Files to be proceesed: $fastqfiles" | tee -a ${outPATH}logfiles/log.txt

for file in $fastqfiles
	do
		ext=`echo $(basename $file) | cut -d "." -f 2` # generated to see file type
		prefix=`echo $(basename $file) | cut -d "." -f 1`  #creates a prefix for each fastq file that is analyzed

#		if [ $ext=="fastq" ]
#		then
#
#		mkdir ${outPATH}${prefix}  #folder for the fastq file / sample
		cd ${outPATH}$prefix
#
#		echo -e "Starting analysis on $(basename $file) ..."
#
#
#		echo "Generating QC reports of unprocessed $(basename $file)"
#		mkdir ./unprocessed_data_qc/
#		fastqc $file -o ./unprocessed_data_qc  2>&1
#
#		echo "Clipping adapter sequences..."
#		fastx_clipper -Q 33 -a $adapter -i $file -o ./${prefix}_clipped.fastq 2>&1
#
#
#		echo "Trimming low quality bases..."
#		fastq_quality_trimmer -Q33 -t 32 -l 30 -i ./${prefix}_clipped.fastq -o ./${prefix}_preprocessed.fastq 2>&1
#
#		echo "Generating QC reports of the preprocessed $(basename $file)..."
#		mkdir ./preprocessed_data_qc
#		fastqc ${prefix}_preprocessed.fastq -o ./preprocessed_data_qc/ 2>&1
#
#		echo "Aligning to Human genome (hg19) ..."
#		bowtie -S -v0 -m1 -t -q $hg19index ${prefix}_preprocessed.fastq ./${prefix}.sam 2>&1
#		cat ./${prefix}.sam | head -n 27  > ./${prefix}_chr12.sam # this appends the sam header to the chromosome 12 only sam file
#
#		echo "Subsetting the data to chromosome 12 aligned reads"
#		grep chr12 ./${prefix}.sam >> ./${prefix}_chr12.sam
#		echo "Alignment to human genome (hg19) complete for $(basename $file)!"
#
#		echo "Generating BAM file..."
#		samtools view -S -b ${prefix}_chr12.sam > ${prefix}_chr12.bam  # sam to bam conversion
#
#		samtools sort -l 9 -n  ${prefix}_chr12.bam -T ${prefix} -o ${prefix}_chr12.sorted.bam  #this sorts the bam file so that it occupies less space
#		echo "BAM file generated!"
#
#		echo "Generating BED file..."
#		bedtools bamtobed -i ${prefix}_chr12.sorted.bam > ${prefix}_chr12.bed
#		echo "BED file generated!"
#		echo "Sorting BED file"
#		sortBed -i ${prefix}_chr12.bed > ${prefix}_chr12_sorted.bed
#		echo "Generating bedgraph file..."
#		bedtools genomecov -ibam ${prefix}_chr12.sorted.bam -bg > ${prefix}.bedgraph #generates the bedgraph from bam directly
#		echo "Bedgraph file generated!"

#==============================================================================================================================================================
# Comment this section if you want to keep all these files
#==============================================================================================================================================================
#		echo "Cleaning up temporary files"
#		rm ${prefix}_clipped.fastq # partially processed file
#		rm ${prefix}_preprocessed.fastq # fully preprocessed fastq file - it occupies a lot of space
#		rm ${prefix}.sam # large file as well, virtually useless once converted to bam
#		rm ${prefix}_chr12.sam # subset of the file above, only for the assigned chromosome
#		rm ${prefix}_chr12.bam # unsorted bam file
#		rm ${prefix}_chr12.bed # unsorted bed file

#==============================================================================================================================================================

#		echo "Preparing bedgraphs for Genome Browser..."
#		if [ $prefix=="treatA_chip_rep1" ] || [ $prefix=="treatA_chip_rep2" ]
#		then
#			awk -v NAME="$prefix" 'BEGIN { print "browser position chr12:5,289,521-5,291,937"
#			print "track type=bedGraph name=\""NAME"\" description=\""NAME"\" visibility=full windowingFunction=maximum color=0,0,125"}
#			{print $0}' ${prefix}.bedgraph > ${prefix}_header.bedgraph
#		elif [ $prefix=="treatAB_chip_rep1" ] || [ $prefix=="treatAB_chip_rep2" ]
#		then
#			awk -v NAME=$prefix 'BEGIN { print "browser position chr12:5,289,521-5,291,937"
#			print "track type=bedGraph name=\""NAME"\" description=\""NAME"\" visibility=full windowingFunction=maximum color=125,0,125"}
#			{print $0}' ${prefix}.bedgraph > ${prefix}_header.bedgraph
#		elif [ $prefix=="Input" ]
#		then
#			awk -v NAME=$prefix 'BEGIN { print "browser position chr12:5,289,521-5,291,937"
#			print "track type=bedGraph name=\""NAME"\" description=\""NAME"\" visibility=full windowingFunction=maximum color=125,0,0"}
#			{print $0}' ${prefix}.bedgraph > ${prefix}_header.bedgraph
#		fi
#		echo "Genome Browser bedgraphs generated!"

#=============================================================================================================================================================
# This part is for calling peaks using MACS. After peak calling, it shifts the peaks by half the "d" value that the pdf reports and creates files for genome
# browser use by adding the BED headers. The genome browser will display a region 300 nts upstream and downstream of the top peak found in chromosome 12.
#For MEME and FIMO usage, one should use the top summits from the entire set of chromosomes (w/ file provided on the server)
#=============================================================================================================================================================
		cd ..

		echo "Calling peaks for Chromosome 12 using MACS"
		if [ -s ${outPATH}peaks ]
		then
		cd peaks
		else
		mkdir peaks
		cd peaks
		fi

		if [ $prefix != "Input" ]
		then
			macs14 -t ${outPATH}${prefix}/${prefix}_chr12.sorted.bam -c ${outPATH}Input/Input_chr12.sorted.bam -f BAM -n ${prefix} -g 133851895
			Rscript ${prefix}_model.r
			peakshift=`grep "legend" ${prefix}_model.r | tail -n 1 | cut -d "=" -f2 | cut -d "\"" -f1`

			top_peak=`sort -k5nr ${prefix}_summits.bed | head -1 | cut -f2`
			browser_start=$(($top_peak - 300))
			browser_end=$(($top_peak + 300))

			echo "Shifting peaks by $peakshift"
			awk -v d=$peakshift '{printf ("%s\t%s\t%s\t%s\t%s\n", $1, $2 + (d/2), $3 - (d/2), $4, $5)}' ${prefix}_peaks.bed > ${prefix}_peaks_shifted.bed

			echo "Generating UCSC BED files with headers for peaks and summits"

			awk -v NAME=${prefix}_peaks -v browser_start=$browser_start -v browser_end=$browser_end 'BEGIN { print "browser position chr12:("browser_start")-("browser_end")"
			print "track type=bed name=\""NAME"\" description=\""NAME"\" visibility=squish autoScale=on colorByStrand=\"255,0,0 0,0,255\""}
			{ print $0}' ${prefix}_peaks_shifted.bed > ${prefix}_peaks_shifted_header.bed

			awk -v NAME=${prefix}_summits -v browser_start=$browser_start -v browser_end=$browser_end 'BEGIN { print "browser position chr12:("browser_start")-("browser_end")"
                        print "track type=bed name=\""NAME"\" description=\""NAME"\" visibility=squish autoScale=on colorByStrand=\"255,0,0 0,0,255\""}
                        { print $0}' ${prefix}_summits.bed > ${prefix}_summits_header.bed

#==============================================================================================================================================================
# This part intersects the datasets to report:
# 1. High confidence peaks between replicates
# 2. Peaks specific only to treatment A or only to treatment A+B
#==============================================================================================================================================================

			sample=`echo $prefix | cut -d "_" -f1,2`
			if [ -s ${sample}_rep1_peaks_shifted.bed ] && [ -s ${sample}_rep2_peaks_shifted.bed ]
			then
			echo "Finding high confidence peaks between replicates"
			bedtools intersect -a ${sample}_rep1_peaks_shifted.bed -b ${sample}_rep2_peaks_shifted.bed > ${sample}_peaks_highconf.bed

# the next statement is a bit iffy because it needs the other sample high confidence peaks bed file, and it depends on the order the files are processed; should work though

				if [ $sample=="treatA_chip" ] && [ -s treatAB_chip_peaks_highconf.bed ]
				then
				bedtools intersect -v -a ${sample}_peaks_highconf.bed -b treatAB_chip_peaks_highconf.bed > ${sample}_only_peaks.bed
				elif [ $sample=="treatAB_chip"] && [ -s treatA_chip_peaks_highconf.bed ]
				then
				bedtools intersect -v -a ${sample}_peaks_highconf.bed -b treatA_chip_peaks_highconf.bed > ${sample}_only_peaks.bed
				fi
			fi
		fi
	done | tee -a ${outPATH}logfiles/log.txt


#==============================================================================================================================================================
# Unloading of required modules:

module unload bowtie2/2.3.1/
module unload fastqc/0.11.5/
module unload samtools/1.3.1/
module unload BedTools/2.26.0/
