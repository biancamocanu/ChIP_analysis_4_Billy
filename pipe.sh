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
chr12="/tempdata3/MCB5430/genomes/hg19/fasta/chr12.fa"
hg19="/tempdata3/MCB5430/genomes/hg19/fasta/hg19.fa"
TSSbackground="/tempdata3/MCB5430/midterm/midterm/hg19_unique_TSSonly_bkgrnd.txt"
outPATH="/home/bim16102/midterm/processed_data/"
jaspar_meme="/tempdata3/MCB5430/TF_db/jaspar.meme"
adapter="GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTGAAA"
summits_highconf=$(find /tempdata3/MCB5430/midterm/midterm/peaks/ -maxdepth 1 -type f)
fastqfiles=$(find ${inPATH} -maxdepth 1 -type f)
#===================================================================================================================================================

if [ -s $outPATH ]
	then
		cd ${outPATH}
		mkdir ${outPATH}logfiles
		touch ${outPATH}logfiles/log.txt
		echo "$outPATH directory already exists" | tee -a ${outPATH}logfiles/log.txt

	else
		mkdir ${outPATH}
		cd ${outPATH}
		mkdir ${outPATH}logfiles
		touch ${outPATH}logfiles/log.txt
		echo "New directory created: ${outPATH}" | tee -a ${outPATH}logfiles/log.txt
fi


echo -e "Files to be proceesed: $fastqfiles" | tee -a ${outPATH}logfiles/log.txt

for file in $fastqfiles
	do
		ext=`echo $(basename $file) | cut -d "." -f 2` # generated to see file type
		prefix=`echo $(basename $file) | cut -d "." -f 1`  #creates a prefix for each fastq file that is analyzed

		if [ $ext=="fastq" ]
		then

		mkdir ${outPATH}${prefix}  #folder for the fastq file / sample
		cd ${outPATH}$prefix

		echo -e "Starting analysis on $(basename $file) ..."


		echo "Generating QC reports of unprocessed $(basename $file)"
		mkdir ./unprocessed_data_qc/
		fastqc $file -o ./unprocessed_data_qc  2>&1

		echo "Clipping adapter sequences..."
		fastx_clipper -Q 33 -a $adapter -i $file -o ./${prefix}_clipped.fastq 2>&1


		echo "Trimming low quality bases..."
		fastq_quality_trimmer -Q33 -t 32 -l 30 -i ./${prefix}_clipped.fastq -o ./${prefix}_preprocessed.fastq 2>&1

		echo "Generating QC reports of the preprocessed $(basename $file)..."
		mkdir ./preprocessed_data_qc
		fastqc ${prefix}_preprocessed.fastq -o ./preprocessed_data_qc/ 2>&1

		echo "Aligning to Human genome (hg19) ..."
		bowtie -S -v0 -m1 -t -q $hg19index ${prefix}_preprocessed.fastq ./${prefix}.sam 2>&1
		cat ./${prefix}.sam | head -n 27  > ./${prefix}_chr12.sam # this appends the sam header to the chromosome 12 only sam file

		echo "Subsetting the data to chromosome 12 aligned reads"
		grep chr12 ./${prefix}.sam >> ./${prefix}_chr12.sam
		echo "Alignment to human genome (hg19) complete for $(basename $file)!"

		echo "Generating BAM file..."
		samtools view -S -b ${prefix}_chr12.sam > ${prefix}_chr12.bam  # sam to bam conversion

		samtools sort -l 9 -n  ${prefix}_chr12.bam -T ${prefix} -o ${prefix}_chr12.sorted.bam  #this sorts the bam file so that it occupies less space
		echo "BAM file generated!"

		echo "Generating BED file..."
		bedtools bamtobed -i ${prefix}_chr12.sorted.bam > ${prefix}_chr12.bed
		echo "BED file generated!"
		echo "Sorting BED file"
		sortBed -i ${prefix}_chr12.bed > ${prefix}_chr12_sorted.bed
		echo "Generating bedgraph file..."
		bedtools genomecov -ibam ${prefix}_chr12.sorted.bam -bg > ${prefix}.bedgraph #generates the bedgraph from bam directly
		echo "Bedgraph file generated!"

#==============================================================================================================================================================
# Comment this section if you want to keep all these files
#==============================================================================================================================================================
		echo "Cleaning up temporary files"
		rm ${prefix}_clipped.fastq # partially processed file
		rm ${prefix}_preprocessed.fastq # fully preprocessed fastq file - it occupies a lot of space
		rm ${prefix}.sam # large file as well, virtually useless once converted to bam
		rm ${prefix}_chr12.sam # subset of the file above, only for the assigned chromosome
		rm ${prefix}_chr12.bam # unsorted bam file
		rm ${prefix}_chr12.bed # unsorted bed file

#==============================================================================================================================================================

		echo "Preparing bedgraphs for Genome Browser..."
		if [ $prefix=="treatA_chip_rep1" ] || [ $prefix=="treatA_chip_rep2" ]
		then
			awk -v NAME="$prefix" 'BEGIN { print "browser position chr12:5,289,521-5,291,937"
			print "track type=bedGraph name=\""NAME"\" description=\""NAME"\" visibility=full windowingFunction=maximum color=0,0,125"}
			{print $0}' ${prefix}.bedgraph > ${prefix}_header.bedgraph
		elif [ $prefix=="treatAB_chip_rep1" ] || [ $prefix=="treatAB_chip_rep2" ]
		then
			awk -v NAME=$prefix 'BEGIN { print "browser position chr12:5,289,521-5,291,937"
			print "track type=bedGraph name=\""NAME"\" description=\""NAME"\" visibility=full windowingFunction=maximum color=125,0,125"}
			{print $0}' ${prefix}.bedgraph > ${prefix}_header.bedgraph
		elif [ $prefix=="Input" ]
		then
			awk -v NAME=$prefix 'BEGIN { print "browser position chr12:5,289,521-5,291,937"
			print "track type=bedGraph name=\""NAME"\" description=\""NAME"\" visibility=full windowingFunction=maximum color=125,0,0"}
			{print $0}' ${prefix}.bedgraph > ${prefix}_header.bedgraph
		fi
		echo "Genome Browser bedgraphs generated!"

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

# the next statement is a bit iffy because it needs the other sample high confidence peaks bed file, and it depends on the order the files are processed, but works

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
cd $outPATH
echo "Generating gene lists" | tee -a ${outPATH}logfiles/log_geneLists.txt
mkdir geneLists
cd geneLists


#==============================================================================================================================================================
# This part processes the hg19_gencode_ENSG_geneID.bed file to retrieve the TSS, promoters, genes and intergenic regions for chromosome 12
#==============================================================================================================================================================


# Retrieving only chromosome 12 entries:

echo "Retrieving chr12 entries"  | tee -a ${outPATH}logfiles/log_geneLists.txt
cat $gencode | grep "chr12" > ./gencode_ENSG_geneID_chr12.txt

# Retrieving TSS - if the gene is on the + strand, start of gene is $2 => subtract 500 from $2 (start of TSS) and add 500 to $2 (end of TSS)
# If the gene is on the - strand, the start of the gene is $3 => add 500 to $3 (start of TSS) and subtract 500 to $3 (end of TSS)
# In these BED files, the smaller coordinate is always in col. 2, regardless of the strand

echo "Creating TSS only file" | tee -a ${outPATH}logfiles/log_geneLists.txt
awk '{if($6=="+" && $2<$3)
printf ("%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2 - 500, $2 + 500, $4, $5, $6);
else if($6=="-" && $2<$3)
printf ("%s\t%s\t%s\t%s\t%s\t%s\n", $1, $3 - 500, $3 + 500, $4, $5, $6)
}' gencode_ENSG_geneID_chr12.txt > gencode_ENSG_geneID_chr12_TSS.bed
bedtools getfasta -name -fi $hg19 -bed gencode_ENSG_geneID_chr12_TSS.bed -fo gencode_ENSG_geneID_chr12_TSS.fasta

# Retrieving genes - if the gene is on the + strand, gene region starts at $2 + 501 and ends at the $3 + 1000
# if the gene is on the - strand, the gene ends at $2 - 1000, starts at $3 - 501
# In these BED files, the smaller coordinate is always in col. 2 regardless of the strand

echo "Creating genes only file" | tee -a ${outPATH}logfiles/log_geneLists.txt
awk '{if($6=="+" && $2<$3)
printf ("%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2 + 501, $3 + 1000, $4, $5, $6);
else if($6=="-" && $2<$3)
printf ("%s\t%s\t%s\t%s\t%s\t%s\n", $1, $2 - 1000, $3 - 501, $4, $5, $6)
}' gencode_ENSG_geneID_chr12.txt > gencode_ENSG_geneID_chr12_genes.bed
bedtools getfasta -name -fi $hg19 -bed gencode_ENSG_geneID_chr12_genes.bed -fo gencode_ENSG_geneID_chr12_genes.fasta


# Intergenic regions - basically if it's not part of the first two - intersect with chromosome 12 file and take the complement

echo "Creating IGS only file" | tee -a ${outPATH}logfiles/log_geneLists.txt
#creating a temporary file with both TSS and genes and creating the chromosome 12 only file size

grep chr12 $hg19chromInfo > ./chr12Info.txt

cat ./gencode_ENSG_geneID_chr12_TSS.bed >> gencode_ENSG_geneID_chr12_genesandTSS.bed
cat ./gencode_ENSG_geneID_chr12_genes.bed >> gencode_ENSG_geneID_chr12_genesandTSS.bed

# intersecting the file with chromosome 12
bedtools sort -i gencode_ENSG_geneID_chr12_genesandTSS.bed > gencode_ENSG_geneID_chr12_genesandTSS.sorted.bed
bedtools complement -i gencode_ENSG_geneID_chr12_genesandTSS.sorted.bed -g chr12Info.txt > gencode_ENSG_geneID_chr12_IGS.bed
rm gencode_ENSG_geneID_chr12_genesandTSS.*
bedtools getfasta -name -fi $hg19 -bed gencode_ENSG_geneID_chr12_IGS.bed -fo gencode_ENSG_geneID_chr12_IGS.fasta


#==============================================================================================================================================================
# This part analyzes the chromosome 12 summits from treatments A and A+B in order to see the distribution of the peaks that fall in the TSS, genes and IGS regions
# The summits files are the ones provided on /tempdata3/MCB5430/midterm/midterm/peaks folder which are the from the entire genome
# Also generates fasta files for MEME motif analysis and finds the MEME motifs
#==============================================================================================================================================================

for file in $summits_highconf
	do
		ext=`echo $(basename $file) | cut -d "." -f 2` # generated to see file type
		prefix=`echo $(basename $file) | cut -d "." -f 1`  #creates a prefix for each bed file that is analyzed
		echo "Starting analysis for high confidence summits for $prefix"
		grep chr12 $file > ${prefix}_chr12.bed

		echo "Examining distribution in TSS, IGS, genes..."
		if [ $ext=="bed" ]
		then
		bedtools coverage -a ${prefix}_chr12.bed -b gencode_ENSG_geneID_chr12_TSS.bed  > ${prefix}_chr12_inTSS.bed
		bedtools coverage -a ${prefix}_chr12.bed -b gencode_ENSG_geneID_chr12_IGS.bed > ${prefix}_chr12_inIGS.bed
		bedtools coverage -a ${prefix}_chr12.bed -b gencode_ENSG_geneID_chr12_genes.bed > ${prefix}_chr12_ingenes.bed

		awk '{if($9!="0.0000000")
		print $0}' ${prefix}_chr12_inTSS.bed > ${prefix}_chr12_inTSS_nozero.bed

		awk '{if($9!="0.0000000")
		print $0}' ${prefix}_chr12_ingenes.bed > ${prefix}_chr12_ingenes_nozero.bed

		awk '{if($9!="0.0000000")
		print $0}' ${prefix}_chr12_inIGS.bed > ${prefix}_chr12_inIGS_nozero.bed

		fi

		echo "Retrieving top 200 peaks from the entire genome"
		sort -k5nr $file | head -n 200 > ${prefix}_top200.bed
		bedtools slop -i ${prefix}_top200.bed -g $hg19chromInfo -b 50 > ${prefix}_top200_100bp.bed
		bedtools getfasta -name -fi $hg19 -bed ${prefix}_top200_100bp.bed -fo ${prefix}_top200_100bp.fasta

		echo "Finding motifs with MEME for $prefix"
		if [ ${prefix}=="treatA_summits" ]
		then
			meme ${prefix}_top200_100bp.fasta -oc ${prefix}_meme_OUT_FOLDER -bfile $TSSbackground -dna -nmotifs 2 -minw 10 -maxw 18 -revcomp -mod anr
		else
			meme ${prefix}_top200_100bp.fasta -oc ${prefix}_meme_OUT_FOLDER -bfile $TSSbackground -dna -nmotifs 2 -minw 12 -maxw 14 -revcomp -mod anr
		fi
#==============================================================================================================================================================
# This part generates fasta sequences for the chromosome 12 TSS, genes and IGS and returns how many of each display the motifs identified with MEME
# (this is done using MAST). It also scans the entire chromosome 12 for motif occurrences (regardless of them being in peaks or not, using FIMO)
#==============================================================================================================================================================
# The starting files are the three summits files with peaks in TSS, IGS and genes. They need to be expanded 50 bps each way and then converted to multi fasta


		echo "Examining the motif occurrence within chromosome 12 IGS/TSS/genes summits"
		echo "Generating fasta files for each region"
		bedtools slop -i ${prefix}_chr12_inIGS_nozero.bed -g chr12Info.txt -b 50 > ${prefix}_chr12_inIGS_100bp.bed
		bedtools slop -i ${prefix}_chr12_inTSS_nozero.bed -g chr12Info.txt -b 50 > ${prefix}_chr12_inTSS_100bp.bed
		bedtools slop -i ${prefix}_chr12_ingenes_nozero.bed -g chr12Info.txt -b 50 > ${prefix}_chr12_ingenes_100bp.bed

		bedtools getfasta -name -fi $hg19 -bed ${prefix}_chr12_inIGS_100bp.bed -fo ${prefix}_chr12_inIGS_100bp.fasta
		bedtools getfasta -name -fi $hg19 -bed ${prefix}_chr12_inTSS_100bp.bed -fo ${prefix}_chr12_inTSS_100bp.fasta
		bedtools getfasta -name -fi $hg19 -bed ${prefix}_chr12_ingenes_100bp.bed -fo ${prefix}_chr12_ingenes_100bp.fasta

#MAST syntax for each file:
		echo "Searching for MEME motifs in chromosome 12 peaks"
		mast ${prefix}_meme_OUT_FOLDER/meme.txt -hit_list ${prefix}_chr12_inIGS_100bp.fasta -oc ${prefix}_IGS_mast_OUT_FOLDER > ${prefix}_mast_hits.txt

		mast ${prefix}_meme_OUT_FOLDER/meme.txt -hit_list ${prefix}_chr12_inTSS_100bp.fasta -oc ${prefix}_TSS_mast_OUT_FOLDER > ${prefix}_mast_hits.txt

		mast ${prefix}_meme_OUT_FOLDER/meme.txt -hit_list ${prefix}_chr12_ingenes_100bp.fasta -oc ${prefix}_genes_mast_OUT_FOLDER > ${prefix}_mast_hits.txt

#FIMO

		echo "Generating chromosome 12 background file for FIMO"
		fasta-get-markov $chr12 > chr12_bkgrnd.txt

		echo "Using FIMO on Chromosome 12 (whole)"
		fimo --oc ${prefix}_Chr12_all_fimo_OUT_FOLDER --bgfile chr12_bkgrnd.txt ${prefix}_meme_OUT_FOLDER/meme.txt $chr12
		echo "Using FIMO on Chromosome 12 IGS sequences"
		fimo --oc ${prefix}_IGS_fimo_OUT_FOLDER --bgfile chr12_bkgrnd.txt ${prefix}_meme_OUT_FOLDER/meme.txt gencode_ENSG_geneID_chr12_IGS.fasta
		echo "Using FIMO on Chromosome 12 TSS sequences"
		fimo --oc ${prefix}_TSS_fimo_OUT_FOLDER --bgfile chr12_bkgrnd.txt ${prefix}_meme_OUT_FOLDER/meme.txt gencode_ENSG_geneID_chr12_TSS.fasta
		echo "Using FIMO on Chromsome 12 gene encoding sequences"
		fimo --oc ${prefix}_genes_fimo_OUT_FOLDER --bgfile chr12_bkgrnd.txt ${prefix}_meme_OUT_FOLDER/meme.txt gencode_ENSG_geneID_chr12_genes.fasta
		echo "FIMO analysis done!"

		echo "Reformatting FIMO outputs to .bed"

		awk 'NR>1 {printf("%s\t%s\t%s\t%s\t%s\t%s\n", $2, $3, $4, $9, $7, $5)
		}' ${prefix}_Chr12_all_fimo_OUT_FOLDER/fimo.txt > ${prefix}_Chr12_all_fimo_OUT_FOLDER/fimo_chr12.bed

		awk 'NR>1 {printf("%s\t%s\t%s\t%s\t%s\t%s\n", $2, $3, $4, $9, $7, $5)
		}' ${prefix}_IGS_fimo_OUT_FOLDER/fimo.txt > ${prefix}_IGS_fimo_OUT_FOLDER/fimo_IGS.bed

		awk 'NR>1 {printf("%s\t%s\t%s\t%s\t%s\t%s\n", $2, $3, $4, $9, $7, $5)
		}' ${prefix}_TSS_fimo_OUT_FOLDER/fimo.txt > ${prefix}_TSS_fimo_OUT_FOLDER/fimo_TSS.bed

		awk 'NR>1 {printf("%s\t%s\t%s\t%s\t%s\t%s\n", $2, $3, $4, $9, $7, $5)
		}' ${prefix}_genes_fimo_OUT_FOLDER/fimo.txt > ${prefix}_genes_fimo_OUT_FOLDER/fimo_genes.bed

		echo "Looking in JASPAR MEME databases (tomtom)"
		tomtom -eps -m 1 -o ${prefix}_tomtom_OUT ${prefix}_meme_OUT_FOLDER/meme.txt $jaspar_meme


	done | tee -a ${outPATH}logfiles/log_geneLists.txt


#==============================================================================================================================================================
# Unloading of required modules:

module unload bowtie2/2.3.1/
module unload fastqc/0.11.5/
module unload samtools/1.3.1/
module unload BedTools/2.26.0/
