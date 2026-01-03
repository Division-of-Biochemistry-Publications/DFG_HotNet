#!/bin/bash

WD=$(pwd) #working directory
DD=${WD}/01_Raw_data #defined different sub folders‚
QC1=${WD}/02_FASTQC
TD=${WD}/03_Trimmed_Data
QC2=${WD}/04_QC_Trimmed
MD=${WD}/00_Meta_Data
MAPD=${WD}/05_Mapped_Data
FCD=${WD}/07_FeatureCounts

REF_BASE=${MD}/01_Referenzgenom_Potato_v6.1/Stuberosum/v6.1
STAR_GENOME=${REF_BASE}/STAR_GENOME
REF_FASTA=${REF_BASE}/assembly/Stuberosum_686_v6.1.fa
REF_ANN=${REF_BASE}/annotation/Stuberosum_686_v6.1.gene_exons.gff3
ADAPTERS=${MD}/adapter.fa

threads=11 # adjust threads according to your computational capacity


##############################################################################
# Quality Control with FastQC BEFORE Trimming
##############################################################################

fastqc -o ${QC1} -t ${threads} ${DD}/*.fq.gz 

multiqc \
 --outdir ${QC1}/01_FASTQC_Summary \
 ${QC1} 


##############################################################################
# Trimming (BBDUK for adapter trimming and quality)
##############################################################################

cd ${DD}

for i in *_1.fq.gz; do
  SAMPLE=$(echo ${i} | sed "s/_1\.fq\.gz//")
  
echo ${SAMPLE} "right trimming"

bbduk\
 in=${DD}/${SAMPLE}_1.fq.gz\
 in2=${DD}/${SAMPLE}_2.fq.gz\
 out=${TD}/Clean/${SAMPLE}_1_clean.fq.gz\
 out2=${TD}/Clean/${SAMPLE}_2_clean.fq.gz\
 ref=${ADAPTERS}\
 ktrim=l\
 k=23\
 mink=11\
 hdist=1\
 tpe=t\
 tbo=t\
 t=${threads}


echo ${SAMPLE} "left trimming"

bbduk\
 in=${TD}/Clean/${SAMPLE}_1_clean.fq.gz\
 in2=${TD}/Clean/${SAMPLE}_2_clean.fq.gz\
 out=${TD}/Clean2/${SAMPLE}_1_clean2.fq.gz\
 out2=${TD}/Clean2/${SAMPLE}_2_clean2.fq.gz\
 ref=${ADAPTERS}\
 ktrim=r\
 k=23\
 mink=11\
 hdist=1\
 tpe=t\
 tbo=t\
 t=${threads}


echo ${SAMPLE} "quality trimming"

bbduk\
 in=${TD}/Clean2/${SAMPLE}_1_clean2.fq.gz\
 in2=${TD}/Clean2/${SAMPLE}_2_clean2.fq.gz\
 out=${TD}/${SAMPLE}_1_t.fq.gz\
 out2=${TD}/${SAMPLE}_2_t.fq.gz\
 qtrim=rl\
 trimq=30\
 minlength=35\
 minavgquality=30\
 t=${threads}
 
done


rm -r ${TD}/Clean #remove temporary folders (Clean + Clean2) to save storage space
rm -r ${TD}/Clean2



##############################################################################
# Quality Control AFTER Trimming
##############################################################################

echo "quality control using fastQC on trimmed data"

fastqc -o ${QC2} -t ${threads} ${TD}/*.fq.gz

‚
echo "quality control using multiQC on trimmed data FASTQC"

multiqc \
 --outdir ${QC2}/01_FASTQC_Summary\
 ${QC2}



##############################################################################
# Mapping using STAR
##############################################################################

mkdir -p ${STAR_GENOME}

STAR\
 --runMode genomeGenerate\
 --genomeSAindexNbases=13\
 --runThreadN=${threads}\
 --genomeDir ${STAR_GENOME}\
 --genomeFastaFiles ${REF_FASTA}\
 --sjdbGTFfile ${REF_ANN}\
  --sjdbOverhang=149\
 --sjdbGTFtagExonParentTranscript=ID\
 --sjdbGTFfeatureExon=exon\
 --sjdbGTFtagExonParentGene=Parent

cd ${TD}

for i in *_1_t.fq.gz;do
  SAMPLE=$(echo ${i} | sed "s/_1_t\.fq\.gz//")  
echo ${SAMPLE} "RUN 1"

STAR\
 --runThreadN=${threads}\
 --genomeDir ${STAR_GENOME}\
 --readFilesIn ${TD}/${SAMPLE}_1_t.fq.gz ${TD}/${SAMPLE}_2_t.fq.gz\
 --readFilesCommand zcat\
 --outSAMtype BAM SortedByCoordinate\
 --alignIntronMax=30000\
 --alignMatesGapMax=30000\
 --outReadsUnmapped=Fastx\
 --outFilterMultimapNmax=1000\
 --outFileNamePrefix ${MAPD}/tmp/${SAMPLE}_  

done

for i in *_1_t.fq.gz;do
  SAMPLE=$(echo ${i} | sed "s/_1_t\.fq\.gz//")  
echo ${SAMPLE} "RUN 2"
STAR\
 --runThreadN=${threads}\
 --genomeDir ${STAR_GENOME}
 --readFilesIn ${TD}/${SAMPLE}_1_t.fq.gz ${TD}/${SAMPLE}_2_t.fq.gz\
 --readFilesCommand zcat\
 --outSAMtype BAM SortedByCoordinate\
 --alignIntronMax=30000\
 --alignMatesGapMax=30000\
 --outReadsUnmapped=Fastx\
 --outFilterMultimapNmax=1000\
 --outMultimapperOrder=Random\
 --runRNGseed=1234\
 --outFileNamePrefix ${MAPD}/${SAMPLE}_\
 --sjdbFileChrStartEnd ${MAPD}/tmp/*_SJ.out.tab

samtools index\
 ${MAPD}/${SAMPLE}_Aligned.sortedByCoord.out.bam\
 ${MAPD}/${SAMPLE}_Aligned.sortedByCoord.out.bai

done

rm -r ${MAPD}/tmp #remove temporary folder to save storage space



##############################################################################
# FeatureCounts
##############################################################################

mkdir -p ${FCD}

for feature in gene;do
featureCounts\
 -t gene\
 -g ID\
 -T ${feature}\
 -p --countReadPairs\
 -a ${REF_ANN}\
 -o ${FCD}/07_FeatureCounts_${feature}_UNIQUE.txt\
 ${MAPD}/*_Aligned.sortedByCoord.out.bam
 
done



exit
