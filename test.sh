###用于运行STAR（Spliced Transcripts Alignment to a Reference）的命令行命令，用于将RNA测序数据比对到人类基因组GRCh38.p14的索引上。

STAR --runThreadN 24 --genomeDir /home/zhushunxin/reference/GRCh38.p14/index/star/ --readFilesIn /home/zhushunxin/analysis/rG4-seq/2.trim/SRR14510889_1.barcode.trim.fastq.gz /home/zhushunxin/analysis/rG4-seq/2.trim/SRR14510889_2.barcode.trim.fastq.gz --outFileNamePrefix ./tmp/SRR14510889. --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFilterMismatchNoverReadLmax 0.1 --outFilterMismatchNmax 999 --alignEndsType Extend5pOfReads12 --alignSoftClipAtReferenceEnds No --outSAMtype None --readFilesCommand zcat


STAR --runThreadN 24 --genomeDir /home/zhushunxin/reference/GRCh38.p14/index/star/ --readFilesIn /home/zhushunxin/analysis/rG4-seq/2.trim/SRR14510890_1.barcode.trim.fastq.gz /home/zhushunxin/analysis/rG4-seq/2.trim/SRR14510890_2.barcode.trim.fastq.gz --outFileNamePrefix ./tmp/SRR14510890. --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFilterMismatchNoverReadLmax 0.1 --outFilterMismatchNmax 999 --alignEndsType Extend5pOfReads12 --alignSoftClipAtReferenceEnds No --outSAMtype None --readFilesCommand zcat


STAR --runThreadN 24 --genomeDir /home/zhushunxin/reference/GRCh38.p14/index/star/ --readFilesIn /home/zhushunxin/analysis/rG4-seq/2.trim/SRR14510893_1.barcode.trim.fastq.gz /home/zhushunxin/analysis/rG4-seq/2.trim/SRR14510893_2.barcode.trim.fastq.gz --outFileNamePrefix ./tmp/SRR14510893. --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFilterMismatchNoverReadLmax 0.1 --outFilterMismatchNmax 999 --alignEndsType Extend5pOfReads12 --alignSoftClipAtReferenceEnds No --outSAMtype None --readFilesCommand zcat


STAR --runThreadN 24 --genomeDir /home/zhushunxin/reference/GRCh38.p14/index/star/ --readFilesIn /home/zhushunxin/analysis/rG4-seq/2.trim/SRR14510894_1.barcode.trim.fastq.gz /home/zhushunxin/analysis/rG4-seq/2.trim/SRR14510894_2.barcode.trim.fastq.gz --outFileNamePrefix ./tmp/SRR14510894. --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFilterMismatchNoverReadLmax 0.1 --outFilterMismatchNmax 999 --alignEndsType Extend5pOfReads12 --alignSoftClipAtReferenceEnds No --outSAMtype None --readFilesCommand zcat

cat ./star1/*.SJ.out.tab > pass1.SJ.out.tab

STAR --runThreadN 24 --genomeDir /home/zhushunxin/reference/GRCh38.p14/index/star/ --readFilesIn /home/zhushunxin/analysis/rG4-seq/2.trim/SRR14510889_1.barcode.trim.fastq.gz /home/zhushunxin/analysis/rG4-seq/2.trim/SRR14510889_2.barcode.trim.fastq.gz --sjdbFileChrStartEnd ./pass1.SJ.out.tab --outFileNamePrefix ./SRR14510889. --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFilterMismatchNoverReadLmax 0.1 --outFilterMismatchNmax 999 --alignEndsType Extend5pOfReads12 --alignSoftClipAtReferenceEnds No --outSAMtype BAM Unsorted --outReadsUnmapped Fastx --outBAMcompression 10 --outSAMattrIHstart 0 --outSAMattributes All --readFilesCommand zcat --outSAMattrRGline ID:1 LB:library1 SM:sample1 PL:illumina PU:hiseqX


STAR --runThreadN 24 --genomeDir /home/zhushunxin/reference/GRCh38.p14/index/star/ --readFilesIn /home/zhushunxin/analysis/rG4-seq/2.trim/SRR14510890_1.barcode.trim.fastq.gz /home/zhushunxin/analysis/rG4-seq/2.trim/SRR14510890_2.barcode.trim.fastq.gz --sjdbFileChrStartEnd ./pass1.SJ.out.tab --outFileNamePrefix ./SRR14510890. --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFilterMismatchNoverReadLmax 0.1 --outFilterMismatchNmax 999 --alignEndsType Extend5pOfReads12 --alignSoftClipAtReferenceEnds No --outSAMtype BAM Unsorted --outReadsUnmapped Fastx --outBAMcompression 10 --outSAMattrIHstart 0 --outSAMattributes All --readFilesCommand zcat --outSAMattrRGline ID:1 LB:library1 SM:sample1 PL:illumina PU:hiseqX


STAR --runThreadN 24 --genomeDir /home/zhushunxin/reference/GRCh38.p14/index/star/ --readFilesIn /home/zhushunxin/analysis/rG4-seq/2.trim/SRR14510893_1.barcode.trim.fastq.gz /home/zhushunxin/analysis/rG4-seq/2.trim/SRR14510893_2.barcode.trim.fastq.gz --sjdbFileChrStartEnd ./pass1.SJ.out.tab --outFileNamePrefix ./SRR14510893. --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFilterMismatchNoverReadLmax 0.1 --outFilterMismatchNmax 999 --alignEndsType Extend5pOfReads12 --alignSoftClipAtReferenceEnds No --outSAMtype BAM Unsorted --outReadsUnmapped Fastx --outBAMcompression 10 --outSAMattrIHstart 0 --outSAMattributes All --readFilesCommand zcat --outSAMattrRGline ID:1 LB:library1 SM:sample1 PL:illumina PU:hiseqX


STAR --runThreadN 24 --genomeDir /home/zhushunxin/reference/GRCh38.p14/index/star/ --readFilesIn /home/zhushunxin/analysis/rG4-seq/2.trim/SRR14510894_1.barcode.trim.fastq.gz /home/zhushunxin/analysis/rG4-seq/2.trim/SRR14510894_2.barcode.trim.fastq.gz --sjdbFileChrStartEnd ./pass1.SJ.out.tab --outFileNamePrefix ./SRR14510894. --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFilterMismatchNoverReadLmax 0.1 --outFilterMismatchNmax 999 --alignEndsType Extend5pOfReads12 --alignSoftClipAtReferenceEnds No --outSAMtype BAM Unsorted --outReadsUnmapped Fastx --outBAMcompression 10 --outSAMattrIHstart 0 --outSAMattributes All --readFilesCommand zcat --outSAMattrRGline ID:1 LB:library1 SM:sample1 PL:illumina PU:hiseqX

#!/bin/bash
#PBS -N rG4_sort
#PBS -l ncpus=12
#PBS -l mem=64G
#PBS -j oe
#PBS -o ./pbs_info

###此脚本主要用于处理BAM文件，包括排序、索引和过滤，以便进行后续的分析或可视化操作。

source activate structure
cd ~/analysis/rG4-seq/3.star_bam

for bam in `ls ./*.Aligned.out.bam`
do
name=`basename ${bam}`
#echo "
samtools sort ${bam} -o ./sort_bam/${name%%.*}.Aligned.sorted.bam -@ 24
samtools index ./sort_bam/${name%%.*}.Aligned.sorted.bam -o ./sort_bam/${name%%.*}.Aligned.sorted.bam.bai -@ 24
#"
done
for bam2 in `ls ./sort_bam/*.Aligned.sorted.bam`
do
name2=`basename ${bam2}`
#echo"
samtools view -q 255 -F3852 -bh ./sort_bam/${name2%%.*}.Aligned.sorted.bam > ./sort_bam/${name2%%.*}.Aligned.sorted.uniquely_aligned.bam
#"
done

#!/bin/bash
#PBS -N rG4_align
#PBS -l ncpus=12
#PBS -l mem=64G
#PBS -j oe
#PBS -o ./pbs_info

source activate structure
cd ~/analysis/rG4-seq/3.star_bam

###自动化地对多个FASTQ文件进行STAR比对，并将比对的结果存储在tmp目录中。

mkdir tmp
GENOMEDIR=~/reference/GRCh38.p14/index/star/
for seq in `ls ~/analysis/rG4-seq/2.trim/*_1.barcode.trim.fastq.gz`
do
#echo "
STAR --runThreadN 24 --genomeDir ${GENOMEDIR} \
--readFilesIn ${seq%_*}_1.barcode.trim.fastq.gz ${seq%_*}_2.barcode.trim.fastq.gz \
--outFileNamePrefix ./tmp/$(basename ${seq%_*}). \
--outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
--outFilterMismatchNoverReadLmax 0.1 --outFilterMismatchNmax 999 \
--alignEndsType Extend5pOfReads12 --alignSoftClipAtReferenceEnds No \
--outSAMtype None --readFilesCommand zcat
#"
done

###对一组FASTQ文件执行STAR比对，并在比对的过程中引用了之前生成的pass1.SJ.out.tab文件。

cat ./star1/*.SJ.out.tab > pass1.SJ.out.tab
for seq2 in `ls ~/analysis/rG4-seq/2.trim/*_1.barcode.trim.fastq.gz`
do
#echo "
STAR --runThreadN 24 --genomeDir ${GENOMEDIR} \
--readFilesIn ${seq2%_*}_1.barcode.trim.fastq.gz ${seq2%_*}_2.barcode.trim.fastq.gz \
--sjdbFileChrStartEnd ./pass1.SJ.out.tab \
--outFileNamePrefix ./$(basename ${seq2%_*}). \
--outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
--outFilterMismatchNoverReadLmax 0.1 --outFilterMismatchNmax 999 \
--alignEndsType Extend5pOfReads12 --alignSoftClipAtReferenceEnds No \
--outSAMtype BAM Unsorted --outReadsUnmapped Fastx --outBAMcompression 10 \
--outSAMattrIHstart 0 --outSAMattributes All --readFilesCommand zcat \
--outSAMattrRGline ID:1 LB:library1 SM:sample1 PL:illumina PU:hiseqX
#"
done

###每一行mv命令都用于将文件从旧的文件名更改为新的文件名。

#!/bin/bash
mv SRR14510889.Aligned.sorted.uniquely_aligned.dedup.bam K-rep2.sorted.uniquely_aligned.dedup.bam
mv SRR14510890.Aligned.sorted.uniquely_aligned.dedup.bam K-rep1.sorted.uniquely_aligned.dedup.bam
mv SRR14510893.Aligned.sorted.uniquely_aligned.dedup.bam Li-rep2.sorted.uniquely_aligned.dedup.bam
mv SRR14510894.Aligned.sorted.uniquely_aligned.dedup.bam Li-rep1.sorted.uniquely_aligned.dedup.bam

#!/bin/bash
#PBS -N rG4_dedup
#PBS -l ncpus=5
#PBS -l mem=40G
#PBS -j oe
#PBS -o ./pbs_info

###在给定的BAM文件上运行umi_tools dedup命令，以去除PCR重复的reads，并将去重后的结果保存为新的BAM文件。

source activate Seq
cd ~/analysis/rG4-seq/4.dedup 

for bam in `ls ~/analysis/rG4-seq/3.star_bam/sort_bam/*.Aligned.sorted.uniquely_aligned.bam`
do
id=$(basename ${bam})
umi_tools dedup \
--stdin=${bam} \
--stdout=./$(basename ${id%%.*}).Aligned.sorted.uniquely_aligned.dedup.bam \
--chimeric-pairs=discard --unpaired-reads=discard --paired --ignore-tlen \
--buffer-whole-contig --mapping-quality 255
done

#!/bin/bash
#PBS -N rG4_bamcompare
#PBS -l ncpus=12
#PBS -l mem=40G
#PBS -j oe
#PBS -o ./bamcompare_info

###行两个样本之间的比对，计算log2比对，并将结果以bigwig格式保存。

source activate chip-seq
cd ~/analysis/rG4-seq/6.bigwig

bamCompare -b1 /home/zhushunxin/analysis/rG4-seq/4.dedup/K-rep1.sorted.uniquely_aligned.dedup.bam -b2 /home/zhushunxin/analysis/rG4-seq/4.dedup/Li-rep1.sorted.uniquely_aligned.dedup.bam -o rep1_log2.bw -of bigwig -bs 10 -p 24 --effectiveGenomeSize 2913022398 --normalizeUsing BPM --operation log2 --skipZeroOverZero --scaleFactorsMethod None

bamCompare -b1 /home/zhushunxin/analysis/rG4-seq/4.dedup/K-rep2.sorted.uniquely_aligned.dedup.bam -b2 /home/zhushunxin/analysis/rG4-seq/4.dedup/Li-rep2.sorted.uniquely_aligned.dedup.bam -o rep2_log2.bw -of bigwig -bs 10 -p 24 --effectiveGenomeSize 2913022398 --normalizeUsing BPM --operation log2 --skipZeroOverZero --scaleFactorsMethod None
#!/bin/bash
#PBS -N rG4_bamcoverage
#PBS -l ncpus=12
#PBS -l mem=40G
#PBS -j oe
#PBS -o ./bamcoverage_info

source activate chip-seq
cd ~/analysis/rG4-seq/6.bigwig

bamCoverage -b /home/zhushunxin/analysis/rG4-seq/4.dedup/K-rep1.sorted.uniquely_aligned.dedup.bam -o ./K-rep1.bw --effectiveGenomeSize 2913022398 -p 24 -bs 10 --normalizeUsing RPGC -of bigwig
bamCoverage -b /home/zhushunxin/analysis/rG4-seq/4.dedup/K-rep2.sorted.uniquely_aligned.dedup.bam -o ./K-rep2.bw --effectiveGenomeSize 2913022398 -p 24 -bs 10 --normalizeUsing RPGC -of bigwig
bamCoverage -b /home/zhushunxin/analysis/rG4-seq/4.dedup/Li-rep1.sorted.uniquely_aligned.dedup.bam -o ./Li-rep1.bw --effectiveGenomeSize 2913022398 -p 24 -bs 10 --normalizeUsing RPGC -of bigwig
bamCoverage -b /home/zhushunxin/analysis/rG4-seq/4.dedup/Li-rep2.sorted.uniquely_aligned.dedup.bam -o ./Li-rep2.bw --effectiveGenomeSize 2913022398 -p 24 -bs 10 --normalizeUsing RPGC -of bigwig

#!/bin/bash
#PBS -N rG4_deeptools
#PBS -l ncpus=12
#PBS -l mem=40G
#PBS -j oe
#PBS -o ./deeptools_info


# for deeptools metagene
# mkdir matrix && pdf

source activate chip-seq
cd ~/analysis/rG4-seq/7.deeptools

###生成基因表达矩阵，并使用该矩阵绘制基因表达剖面图，以便分析基因在不同样本中的表达情况。

prefix=rG4_point

computeMatrix reference-point \
-a 500 \
-b 500 \
--referencePoint TSS \
-S /home/zhushunxin/analysis/rG4-seq/6.bigwig/K-rep1.bw /home/zhushunxin/analysis/rG4-seq/6.bigwig/Li-rep1.bw \
-R /home/zhushunxin/analysis/rG4-seq/5.rG4_seeker/rG4.bed \
-p 24 \
-o ./matrix/${prefix}.mtx

plotProfile -m ./matrix/${prefix}.mtx \
-o ./pdf/${prefix}.pdf \
--dpi 300 \
--plotHeight 7 --plotWidth 10 --plotFileFormat pdf --perGroup \
--refPointLabel "rG4 peaks" \
-T rG4 \
--colors red blue 

#!/bin/bash
#PBS -N rG4_deeptools
#PBS -l ncpus=12
#PBS -l mem=40G
#PBS -j oe
#PBS -o ./deeptools_info


# for deeptools metagene
# mkdir matrix && pdf

source activate chip-seq
cd ~/analysis/rG4-seq/7.deeptools

###生成基因表达矩阵，以及绘制基因表达剖面图，以比较基本基因和lncRNA基因的表达情况。

prefix=basic_lnc

computeMatrix scale-regions \
--metagene \
-m 1000 \
-a 0 \
-b 0 \
-p 24 \
-S /home/zhushunxin/analysis/rG4-seq/6.bigwig/K-rep1.bw /home/zhushunxin/analysis/rG4-seq/6.bigwig/Li-rep1.bw \
-R /home/zhushunxin/reference/GRCh38.p14/gtf/gencode.v44.basic.annotation.gtf /home/zhushunxin/reference/GRCh38.p14/gtf/gencode.v44.long_noncoding_RNAs.gtf \
--transcriptID transcript \
--exonID exon \
--transcript_id_designator transcript_id \
-o ./matrix/${prefix}.mtx

plotProfile -m ./matrix/${prefix}.mtx \
-o ./pdf/${prefix}.pdf \
--dpi 300 \
--plotHeight 7 --plotWidth 10 --plotFileFormat pdf --perGroup \
-T "basic genes && lncRNA genes" \
--yMin -0.02 \
--yMax 0.02 \
--colors red blue
