read1=$1
read2=$2
sample=$3
output=$4
genome=$5
threads=$6
GATK=$7

mkdir -p $output

#activate conda
eval "$(conda shell.bash hook)"

## trim with cutadapt
conda activate cutadapt
cutadapt \
    -q 18 \
    --minimum-length 25 \
    -o $output/${sample}.trimmed1.fq.gz \
    -p $output/${sample}.trimmed2.fq.gz \
    -b ACACTCTTTCCCTACACGACGCTCTTCCGATC \
    -B CAAGCAGAAGACGGCATACGAGAT \
    -O 15 \
    -n 3 \
    --cores=$threads \
    $read1 $read2

conda deactivate

## map & sort
conda activate bwa-mem2
bwa-mem2 mem \
    -t $threads \
    -M \
    -R "@RG\tID:$sample\tSM:$sample\tPL:illumina\tLB:lib1" \
    ${genome}.fa.gz \
    $output/${sample}.trimmed1.fq.gz \
    $output/${sample}.trimmed2.fq.gz |
    samtools view -@ $threads -Sbh -q 20 -F 0x100 - |
    samtools sort -@ $threads - >$output/${sample}.sorted_merged.bam
conda deactivate

## deduplicate
conda activate picard
picard MarkDuplicates \
    REMOVE_DUPLICATES=true \
    I=$output/${sample}.sorted_merged.bam \
    O=$output/${sample}.dedup.bam \
    M=$output/${sample}.mark_duplicates_report.txt \
    VALIDATION_STRINGENCY=SILENT
conda deactivate

## realign around indels
conda activate bwa-mem2
samtools index $output/${sample}.dedup.bam

java -jar $GATK -T RealignerTargetCreator \
    -nt $threads \
    -R ${genome}.fa \
    -I $output/${sample}.dedup.bam \
    -o $output/${sample}.intervals

java -jar $GATK \
    -T IndelRealigner \
    -R ${genome}.fa \
    -I $output/${sample}.dedup.bam \
    -targetIntervals $output/${sample}.intervals \
    -o $output/${sample}.RG.bam

conda deactivate
