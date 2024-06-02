################# ANALYSIS PIPELINE ##################

### define working directory
WD=/media/inter/mkapun/projects/InvChapter

## (1) install dependencies
sh ${WD}/shell/dependencies

## (2) Get information of individual sequencing data and isolate samples with known inversion status
mkdir ${WD}/data
cd ${WD}/data

### download Excel table
wget http://johnpool.net/TableS1_individuals.xls

### process table and generate input files for downstream analyses
Rscript ${WD}/scripts/ReadXLS.r ${WD}

### Define arrays with the inverions names, chromosome, start and end breakpoints; These data will be reused in the whole pipleine for the sequential analysis and visulaization of both focal inversions
DATA=("IN2Lt" "IN3RP")
Chrom=("2L" "3R")
Start=(2225744 16432209)
End=(13154180 24744010)

## (3) Get read data from SRA
mkdir ${WD}/data/reads
mkdir ${WD}/shell/reads
conda activate sra-tools

### loop over both inversions
for index in ${!DATA[@]}; do
    INVERSION=${DATA[index]}
    while
        IFS=',' read -r ID SRR Inv
    do
        if [[ -f ${WD}/data/reads/${ID}_1.fastq.gz ]]; then
            continue
        fi

        echo """
        ## download reads and convert to FASTQ files
        fasterq-dump \
            --split-3 \
            -o ${ID} \
            -O ${WD}/data/reads \
            -e 8 \
            -f \
            -p \
            ${SRR}
        ## compress data
        gzip ${WD}/data/reads/${ID}*
        """ >${WD}/shell/reads/${ID}.sh
        sh ${WD}/shell/reads/${ID}.sh
    done <${WD}/data/${INVERSION}.txt
done

## (4) map reads
### obtain Drosophila reference from FlyBase
cd ${WD}/data
wget -O dmel-6.57.fa.gz http://ftp.flybase.net/genomes/Drosophila_melanogaster/current/fasta/dmel-all-chromosome-r6.57.fasta.gz

### index the reference
conda activate bwa-mem2
bwa-mem2 index dmel-6.57.fa.gz
gunzip -c dmel-6.57.fa.gz >dmel-6.57.fa
samtools faidx dmel-6.57.fa
samtools dict dmel-6.57.fa >dmel-6.57.dict
conda deactivate

### trim & map & sort & remove duplicates & realign around indels
for index in ${!DATA[@]}; do
    INVERSION=${DATA[index]}
    while
        IFS=',' read -r ID SRR Inv
    do
        if [[ ${ID} == "Stock ID" || -f ${WD}/mapping/${ID}_RG.bam ]]; then
            continue
        fi

        sh ${WD}/shell/mapping.sh \
            ${WD}/data/reads/${ID}_1.fastq.gz \
            ${WD}/data/reads/${ID}_2.fastq.gz \
            ${ID} \
            ${WD}/mapping \
            ${WD}/data/dmel-6.57 \
            100 \
            ${WD}/scripts/gatk/GenomeAnalysisTK.jar
    done <${WD}/data/${INVERSION}.txt
done

## (5) SNP calling using freebayes with 100 threads
for index in ${!DATA[@]}; do

    INVERSION=${DATA[index]}
    while
        IFS=',' read -r ID SRR Inv
    do
        if [[ ${ID} == "Stock ID" ]]; then
            continue
        fi

        mkdir -p ${WD}/results/SNPs_${INVERSION}

        ### store the PATHs to all BAM files in a text, which will be used as the input for FreeBayes
        echo ${WD}/mapping/${ID}_RG.bam >>${WD}/mapping/BAMlist_${INVERSION}.txt

    done <${WD}/data/${INVERSION}.txt

    conda activate freebayes

    ### run FreeBayes in parallel by splitting the reference genome in chuncks of 100,000bps and and use GNU parallel for multithreading. I am using 100 threads. Please adjust to your system.
    freebayes-parallel \
        <(fasta_generate_regions.py \
            ${WD}/data/dmel-6.57.fa.fai \
            100000) \
        100 \
        -f ${WD}/data/dmel-6.57.fa \
        -L ${WD}/mapping/BAMlist_${INVERSION}.txt \
        --ploidy 1 |
        gzip >${WD}/results/SNPs_${INVERSION}/SNPs_${INVERSION}.vcf.gz
    conda deactivate
done

## (6) calculate FST between karyotypes

### make input files for STD and INV samples
for i in IN2Lt; do
    mkdir ${WD}/data/${INVERSION}
    output_dir=${WD}/data/${INVERSION}
    ### split file with sample IDs based on Inversions status
    awk -F',' '
    {
        filename = $3 ".csv"
        filepath = "'$output_dir'/" filename
        if (filename == ".csv") next
        print $1 >> filepath
    }
    ' ${WD}/data/${INVERSION}.txt

    # ### filter VCF for biallelic SNPs
    conda activate vcftools

    vcftools --gzvcf ${WD}/results/SNPs_${INVERSION}/SNPs_${INVERSION}.vcf.gz \
        --min-alleles 2 \
        --max-alleles 2 \
        --remove-indels \
        --recode \
        --out ${WD}/results/SNPs_${INVERSION}/SNPs_${INVERSION}

    gzip ${WD}/results/SNPs_${INVERSION}/SNPs_${INVERSION}.recode.vcf

    ## convert haploid VCF to diploid
    python ${WD}/scripts/hap2dip.py \
        --input ${WD}/results/SNPs_${INVERSION}/SNPs_${INVERSION}.recode.vcf.gz \
        --output ${WD}/results/SNPs_${INVERSION}/SNPs_${INVERSION}.recode_dip.vcf.gz

    for karyo in INV ST; do

        ### calculate PI
        vcftools --gzvcf ${WD}/results/SNPs_${INVERSION}/SNPs_${INVERSION}.recode_dip.vcf.gz \
            --keep ${WD}/data/${INVERSION}/${karyo}.csv \
            --window-pi 200000 \
            --out ${WD}/results/SNPs_${INVERSION}/${INVERSION}_${karyo}_pi
    done

    ## combine pi of INV and ST chromosomes
    awk 'NR ==1 {print $0"\tType"}' ${WD}/results/SNPs_${INVERSION}/${INVERSION}_INV_pi.windowed.pi >${WD}/results/SNPs_${INVERSION}/${INVERSION}_pi.tsv
    awk 'NR>1  {print $0"\tINV"}' ${WD}/results/SNPs_${INVERSION}/${INVERSION}_INV_pi.windowed.pi >>${WD}/results/SNPs_${INVERSION}/${INVERSION}_pi.tsv
    awk 'NR>1  {print $0"\tST"}' ${WD}/results/SNPs_${INVERSION}/${INVERSION}_ST_pi.windowed.pi >>${WD}/results/SNPs_${INVERSION}/${INVERSION}_pi.tsv

done

### plot PI as Manhattan Plots
for index in ${!DATA[@]}; do

    INVERSION=${DATA[index]}
    St=${Start[index]}
    En=${End[index]}
    Ch=${Chrom[index]}

    Rscript ${WD}/scripts/Plot_pi.r \
        ${INVERSION} \
        ${Ch} \
        ${St} \
        ${En} \
        ${WD}

done

## (6) calculate FST between karyotypes

### make input files for STD and INV samples
for index in ${!DATA[@]}; do
    INVERSION=${DATA[index]}

    conda activate vcftools

    ## calculate FST per SNP
    vcftools --gzvcf ${WD}/results/SNPs_${INVERSION}/SNPs_${INVERSION}.recode_dip.vcf.gz \
        --weir-fst-pop ${WD}/data/${INVERSION}/INV.csv \
        --weir-fst-pop ${WD}/data/${INVERSION}/ST.csv \
        --out ${WD}/results/SNPs_${INVERSION}/${INVERSION}.fst

    ## calculate FST
    vcftools --gzvcf ${WD}/results/SNPs_${INVERSION}/SNPs_${INVERSION}.recode_dip.vcf.gz \
        --weir-fst-pop ${WD}/data/${INVERSION}/INV.csv \
        --weir-fst-pop ${WD}/data/${INVERSION}/ST.csv \
        --fst-window-size 100000 \
        --out ${WD}/results/SNPs_${INVERSION}/${INVERSION}_window.fst

    conda deactivate
done

for index in ${!DATA[@]}; do

    INVERSION=${DATA[index]}
    St=${Start[index]}
    En=${End[index]}
    Ch=${Chrom[index]}

    ### plot FST as Manhattan Plots
    Rscript ${WD}/scripts/Plot_fst.r \
        ${INVERSION} \
        ${Ch} \
        ${St} \
        ${En} \
        ${WD}

done

## (6) obtain diagnostic SNPs for inversion

for index in ${!DATA[@]}; do

    INVERSION=${DATA[index]}
    St=${Start[index]}
    En=${End[index]}
    Ch=${Chrom[index]}

    BP="${Ch},${St},${En}"
    gunzip -c ${WD}/results/SNPs_${INVERSION}/SNPs_${INVERSION}.recode.vcf.gz |
        awk -v Ch=${Ch} '$1~/^#/|| $1 == Ch' |
        python ${WD}/scripts/DiagnosticSNPs.py \
            --range 200000 \
            --breakpoints ${BP} \
            --input - \
            --output ${WD}/results/SNPs_${INVERSION}/${INVERSION} \
            --MinCov 10 \
            --Variant ${WD}/data/${INVERSION}.txt
done
## (7) estimate inversion frequency in PoolSeq data

###  download sripts
cd ${WD}/scripts
wget https://raw.githubusercontent.com/DEST-bio/DESTv2_data_paper/main/16.Inversions/scripts/VCF2sync.py
wget https://raw.githubusercontent.com/DEST-bio/DESTv2_data_paper/main/16.Inversions/scripts/inversion-freqs.py
wget https://raw.githubusercontent.com/DEST-bio/DESTv2_data_paper/main/16.Inversions/scripts/overlap_in_SNPs.py
cp /media/inter/mkapun/projects/DESTv2_data_paper/16.Inversions/scripts/overlap_in_SNPs.py .
### download VCF file and metadata for DEST dataset
cd ${WD}/data
wget -O DEST.vcf.gz http://berglandlab.uvadcos.io/vcf/dest.all.PoolSNP.001.50.3May2024.ann.vcf.gz
wget -O meta.csv https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_3May2024.csv

### convert VCF to SYNC file format
conda activate parallel
gunzip -c ${WD}/data/DEST.vcf.gz |
    parallel \
        --jobs 200 \
        --pipe \
        -k \
        --cat python3 ${WD}/scripts/VCF2sync.py \
        --input {} |
    gzip >${WD}/data/DEST.sync.gz

### Get positions at inversion specific marker SNPs
for index in ${!DATA[@]}; do

    INVERSION=${DATA[index]}
    gunzip -c ${WD}/data/DEST.sync.gz |
        parallel \
            --pipe \
            --jobs 20 \
            -k \
            --cat python3 ${WD}/scripts/overlap_in_SNPs.py \
            --source ${WD}/results/SNPs_${INVERSION}/${INVERSION}_diag.txt \
            --target {} \
            >${WD}/data/DEST_${INVERSION}.sync
done

### convert diagnostic SNP file to match prerequisites for inv script
for index in ${!DATA[@]}; do

    INVERSION=${DATA[index]}
    cut -f1-3 ${WD}/results/SNPs_${INVERSION}/${INVERSION}_diag.txt |
        awk -v INV=${INVERSION} 'NR>1{print INV"\t"$0}' \
            >${WD}/results/SNPs_${INVERSION}/${INVERSION}_diag.markers
done

NAMES=$(gunzip -c ${WD}/data/DEST.vcf.gz | head -150 | awk '/^#C/' | cut -f10- | tr '\t' ',')

# Calculate average frequencies for marker SNPs
for index in ${!DATA[@]}; do

    INVERSION=${DATA[index]}
    Ch=${Chrom[index]}

    python3 ${WD}/scripts/inversion-freqs.py \
        ${WD}/results/SNPs_${INVERSION}/${INVERSION}_diag.markers \
        ${WD}/data/DEST_${INVERSION}.sync \
        $NAMES \
        >${WD}/results/SNPs_${INVERSION}/${INVERSION}.af

    gunzip -c ${WD}/data/DEST.vcf.gz |
        awk -v Ch=${Ch} '$1~/^#/|| $1 == Ch' |
        python3 ${WD}/scripts/AFbyAllele.py \
            --input - \
            --diag ${WD}/results/SNPs_${INVERSION}/${INVERSION}_diag.txt \
            >${WD}/results/SNPs_${INVERSION}/${INVERSION}_pos.af
done

### generate plots for each population
for index in ${!DATA[@]}; do

    INVERSION=${DATA[index]}
    Rscript ${WD}/scripts/Plot_InvMarker.r \
        ${INVERSION} \
        ${WD}
done

## (8) now subset to European and North American datasets and calculate AFs

### Split metadata by continent

sed -i "s/'//g" ${WD}/data/meta.csv

awk -F "," '$6 =="Europe" {print $1}' ${WD}/data/meta.csv >${WD}/data/Europe.ids
awk -F "," '$6 =="North_America" {print $1}' ${WD}/data/meta.csv >${WD}/data/NorthAmerica.ids
awk -F "," '$(NF-7) !="Pass" || $(NF-9)<15 {print $1"\t"$(NF-7)"\t"$(NF-9)}' ${WD}/data/meta.csv >${WD}/data/REMOVE.ids

### subset the VCF file to only (1) contain only European data (2) remove problematic populations (based on DEST recommendations), remove (3) populations with < 15-fold average read depth, (4) only retain bilallic SNPs, (5) subsample to 50,000 randomly drawn genome-wide SNPs and (6) convert the allele counts to frequencies and weights (read-depths).

mkdir ${WD}/results/SNPs

conda activate vcftools
pigz -dc ${WD}/data/DEST.vcf.gz |
    awk '$0~/^\#/ || length($5)==1' |
    vcftools --vcf - \
        --keep ${WD}/data/Europe.ids \
        --remove ${WD}/data/REMOVE.ids \
        --recode \
        --stdout |
    grep -v "\./\." |
    python ${WD}/scripts/SubsampleVCF.py \
        --input - \
        --snps 50000 |
    python ${WD}/scripts/vcf2af.py \
        --input - \
        --output ${WD}/results/SNPs/Europe

### do the same for North American samples
pigz -dc ${WD}/data/DEST.vcf.gz |
    awk '$0~/^\#/ || length($5)==1' |
    vcftools --vcf - \
        --keep ${WD}/data/NorthAmerica.ids \
        --remove ${WD}/data/REMOVE.ids \
        --recode \
        --stdout |
    grep -v "\./\." |
    python ${WD}/scripts/SubsampleVCF.py \
        --input - \
        --snps 50000 |
    python ${WD}/scripts/vcf2af.py \
        --input - \
        --output ${WD}/results/SNPs/NorthAmerica

## (9) calculate SNP-wise logistic regressions testing for associations between SNP allele frequencies and inversion frequencies to test for linkage between SNPs and the inversion for Europe and North America

for index in ${!DATA[@]}; do

    INVERSION=${DATA[index]}
    St=${Start[index]}
    En=${End[index]}
    Ch=${Chrom[index]}

    Rscript ${WD}/scripts/PlotInvLD.r \
        ${INVERSION} \
        ${Ch} \
        ${St} \
        ${En} \
        ${WD}

done

## (10) The influence of Inversions on population structure

### use PCA to test for patterns inside and outside the genomic region spanned by an inversion
for index in ${!DATA[@]}; do

    INVERSION=${DATA[index]}
    St=${Start[index]}
    En=${End[index]}
    Ch=${Chrom[index]}

    Rscript ${WD}/scripts/PCA_Inv.r \
        ${INVERSION} \
        ${Ch} \
        ${St} \
        ${En} \
        ${WD}

done

### does the Inv Frequency influence the PCA results?
for index in ${!DATA[@]}; do

    INVERSION=${DATA[index]}
    Rscript ${WD}/scripts/Plot_PCAInvFreq.r \
        ${INVERSION} \
        ${WD}
done

## (11) test for clinality of inversion frequency
for index in ${!DATA[@]}; do

    INVERSION=${DATA[index]}
    Rscript ${WD}/scripts/Plot_Clinality.r \
        ${INVERSION} \
        ${WD}
done

### Test if clinality due to demography or potentially adaptive
for index in ${!DATA[@]}; do

    INVERSION=${DATA[index]}
    St=${Start[index]}
    En=${End[index]}
    Ch=${Chrom[index]}

    Rscript ${WD}/scripts/LFMM.r \
        ${INVERSION} \
        ${Ch} \
        ${St} \
        ${En} \
        ${WD}

done

## copy figures to output folder
mkdir /media/inter/mkapun/projects/InvChapter/output

cp /media/inter/mkapun/projects/InvChapter/results/SNPs*/*.png /media/inter/mkapun/projects/InvChapter/output
cp /media/inter/mkapun/projects/InvChapter/results/SNPs_*/LFMM_*/*.png /media/inter/mkapun/projects/InvChapter/output
cp /media/inter/mkapun/projects/InvChapter/results/SNPs_*/LDwithSNPs/*.png /media/inter/mkapun/projects/InvChapter/output

pandoc -f markdown \
    -t latex \
    -o /media/inter/mkapun/projects/InvChapter/README.pdf \
    /media/inter/mkapun/projects/InvChapter/README.md

pandoc -f markdown \
    -t docx \
    -o /media/inter/mkapun/projects/InvChapter/README.docx \
    /media/inter/mkapun/projects/InvChapter/README.md
