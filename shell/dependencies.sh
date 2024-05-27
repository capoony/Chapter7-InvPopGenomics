#activate conda
eval "$(conda shell.bash hook)"

conda create -y --name sra-tools -c bioconda sra-tools

## install CutAdapt
conda create -y --name cutadapt -c "bioconda/label/cf201901" cutadapt

conda create -y --name bwa-mem2 -c bioconda bwa-mem2
conda activate bwa-mem2
conda install -y -c conda-forge -c bioconda samtools
conda deactivate

conda create -y --name picard -c bioconda picard

cd ~/InvChapter/scripts
wget https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2 &&
    tar -xvf GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2 &&
    mv GenomeAnalysisTK-3.8-1-0-gf15c1c3ef gatk &&
    rm GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2

conda create -y --name freebayes -c bioconda freebayes=1.3.6

conda create -y --name vcftools -c bioconda vcftools

conda create -y --name parallel -c bioconda parallel
