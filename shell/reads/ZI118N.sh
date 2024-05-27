
        ## download reads and convert to FASTQ files
        fasterq-dump             --split-3             -o ZI118N             -O /media/inter/mkapun/projects/InvChapter/data/reads             -e 8             -f             -p             SRR654664
        ## compress data
        gzip /media/inter/mkapun/projects/InvChapter/data/reads/ZI118N*
        
