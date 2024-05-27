
        ## download reads and convert to FASTQ files
        fasterq-dump             --split-3             -o ZI167             -O /media/inter/mkapun/projects/InvChapter/data/reads             -e 8             -f             -p             SRR654553
        ## compress data
        gzip /media/inter/mkapun/projects/InvChapter/data/reads/ZI167*
        
