
        ## download reads and convert to FASTQ files
        fasterq-dump             --split-3             -o ZI126             -O /media/inter/mkapun/projects/InvChapter/data/reads             -e 8             -f             -p             SRR203336
        ## compress data
        gzip /media/inter/mkapun/projects/InvChapter/data/reads/ZI126*
        
