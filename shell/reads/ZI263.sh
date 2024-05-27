
        ## download reads and convert to FASTQ files
        fasterq-dump             --split-3             -o ZI263             -O /media/inter/mkapun/projects/InvChapter/data/reads             -e 8             -f             -p             SRR210783
        ## compress data
        gzip /media/inter/mkapun/projects/InvChapter/data/reads/ZI263*
        
