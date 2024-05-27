
        ## download reads and convert to FASTQ files
        fasterq-dump             --split-3             -o ZI291             -O /media/inter/mkapun/projects/InvChapter/data/reads             -e 8             -f             -p             SRR202083
        ## compress data
        gzip /media/inter/mkapun/projects/InvChapter/data/reads/ZI291*
        
