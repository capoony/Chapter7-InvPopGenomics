
        ## download reads and convert to FASTQ files
        fasterq-dump             --split-3             -o ZI114N             -O /media/inter/mkapun/projects/InvChapter/data/reads             -e 8             -f             -p             SRR654662
        ## compress data
        gzip /media/inter/mkapun/projects/InvChapter/data/reads/ZI114N*
        
