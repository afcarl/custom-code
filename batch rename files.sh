# this command will rename files after the name of upper directory
$ for subdir in *; do cp $subdir/genes.fpkm_tracking ./genes/$subdir.genes; done;

$ for subdir in *; do cp $subdir/accepted_hits.bam.flagstat ./flagstat/$subdir.accepted_hits.bam.flagstat; done;

$ for subdir in *; do cp $subdir/*paired_end.fastq.gz.sequence_info ./filenames/$subdir.txt; done;