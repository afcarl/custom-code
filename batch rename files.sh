# this command will rename files after the name of upper directory
$ for subdir in *; do cp $subdir/genes.fpkm_tracking ./genes/$subdir.genes; done;

$ for subdir in *; do cp $subdir/accepted_hits.bam.flagstat ./flagstat/$subdir.accepted_hits.bam.flagstat; done;