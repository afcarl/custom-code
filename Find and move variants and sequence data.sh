find varicella/ -name '*.sorted.markdup.realigned.recal.bai' -exec mv {} bams/ \;
find varicella/ -name '*.sorted.markdup.realigned.recal.bam' -exec mv {} bams/ \;
find varicella/ -name '*.final_variants.vcf.gz*' -exec cp {} variants/ \;
