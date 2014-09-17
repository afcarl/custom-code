rm findings.txt
find *.result > findings.txt
sed -i 's/.result//g' findings.txt
Rscript missing.R