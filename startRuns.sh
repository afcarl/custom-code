source /com/extra/R/3.0.1/load.sh
rm x*
split -l 16 commandfile.txt
find /home/michal/molpros/simstudy3/2way -name "x*" >> torun.tmp

NUM_THREADS=1
COMMANDFILE=torun.sh
rm -f ${COMMANDFILE}
while read LIST_FILE
do
echo "qx -p --nodes=1 -c 16 -t 00:59:59 --dir=/home/michal/molpros/simstudy3/2way -i=${LIST_FILE} -i=/home/michal/molpros/simstudy3/2way/dfgTrain_static -i=/home/michal/molpros/simstudy3/2way/dfgEval_static -i=/home/michal/molpros/simstudy3/2way/essentials_SIM.RData -i=/home/michal/molpros/simstudy3/2way/sim_2way_analysis.R -o=*.result dispatch -r ${LIST_FILE}" >> ${COMMANDFILE}
done < torun.tmp
chmod +x torun.sh
./torun.sh
