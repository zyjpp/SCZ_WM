##------------------------------
sigfile=/FDR/list.txt
line_count=$(awk 'END{print NR}' "$sigfile")

for i in $(seq 1 $line_count);do
IDP=`cat $sigfile|awk 'NR=='$i' {print $1}'`
disease=`cat $sigfile|awk 'NR=='$i' {print $2}'`
name=$(echo "$disease" | awk -F'_' '{print $1}')

## change config.txt and runme.m 
new_traitfolder="/FDR/traitfolder/"
new_traitfile1=$IDP".mat"
new_traitname1=$IDP
new_traitfiles=$disease
new_traitnames=$disease
new_outputdir="/FDR/results/"$IDP"_"$name
new_reffile=/FDR/ref/ref9545380_1kgPhase3eur_LDr2p1.mat
new_randprune_n=500
new_stattype="conjfdr"
new_fdrthresh="0.05"
new_exclude="[6 25119106 33854733; 8 7200000 12500000]"
####################
cp -f /pleiofdr-master/config_default.txt "/FDR/results/config_"$IDP"_"$name"_cjfdr.txt"
cp -f /pleiofdr-master/runme.m "/FDR/results/runme_"$IDP"_"$name"_cjfdr.m" 
configfile="/FDR/results/config_"$IDP"_"$name"_cjfdr.txt"

sed -i "s|reffile=.*|reffile=$new_reffile|g" $configfile
sed -i "s|traitfolder=.*|traitfolder=$new_traitfolder|g" $configfile
sed -i "s|traitfile1=.*|traitfile1=$new_traitfile1|g" $configfile
sed -i "s|traitname1=.*|traitname1=$new_traitname1|g" $configfile
sed -i "s|traitfiles=.*|traitfiles={'${new_traitfiles[@]}'}|g" $configfile
sed -i "s|traitnames=.*|traitnames={'${new_traitnames[@]}'}|g" $configfile
sed -i "s|outputdir=.*|outputdir=$new_outputdir|g" $configfile
sed -i "s|stattype=.*|stattype=$new_stattype|g" $configfile
sed -i "s|fdrthresh=.*|fdrthresh=$new_fdrthresh|g" $configfile
sed -i "s|randprune_n=.*|randprune_n=$new_randprune_n|g" $configfile
sed -i "s|exclude_chr_pos=.*|exclude_chr_pos=$new_exclude|g" $configfile

newconfig="config_"$IDP"_"$name"_cjfdr.txt"
sed -i "4s|config=.*|config='$newconfig';end|g" "$resultfile/runme_"$IDP"_"$name"_cjfdr.m"
mv "/FDR/results/runme_"$IDP"_"$name"_cjfdr.m" /pleiofdr-master/ -f

mv /FDR/results/ /pleiofdr-master/ -f
done

##------------------------------
sigfile=/FDR/list.txt
resultfile=/FDR/conjFDR/log
pleiofile=/pleiofdr-master
line_count=$(awk 'END{print NR}' "$sigfile")
for i in $(seq 1 $line_count);do
IDP=`cat $sigfile|awk 'NR=='$i' {print $1}'`
disease=`cat $sigfile|awk 'NR=='$i' {print $2}'`
name=$(echo "$disease" | awk -F'_' '{print $1}')

cmd="args=\"addpath(genpath('$pleiofile/'));runme_"$IDP"_"$name"_cjfdr\"\n/gpfs/chenglan/share/app/imaging/matlab2016b/bin/matlab -nodesktop -nosplash -r \"\$args\""
echo -e $cmd > "$resultfile/a"$IDP"_"$disease"_cfdr.sh"
qsub -q clc2 -e "$resultfile/a"$IDP"_"$disease"_cfdr.err" -o "$resultfile/a"$IDP"_"$disease"_cfdr.out" "$resultfile/a"$IDP"_"$disease"_cfdr.sh"
done


##------------------------------
sigfile=/FDR/list.txt
line_count=$(awk 'END{print NR}' "$sigfile")
for i in $(seq 1 $line_count);do
IDP=`cat $sigfile|awk 'NR=='$i' {print $1}'`
disease=`cat $sigfile|awk 'NR=='$i' {print $2}'`
name=$(echo "$disease" | awk -F'_' '{print $1}')

/anaconda3/envs/ldsc/bin/python2 /python_convert-master/fdrmat2csv.py \
--mat "/FDR/results//"$IDP"_"$name"/result.mat" \
--ref /pleiofdr-master/ref/9545380.ref \
--out "/FDR/results/"$IDP"_"$name"/result.mat.csv"

/anaconda3/envs/ldsc/bin/python2 /python_convert-master/sumstats.py clump \
--clump-field FDR \
--force  \
--plink plink \
--sumstats /FDR/results/"$IDP"_"$name"/result.mat.csv \
--bfile-chr /pleiofdr-master/chr@ \
--exclude-ranges '6:25119106-33854733' '8:7200000-12500000' \
--clump-p1 0.05 \
--out /FDR/results/"$IDP"_"$name"/AFUMA_leadsnp.csv
done


