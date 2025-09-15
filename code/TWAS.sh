sigfile=/TWAS/tissue_list

gwasfolder=/traitfold/
modelfolder=/predixcan/elastic_net_models
outfolder=/TWAS/result

for i in `cat $sigfile`; do
/MetaXcan/software/SPrediXcan.py \
--model_db_path $modelfolder/"en_Brain_$i.db" \
--covariance $modelfolder/"en_Brain_$i.txt.gz" \
--gwas_folder $gwasfolder \
--gwas_file_pattern SCZ.txt \
--snp_column rsid \
--effect_allele_column a1 \
--non_effect_allele_column a2 \
--beta_column beta \
--pvalue_column p \
--output_file "$outfolder/SCZ_"$i".csv"
done