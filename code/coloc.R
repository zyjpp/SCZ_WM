rm(list=ls())
library(coloc)
library(data.table)
library(openxlsx)
library(dplyr)
##gwas
GWAS_SCZ <- fread("/traitfold/SCZ.txt")
colnames(GWAS_SCZ)[colnames(GWAS_SCZ) %in% c("chr","rsid","pos","beta", "se")] <- c("Chr","SNP","BP","b", "SE")
##eqtl
eqtl <- fread("/coloc/eqtl/BrainMeta_cis_eqtl_gene_sig.txt")

maf <- fread("/g1000_eur/g1000_EUR_maf.frq")
eqtl <- left_join(eqtl,maf[,c(2,5)],by="SNP")


##loci
loci <- fread("/coloc/loci/loci_all.txt")
loci_unique <- loci[!duplicated(loci[,c("CHR","MinBP","MaxBP")]),]
loci_unique <- loci_unique[order(loci_unique$CHR,loci_unique$MinBP,loci_unique$MaxBP,decreasing = FALSE),]

##data
newall <- as.data.frame(matrix(NA,0,15))
new <- as.data.frame(matrix(NA,0,17))

for (i in 1:nrow(loci_unique)){
  
  ##eqtl
  ind1 <- which(eqtl$Chr==loci_unique$CHR[i] & eqtl$BP>=loci_unique$MinBP[i] & eqtl$BP<=loci_unique$MaxBP[i])
  eqtl.temp <- eqtl[ind1,]
  ##gwas
  ind2 <- which(GWAS_SCZ$Chr==loci_unique$CHR[i] & GWAS_SCZ$BP>=loci_unique$MinBP[i] & GWAS_SCZ$BP<=loci_unique$MaxBP[i])
  gwas.temp <- GWAS_SCZ[ind2,]
  gwas.temp <- gwas.temp[!duplicated(gwas.temp$SNP),]
  ##merge
  coloc_data <- merge(eqtl.temp, gwas.temp, by=c("Chr","BP","SNP"), all=FALSE, suffixes=c("_eqtl","_gwas"))
  
  if(nrow(coloc_data)!=0){
    gene <- unique(coloc_data$Probe)
    for(j in 1:length(gene)){
      cat("loci: ",i,"gene: ",j,"\n")
      gene1 <- gene[j]
      index <- which(coloc_data$Probe==gene1)
      coloc_data1 <- coloc_data[index,c("SNP","p_eqtl", "b_eqtl", "SE_eqtl", "p_gwas", "b_gwas", "SE_gwas", "MAF")]
      
      coloc_data1 <- na.omit(coloc_data1)
      coloc_data1 <- coloc_data1[!duplicated(coloc_data1$SNP),]
      if(nrow(coloc_data1)!=0){
      ##coloc
      result <- coloc.abf(dataset1=list(snp=coloc_data1$SNP,pvalues=coloc_data1$p_eqtl,beta=coloc_data1$b_eqtl,varbeta=(coloc_data1$SE_eqtl)^2,type="quant", N=2865), dataset2=list(snp=coloc_data1$SNP,pvalues=coloc_data1$p_gwas, beta=coloc_data1$b_gwas,varbeta=(coloc_data1$SE_gwas)^2, type="cc", N=130644),MAF=coloc_data1$MAF)
      PP.H4 <- result$summary[6]
      sumstat <- data.frame(c(loci_unique[i,],Probe=gene1,result$summary))
      
      colnames(newall) <- colnames(sumstat)
      newall <- rbind(newall,sumstat)
      if (PP.H4>0.8){
        o <- order(result$results$SNP.PP.H4,decreasing=TRUE)
        temp <- cbind(sumstat,result$results[o,][1,][,c("snp","SNP.PP.H4")])
        colnames(new) <- colnames(temp)
        new <- rbind(new,temp)
      }
      }
    }    
  }
}


