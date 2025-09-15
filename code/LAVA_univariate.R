arg = commandArgs(T); info = arg[1]; outdir=arg[2];locifile=arg[3];
library(LAVA)
### Read in data
sum.id<-read.table(info,header=T)
loci = read.loci(locifile)
n.loc = nrow(loci)
input = process.input(input.info.file=info,           # input info file
                      sample.overlap.file=NULL,   # sample overlap file (can be set to NULL if there is no overlap)
                      ref.prefix="/gpfsnew/lab/groupYU/members/zhangyujie/SCZ/LAVA/g1000_eur/g1000_eur",                    # reference genotype data prefix
                      phenos=sum.id$phenotype)       # subset of phenotypes listed in the input info file that we want to process

t1=proc.time()
print(paste("Starting LAVA analysis for",n.loc,"loci"))
### Analyse
u=b=list()
for (i in 1:n.loc){
        print(i)
        locus = process.locus(loci[i,], input)             # process locus  
        # It is possible that the locus cannot be defined for various reasons (e.g. too few SNPs), so the !is.null(locus) check is necessary before calling the analysis functions.
        if (!is.null(locus)) {
                loc.info = data.frame(locus = locus$id, chr = locus$chr, start = locus$start, stop = locus$stop, n.snps = locus$n.snps, n.pcs = locus$K)
                # run the univariate tests
                loc.out = run.univ.bivar(locus)
                u[[i]] = cbind(loc.info, loc.out$univ)
                if(!is.null(loc.out$bivar)) b[[i]] = cbind(loc.info, loc.out$bivar)
        }
}
# save the output
write.table(do.call(rbind,u), paste0(outdir,sum.id[1,1],"_",sum.id[2,1],".univ.lava"), row.names=F,quote=F,col.names=T)
write.table(do.call(rbind,b), paste0(outdir,sum.id[1,1],"_",sum.id[2,1],".bivar.lava"), row.names=F,quote=F,col.names=T)

t2=proc.time()
print(paste("time:",t2[3]-t1[3]))
print(paste0("Done! Analysis output written to univ.lava and bivar.lava"))