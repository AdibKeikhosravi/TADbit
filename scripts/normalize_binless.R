#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args)<1) {
	stop("Missing config file.n", call.=FALSE)
} 

source(args[1])

if (exists("fast_binless")) { 
	fast_binless = (fast_binless=='True')
} else {
	fast_binless = F
}
if (exists("chr") && exists("beg") && exists("end")) {
	locus = c(paste(chr),as.integer(beg),as.integer(end))
} else {
	locus = NULL
}
resolution = as.integer(resolution)
name="TADbit norm"

library(binless)
library(foreach)
allcsd = c()
foreach (i=1:length(infiles)) %do% {
	a=examine_dataset(paste(infiles[i]),locus=locus,
			skip=0,nrows=100000, skip.fbm=T, read.len=read_lens[i])
	x = a$pdangling$data[cat=='left']
	y=setNames(aggregate(x$dist, by=list(x$dist), FUN=length),c("pos","val"))
	dleft = y[which.max(y$val),]$pos
	cat("*** detected dangling.L ",dleft,"\n")
	x = a$pdangling$data[cat=='right']
	y=setNames(aggregate(x$dist, by=list(x$dist), FUN=length),c("pos","val"))
	dright = y[which.max(y$val),]$pos
	cat("*** detected dangling.R ",dright,"\n")
	remove(a)
	dset=paste0("dset ",i)
	csd=read_and_prepare(paste0(infiles[i]), dset, "1", locus=locus,
			enzyme=paste(enzyme), name=paste(name,dset), circularize=-1, dangling.L=c(dleft),
			dangling.R=c(dright), read.len=read_lens[i], dmin=1000, replicate=i, skip=0)
	allcsd = c(allcsd,csd)
}
cs=merge_cs_norm_datasets(allcsd, different.decays="none")
ncores=4 #parallelize on so many processors
ngibbs=15 #maximum number of iterations
base.res=resolution #base resolution for the fused lasso signal detection
bg.steps=5 #maximum number of steps where only the background model is fitted
tol=1e-1 #relative tolerance on computed quantities upon convergence
cs <- normalize_binless(cs, ngibbs = ngibbs, ncores = ncores, base.res = base.res, bg.steps = bg.steps, tol = tol)
cs=bin_all_datasets(cs, ncores=ncores)
mat=get_binned_matrices(cs)
write.table(exp(mat$observed)/exp(mat$normalized),file = output_bias,row.names = FALSE,col.names = FALSE,sep = ',')
write.table(exp(cs@par$decay$log_decay),file = output_decay,row.names = FALSE,col.names = FALSE,sep = ',')
cs=detect_binless_interactions(cs, ncores=ncores)
mat=get_binless_interactions(cs)
write.table(mat$signal,file = output_signal,row.names = FALSE,col.names = FALSE,sep = ',')