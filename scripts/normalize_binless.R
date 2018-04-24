#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args)<1) {
	stop("Missing config file.n", call.=FALSE)
} 

source(args[1])

if (exists("fast_binless")) fast_binless = (fast_binless=='True') else fast_binless = F
if (exists("ncores")) ncores=as.integer(ncores) else ncores=4 #parallelize on so many processors
if (!exists("action")) action = 'normalize'

library(binless)

if(action == 'normalize') {

	if (exists("chr") && exists("beg") && exists("end")) {
		locus = c(paste(chr),as.integer(beg),as.integer(end))
	} else {
		locus = NULL
	}
	if (exists("bf_per_kb")) bf_per_kb=as.integer(bf_per_kb) else bf_per_kb=50 
	if (exists("bf_per_decade")) bf_per_decade=as.integer(bf_per_decade) else bf_per_decade=10
	if (exists("bins_per_bf")) bins_per_bf=as.integer(bins_per_bf) else bins_per_bf=10
	if (exists("ngibbs")) ngibbs=as.integer(ngibbs) else ngibbs=15 #maximum number of iterations
	if (exists("bg.steps")) bg.steps=as.integer(bg.steps) else bg.steps=5 #maximum number of steps where only the background model is fitted
	if (exists("iter")) iter=as.integer(iter) else iter=100
	if (exists("tol")) tol=as.numeric(ngibbs) else tol=1e-1 #relative tolerance on computed quantities upon convergence
	
	resolution = as.integer(resolution)
	base.res=resolution #base resolution for the fused lasso signal detection
	name="TADbit norm"
	
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
		dset=paste0("dset_",i)
		csd=read_and_prepare(paste0(infiles[i]), dset, "1", locus=locus,
				enzyme=paste(enzyme), name=paste(name,dset), circularize=-1, dangling.L=c(dleft),
				dangling.R=c(dright), read.len=read_lens[i], dmin=1000, replicate=i, skip=0)
		allcsd = c(allcsd,csd)
	}
	cs=merge_cs_norm_datasets(allcsd, different.decays="none")
	cs <- normalize_binless(cs, ngibbs = ngibbs, ncores = ncores, base.res = base.res, bg.steps = bg.steps, tol = tol,
							bf_per_kb = bf_per_kb, bf_per_decade = bf_per_decade, bins_per_bf = bins_per_bf, iter=iter)
	save(cs,file=paste(rdata))
	cs=bin_all_datasets(cs, ncores=ncores)
	mat=get_binned_matrices(cs)
	#mat=mat[begin1<=as.integer(end) & begin2<=as.integer(end)]
	mat_biasmat = mat$biasmat
	length_dset = length(mat_biasmat)/length(infiles)
	nbins = (as.integer(end)-as.integer(beg))/resolution
	b_nbins = as.integer(floor(sqrt(2*length_dset)))
	if(b_nbins > nbins) {
		seq_bins = c()
		p = 1
		for (i in 1:length(infiles)) {
			for (a in 1:b_nbins) {
				for (b in a:b_nbins) {
					if(a<=nbins && b<=nbins) {
						seq_bins = c(seq_bins,p)
					}
					p=p+1
				}	
			}
		}
		mat_biasmat = mat_biasmat[seq_bins]
	}
	write.table(mat$decaymat[1:nbins],file = output_decay,row.names = FALSE,col.names = FALSE,sep = ',')
	write.table(cs@par$decay$distance[1:nbins],file = output_distance,row.names = FALSE,col.names = FALSE,sep = ',')
	write.table(mat_biasmat,file = output_bias,row.names = FALSE,col.names = FALSE,sep = ',')
} else if(action == 'signal') {
	if (exists("nperf")) nperf=as.integer(nperf) else nperf=50
	if (exists("dataset_index")) dataset_index=as.integer(dataset_index) else dataset_index=0
	load(paste(rdata))
	cs=bin_all_datasets(cs, ncores=ncores)
	cs=detect_binless_interactions(cs, ncores = ncores, nperf = nperf)
	mat=get_binless_interactions(cs)
	mat=mat[begin1>=as.integer(start) & begin2>=as.integer(start) & begin1<=as.integer(end) & begin2<=as.integer(end)]
	mat_signal = mat$signal
	length_dset = length(mat_signal)/cs@experiments[,.N]
	nbins = (as.integer(end)-as.integer(start))/as.integer(resolution)
	b_nbins = as.integer(floor(sqrt(2*length_dset)))
	if(b_nbins > nbins) {
		seq_bins = c()
		p = 1
		for (i in 1:cs@experiments[,.N]) {
			for (a in 1:b_nbins) {
				for (b in a:b_nbins) {
					if(a<=nbins && b<=nbins) {
						seq_bins = c(seq_bins,p)
					}
					p=p+1
				}	
			}
		}
		mat_signal = mat_signal[seq_bins]
	}
	if(cs@experiments[,.N] > 1) {
		mat_signal = mat_signal[(1+length_dset*dataset_index):(length_dset*(dataset_index+1))]
	}
	write.table(mat_signal,file = output_signal,row.names = FALSE,col.names = FALSE,sep = ',')	
} else if(action == 'difference') {
	load(paste(rdata))
	if(cs@experiments[,.N] < 2) {
		cat("*** Dataset only contains one experiment. Use tadbit normalize with multiple jobids\n")
	} else {
		ref = cs@experiments[1,name]
		cs=bin_all_datasets(cs, ncores=ncores)
		cs=detect_binless_differences(cs, ref = ref, ncores=ncores)
		mat=get_binless_differences(cs, ref = ref)
		mat=mat[begin1>=as.integer(start) & begin2>=as.integer(start) & begin1<=as.integer(end) & begin2<=as.integer(end)]
		mat_diff = mat$difference
		length_dset = length(mat_diff)
		nbins = (as.integer(end)-as.integer(start))/as.integer(resolution)
		b_nbins = as.integer(floor(sqrt(2*length_dset)))
		if(b_nbins > nbins) {
		  seq_bins = c()
		  p = 1
		  for (i in 1:cs@experiments[,.N]) {
		    for (a in 1:b_nbins) {
		      for (b in a:b_nbins) {
		        if(a<=nbins && b<=nbins) {
		          seq_bins = c(seq_bins,p)
		        }
		        p=p+1
		      }	
		    }
		  }
		  mat_diff = mat_diff[seq_bins]
		}
		write.table(mat_diff,file = output_signal,row.names = FALSE,col.names = FALSE,sep = ',')
	}
}