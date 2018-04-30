#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args)<1) {
  stop("Missing config file.n", call.=FALSE)
} 

source(args[1])

match_size = function(mat, b_nbins, nbins, dsets) {
  if(b_nbins > nbins) {
    seq_bins = c()
    p = 1
    for (i in 1:dsets) {
      for (a in 1:b_nbins) {
        for (b in a:b_nbins) {
          if(a<=nbins && b<=nbins) {
            seq_bins = c(seq_bins,p)
          }
          p=p+1
        }	
      }
    }
    mat = mat[seq_bins]
  }
  if(b_nbins < nbins) {
    ext_mat = c()
    p = 1
    for (i in 1:dsets) {
      for (a in 1:nbins) {
        for (b in a:nbins) {
          if(a>b_nbins || b>b_nbins) {
            ext_mat = c (ext_mat,0.)
          } else {
            ext_mat = c (ext_mat,mat[p])
          }
          p=p+1
        }	
      }
    }
    mat = ext_mat
  }
  return(mat)
}

get_bin_data = function(csb,resolution)  {
  counts = data.table()
  bin_borders = csb@biases[, seq(min(pos) - 1, max(pos) + 
                                   1, resolution)]
  bins = cut(head(bin_borders, n = length(bin_borders) - 
                    1) + resolution/2, bin_borders, ordered_result = T, 
             right = F, include.lowest = T, dig.lab = 12)
  counts = CJ(name = csb@experiments[, name], bin1 = bins, 
              bin2 = bins)[bin2 >= bin1]
  poscounts = csb@counts[, .(name, bin1 = cut(pos1, bin_borders, 
                                              ordered_result = T, right = F, include.lowest = T, 
                                              dig.lab = 12), bin2 = cut(pos2, bin_borders, ordered_result = T, 
                                                                        right = F, include.lowest = T, dig.lab = 12), observed = contact.close + 
                               contact.far + contact.up + contact.down)][(!is.na(bin1)) & 
                                                                           (!is.na(bin2)), .(observed = sum(observed)), keyby = key(counts)]
  counts = poscounts[counts]
  counts[is.na(observed), `:=`(observed, 0)]
  counts[bin1 == min(bin1) & bin2 == max(bin2), `:=`(observed, 
                                                     observed + as.integer(1))]
  counts = add_bin_begin_and_end(counts)
  counts[, `:=`(distance, begin2 - begin1 + resolution/2)]
  return(counts)
}

get_autocorrelation = function(cssub, chunks,resolution) {
  chunks_norm=foreach (c=2:length(chunks),.errorhandling = 'remove',.combine=rbind) %do% {
    cs = cssub
    cs = zoom_csnorm(cs, chunks[c-1], chunks[c])
    matsub = get_bin_data(cs,resolution=resolution)
    c(chunks[c-1],chunks[c],cor(matsub$observed[-length(matsub$observed)],matsub$observed[-1]), sd(matsub$observed))
  }
}

if (exists("ncores")) ncores=as.integer(ncores) else ncores=4 #parallelize on so many processors
if (!exists("action")) action = 'normalize'

library(data.table)
library(binless)

if(action == 'normalize') {
  
  resolution = as.integer(resolution)
  
  if (exists("chr") && exists("beg") && exists("end")) {
    locus = c(paste(chr),as.integer(beg),as.integer(end))
  } else {
    locus = NULL
  }
  if (exists("fast_binless")) fast_binless = (fast_binless=='True') else fast_binless = F
  if (exists("bf_per_kb")) bf_per_kb=as.integer(bf_per_kb) else bf_per_kb=50 
  if (exists("bf_per_decade")) bf_per_decade=as.integer(bf_per_decade) else bf_per_decade=10
  if (exists("bins_per_bf")) bins_per_bf=as.integer(bins_per_bf) else bins_per_bf=10
  if (exists("ngibbs")) ngibbs=as.integer(ngibbs) else ngibbs=15 #maximum number of iterations
  if (exists("bg.steps")) bg.steps=as.integer(bg.steps) else bg.steps=5 #maximum number of steps where only the background model is fitted
  if (exists("iter")) iter=as.integer(iter) else iter=100
  if (exists("tol")) tol=as.numeric(ngibbs) else tol=1e-1 #relative tolerance on computed quantities upon convergence
  
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
  if(fast_binless) {
    if (exists("chr") && exists("beg") && exists("end")) {
      beg_ini=as.integer(beg)
      end_ini=as.integer(end)
    } else {
      beg_ini=csd@biases[,min(pos)]
      end_ini=csd@biases[,max(pos)]
    }
    chunks=seq.int(from=beg_ini,to=end_ini, by=2000000)
    chunks[length(chunks)] = end_ini
    npar = floor(ncores/4)
    npar_length = ceiling(length(chunks)/npar)
    chunks_par = list()
    cs_arr = list()
    csall = cs
    for(c in 1:npar) {
      cs = csall
      start_bin = ((c-1)*(npar_length-1)+1)
      end_bin = min(length(chunks),c*npar_length)
      chunks_par[[c]] = chunks[start_bin:end_bin]
      cs_arr[[c]] = zoom_csnorm(cs, chunks[start_bin], chunks[end_bin])
    }
    library(doParallel)
    cl<-makeCluster(ncores, type="FORK", outfile="normalize_binless.log")
    registerDoParallel(cl)
    #registerDoParallel(cores=ncores)
    stats=foreach (chunk=chunks_par,csa=cs_arr,.errorhandling = 'remove',.combine=rbind) %dopar% {
      get_autocorrelation(csa,chunk,resolution)
    }
    prmatrix(stats)
    all_chunks_norm=foreach (chunk=chunks_par,csa=cs_arr,.combine=rbind) %dopar% {
      chunks_norm=foreach (c=2:length(chunk),.errorhandling = 'remove',.combine=rbind) %do% {
        csb = csa
        csb = zoom_csnorm(csb, chunk[c-1], chunk[c])
        cat("*** Normalization of submatrix",chunk[c-1],"-",chunk[c],"\n")
        csb <- normalize_binless(csb, ngibbs = ngibbs, ncores = 4, base.res = base.res, bg.steps = bg.steps, tol = tol,
                                 bf_per_kb = bf_per_kb, bf_per_decade = bf_per_decade, bins_per_bf = bins_per_bf, iter=iter)  
        save(csb,file=paste0(dirname(dirname(rdata)),'/',chunk[c-1],"-",chunk[c],'.RData'))
        c(chunks[c-1],chunks[c],csb@par$lambda2,csb@par$alpha[1])
      }
    }
    #stopImplicitCluster()
    stopCluster(cl)
    prmatrix(all_chunks_norm)
    #idx_corr = which.min(abs(corr))
    #cs = csall
    #cs = zoom_csnorm(cs, chunks[idx_corr], chunks[idx_corr+1])
    #cs <- normalize_binless(cs, ngibbs = ngibbs, ncores = ncores, base.res = base.res, bg.steps = bg.steps, tol = tol,
    #		bf_per_kb = bf_per_kb, bf_per_decade = bf_per_decade, bins_per_bf = bins_per_bf, iter=iter)
    #lam2 = cs@par$lambda2
    #alpha = cs@par$alpha[1]
    #cs = csall
    lam2=mean(all_chunks_norm[,3])
    cat("*** Using average lambda2=",lam2,"\n")
    alpha=mean(all_chunks_norm[,4])
    cat("*** Using average alpha=",alpha,"\n")
    mat=binless:::bin_data(cs,resolution=resolution)
    if (exists("nouter")) nouter=as.integer(nouter) else nouter=25
    if (exists("tol_val")) tol_val=as.integer(tol_val) else tol_val=2e-1
    if (exists("bg_steps")) bg_steps=as.integer(bg_steps) else bg_steps=5
    if (exists("free_decay")) free_decay=as.integer(free_decay) else free_decay=10000
    out=binless:::fast_binless(mat, mat[,nlevels(bin1)], lam2, alpha, nouter, tol_val, bg_steps, free_decay)
    mat=as.data.table(out$mat)
    save(out,file=paste(rdata))
    mat_distance = exp(mat$distance)
  } else {
    cs <- normalize_binless(cs, ngibbs = ngibbs, ncores = ncores, base.res = base.res, bg.steps = bg.steps, tol = tol,
                            bf_per_kb = bf_per_kb, bf_per_decade = bf_per_decade, bins_per_bf = bins_per_bf, iter=iter)
    save(cs,file=paste(rdata))
    cs=bin_all_datasets(cs, ncores=ncores)
    mat=get_binned_matrices(cs)
    mat_distance = cs@par$decay$distance
  }
  mat_biasmat = mat$biasmat
  mat_decaymat = mat$decaymat
  #mat=mat[begin1<=as.integer(end) & begin2<=as.integer(end)]
  length_dset = length(mat_biasmat)/length(infiles)
  nbins = (as.integer(end)-as.integer(beg))/resolution
  b_nbins = as.integer(floor(sqrt(2*length_dset)))
  mat_biasmat=match_size(mat_biasmat, b_nbins, nbins, length(infiles)) 
  write.table(mat_decaymat[1:nbins],file = output_decay,row.names = FALSE,col.names = FALSE,sep = ',')
  write.table(mat_distance[1:nbins],file = output_distance,row.names = FALSE,col.names = FALSE,sep = ',')
  write.table(mat_biasmat,file = output_bias,row.names = FALSE,col.names = FALSE,sep = ',')
} else if(action == 'signal') {
  if (exists("nperf")) nperf=as.integer(nperf) else nperf=50
  if (exists("dataset_index")) dataset_index=as.integer(dataset_index) else dataset_index=0
  load(paste(rdata))
  if (exists("cs")) {
    cs=bin_all_datasets(cs, ncores=ncores)
    cs=detect_binless_interactions(cs, ncores = ncores, nperf = nperf)
    mat=get_binless_interactions(cs)
    mat=mat[begin1>=as.integer(start) & begin2>=as.integer(start) & begin1<=as.integer(end) & begin2<=as.integer(end)]
    mat_signal = mat$signal
    dsets = cs@experiments[,.N]
    length_dset = length(mat_signal)/dsets
  }
  if (exists("out")) {
    library('data.table')
    mat=as.data.table(out$mat)
    mat_signal = mat$signal
    dsets = 1
    length_dset = length(mat_signal)
  }
  nbins = (as.integer(end)-as.integer(start))/as.integer(resolution)
  b_nbins = as.integer(floor(sqrt(2*length_dset)))
  mat_signal=match_size(mat_signal, b_nbins, nbins, dsets) 
  if(dsets > 1) {
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
