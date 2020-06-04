#########################################################################
# "A Bayesian model-free approach to combination therapy phase I trials
#  using censored time-to-toxicity data"
#
# R Code for simulating trials under the TITE-PIPE design
#
# Graham Wheeler (CRUK & UCL Cancer Trials Centre, University College London, UK)*
# Michael Sweeting (University of Leicester, UK)
# Adrian Mander (University of Cambridge, UK)
# 
# *Corresponding author: graham.wheeler@ucl.ac.uk
#########################################################################

#######################
# Packages to install #
#######################

install.packages("gtools")
install.packages("ggplot2")
install.packages("mgcv")
install.packages("VGAM")
install.packages("xtable")
library(gtools)
library(ggplot2)
library(mgcv)
library(VGAM)
library(xtable)


#####################
# Functions to load #
#####################

tite.pipe.design<-function (N, S = 1, c = 1, theta, pi = NULL, prior.med = NULL, 
                            prior.ss = NULL, strategy, admis, constraint = NULL, epsilon = NULL, 
                            a = NULL, b = NULL, alternate = FALSE, uppertox.constraint = NULL, 
                            stop = NULL, non.admissible = NULL, seed = NULL,
                            Tmax = 1, lambda = 1, first.full.eval = 0, wait.to.enter = FALSE, failure.type = "uniform", weight.type = "uniform",
                            min.cohort.size = 1, min.patient.type="complete.dose") 
{
  contour.select = "sim"
  reweight = FALSE
  R = 0
  P = 0
  if (wait.to.enter==FALSE & min.cohort.size!=first.full.eval) {
    stop("first.full.eval and min.cohort.size must be equal if wait.to.enter = FALSE")
  }
  if (is.null(prior.med) & is.null(a)) {
    stop("Either `prior.med' and `prior.ss' or `a' and `b' must be specified")
  }
  if (!is.null(prior.med)) {
    I = dim(prior.med)[1]
    J = dim(prior.med)[2]
  }
  else {
    I = dim(a)[1]
    J = dim(a)[2]
  }
  if (!(strategy %in% c("ss", "ss-random"))) 
    stop("strategy must be one of `ss' or `ss-random'")
  if (!(admis %in% c("adjacent", "closest"))) 
    stop("admis must be one of `adjacent' or `closest'")
  if (!is.null(constraint) & !(constraint %in% c("none", "neighbouring", 
                                                 "no.dose.skip", "neighbouring-nodiag", "no.dose.skip-nodiag"))) 
    stop("constraint must be one of `neighbouring', `no.dose.skip', `neighbouring-nodiag' or `no.dose.skip-nodiag'")
  if (is.null(pi)) {
    stop("pi (true DLT probabilities) must be specified")
  }
  if (!is.null(uppertox.constraint)) {
    if ((theta > uppertox.constraint | uppertox.constraint > 
         1)) 
      stop("uppertox.constraint must be a number between theta and 1")
  }
  if (!is.null(epsilon)) {
    if (epsilon < 0 | epsilon > 1) 
      stop("epsilon must be a number between 0 and 1")
  }
  if (ceiling(N/c) != floor(N/c)) 
    stop("Total sample size `N' must be divisible by cohort size `c'")
  if (!is.null(non.admissible)) {
    if (dim(non.admissible)[1] != I | dim(non.admissible)[2] != 
        J) 
      stop(paste("Expecting non.admissible to be an ", 
                 I, "x", J, " matrix", sep = ""))
  }
  if (!is.null(non.admissible)) {
    if (!is.logical(non.admissible)) 
      stop("non.admissible must be a matrix of logicals")
  }
  k = I * J
  matrices <- monotonic.matrices(I, J)
  hsegments <- lapply(matrices, function(m) {
    mp <- rbind(rep(0, J), m, rep(1, J))
    mp[-1, ] - mp[-(I + 2), ]
  })
  vsegments <- lapply(matrices, function(m) {
    mp <- cbind(rep(0, I), m, rep(1, I))
    mp[, -1] - mp[, -(J + 2)]
  })
  doses <- 1
  r.sim <- n.sim <- list()
  arrival.list<-weights.list<-original.arrival.list<-trial.duration.list<-failure.list<-consider.stop<-list() ###
  rec.i.sim <- rec.j.sim <- array(NA, dim = c(S, N/c, doses))
  rec <- matrix(0, ncol = J, nrow = I)
  n.rpII <- vector()
  if (is.null(a) & is.null(b)) {
    prior <- beta.med(prior.med, prior.ss)
    a <- prior$a
    b <- prior$b
  }
  mat.list = uppermat.list = uppermat2.list = n.list = r.list = pi.theta.list = rpII.list = list()
  means <- cdfs <- list()
  mat.lik.list <- list()
  h.lik.list <- v.lik.list <- list()
  dom.list <- admis.list <- list()
  for (s in 1:S) {
    if(seed==TRUE){set.seed(s)} 
    consider.stop.vec<-NULL
    arrivals<-interarrival.time.fn(n = N-1, rate = lambda, Tmax = Tmax)
    arrival.times<-arrivals$S.startend
    failure<-list()
    r <- r.stop <- matrix(0, nrow = I, ncol = J)
    n <- n.stop <- matrix(0, nrow = I, ncol = J)
    p <- pbeta(theta, a, b)
    if (!is.null(uppertox.constraint)) {
      pconstraint <- pbeta(uppertox.constraint, a, b)
    }
    else {
      pconstraint <- NULL
    }
    rec.i = rec.j = matrix(nrow = 0, ncol = doses)
    create <- mtc.create(matrices, p, constraint, pconstraint, 
                         epsilon, admis, rec.i, rec.j, n, contour.select = contour.select, 
                         hsegments = hsegments, vsegments = vsegments, S, 
                         non.admissible, reweight, R, P)
    mat <- create$mat
    pi.theta <- 1/(a + b)
    nxt <- mtc(create$dominant, create$admissible, strategy, 
               rec.i, rec.j, pi.theta, mat, p, alternate, psmooth = create$mat.lik)
    if (S == 1) {
      thetas <- seq(0.05, 0.95, 0.05)
      ps <- lapply(thetas, pbeta, shape1 = a, shape2 = b)
      mtc.nums <- lapply(ps, function(p) unlist(lapply(matrices, 
                                                       function(l) {
                                                         prod((1 - p)[l == 1]) * prod(p[l == 0])
                                                       })))
      mtc.liks <- lapply(mtc.nums, function(mtc.num) mtc.num/sum(mtc.num))
      mat.liks <- lapply(mtc.liks, function(mtc.lik) sapply(1:length(mtc.lik), 
                                                            function(l) {
                                                              matrices[[l]] * mtc.lik[[l]]
                                                            }, simplify = F))
      cdfs[[1]] <- lapply(mat.liks, function(mat.lik) Reduce("+", mat.lik))
      means[[1]] <- 0.95 - 0.05 * Reduce("+", cdfs[[1]])
      mat.list[[1]] = create$mat.mode
      mat.lik.list[[1]] <- create$mat.lik
      h.lik.list[[1]] <- create$h.lik
      v.lik.list[[1]] <- create$v.lik
      uppermat.list[[1]] = create$matupper
      uppermat2.list[[1]] = create$matupper2
      dom.list[[1]] <- create$dominant
      admis.list[[1]] <- create$admissible
      n.list[[1]] = n
      r.list[[1]] = r
      pi.theta.list[[1]] = pi.theta
    }
    rec.i = nxt$rec.i
    rec.j = nxt$rec.j
    for(m in 1:N){ 
      if (doses == 1) {
        if (!is.null(pi)) {
          failure[[m]]<-failure.time.fn(pi[rec.i[m,1], rec.j[m,1]], type=failure.type, Tmax=Tmax, t.i0=arrival.times[m])
          if(wait.to.enter==TRUE){
            if(n[rec.i[m,1], rec.j[m,1]]+1 == min.cohort.size){
              fails<-lapply(seq(m,m - min.cohort.size + 1, by=-1), function(x) failure[[x]])
              fail.inds<-fail.ind.fn(fails, min.cohort.size)
              arrival.inds<-arrivals.fn(arrival.times[m:(m - min.cohort.size + 1)], Tmax)
              fail.index<-which(fail.inds!=0)
              comb.ind.vec<-replace(arrival.inds, fail.index, fail.inds[fail.index])
              index<-which(arrival.times[(m+1):(N+1)]<max(comb.ind.vec))+m # edited colon
              if(arrival.times[m+1] < max(comb.ind.vec)){arrival.times[index]<-max(comb.ind.vec)}
            }
            if(n[rec.i[m,1], rec.j[m,1]]+1 > min.cohort.size){
              fail.inds<-failure[[m]]$time.i*(failure[[m]]$y.i==1)
              arrival.inds<-arrival.times[m]+Tmax
              fail.index<-which(fail.inds!=0)
              comb.ind.vec<-replace(arrival.inds, fail.index, fail.inds[fail.index])
              index<-which(arrival.times[(m+1):(N+1)]<max(comb.ind.vec))+m # edited colon
              if(arrival.times[m+1] < max(comb.ind.vec)){arrival.times[index]<-max(comb.ind.vec)}
            }
          }else{
            if(m<=first.full.eval & n[rec.i[m,1], rec.j[m,1]]+1==min.cohort.size){
              fails<-lapply(seq(m,m - min.cohort.size + 1, by=-1), function(x) failure[[x]])
              fail.inds<-fail.ind.fn(fails, min.cohort.size)
              arrival.inds<-arrivals.fn(arrival.times[m:(m - min.cohort.size + 1)], Tmax)
              fail.index<-which(fail.inds!=0)
              comb.ind.vec<-replace(arrival.inds, fail.index, fail.inds[fail.index])
              index<-which(arrival.times[(m+1):(N+1)]<max(comb.ind.vec))+m
              if(arrival.times[m+1] < max(comb.ind.vec)){arrival.times[index]<-max(comb.ind.vec)}
            }
            if(n[rec.i[m,1], rec.j[m,1]]+1==min.cohort.size & min.patient.type=="complete.dose"){
              fails<-lapply(seq(m,m - min.cohort.size + 1, by=-1), function(x) failure[[x]])
              fail.inds<-fail.ind.fn(fails, min.cohort.size)#
              arrival.inds<-arrivals.fn(arrival.times[m:(m - min.cohort.size + 1)], Tmax)#
              fail.index<-which(fail.inds!=0)#
              comb.ind.vec<-replace(arrival.inds, fail.index, fail.inds[fail.index])
              index<-which(arrival.times[(m+1):(N+1)]<max(comb.ind.vec))+m
              if(arrival.times[m+1] < max(comb.ind.vec)){arrival.times[index]<-max(comb.ind.vec)}
            }
          }
          if(wait.to.enter==TRUE & m==N & failure[[m]]$y.i==0){arrival.times[m+1]<-arrival.times[m]+Tmax}  
          if(m==N){
            end.times<-arrival.times+Tmax
            dltend.times<-sapply(1:N, function(z) ifelse(failure[[z]]$y.i==1, failure[[z]]$time.i, end.times[z]))
            arrival.times[N+1]<-max(dltend.times)}
          weight.mat<-weight.mat.fn(times=arrival.times[1:(m+1)], dlt.tox=sapply(1:m, function(z) failure[[z]]$y.i), dlt.times=sapply(1:m, function(z) failure[[z]]$time.i), weight.type = weight.type, Tmax=Tmax) ###-arrival.times[z]
          dose.mat<-cbind(rec.i, rec.j)
          unique.dose.mat<-uniquecombs(dose.mat)
          if(is.null(dim(unique.dose.mat)[1])==TRUE){unique.dose.mat<-matrix(unique.dose.mat, ncol=2, byrow=T)}
          for(i in 1:dim(unique.dose.mat)[1]){
            index<-which(sapply(1:dim(dose.mat)[1], function(z) identical(dose.mat[z,], unique.dose.mat[i,]))==TRUE)
            r[unique.dose.mat[i,1], unique.dose.mat[i,2]]<-ifelse(length(index)==1, weight.mat[index,m+1], apply(weight.mat[index, ], 2, sum)[m+1])
          }
          n[rec.i[m, 1], rec.j[m, 1]] <- n[rec.i[m, 1], rec.j[m, 1]] + c
        }
        else {
          break
        }
      }
      else {
        r[rec.i[m, 1], rec.j[m, 1]] <- r[rec.i[m, 1], rec.j[m, 1]] + rbinom(1, c/2, pi[rec.i[m, 1], rec.j[m, 1]])
        n[rec.i[m, 1], rec.j[m, 1]] <- n[rec.i[m, 1], rec.j[m, 1]] + c/2
        r[rec.i[m, 2], rec.j[m, 2]] <- r[rec.i[m, 2], rec.j[m, 2]] +  rbinom(1, c/2, pi[rec.i[m, 2], rec.j[m, 2]])
        n[rec.i[m, 2], rec.j[m, 2]] <- n[rec.i[m, 2], rec.j[m, 2]] + c/2
      }
      p <- pbeta(theta, a + r, b + n - r)
      if (!is.null(uppertox.constraint)) {
        pconstraint <- pbeta(uppertox.constraint, a + r, b + n - r)
      }
      else {
        pconstraint <- NULL
      }
      create <- mtc.create(matrices, p, constraint, pconstraint, 
                           epsilon, admis, rec.i, rec.j, n, contour.select = contour.select, 
                           hsegments = hsegments, vsegments = vsegments, 
                           S, non.admissible, reweight, R, P)
      mat <- create$mat
      pi.theta <- 1/(a + b + n)
      if(wait.to.enter==FALSE){
        stop.dat.ind<-which(arrival.times[1:m]!=arrival.times[m+1] & (weight.mat[1:m,m+1]==0 | weight.mat[1:m,m+1]==1))
        rec.i.stop<-matrix(rec.i[stop.dat.ind,],ncol=1)
        rec.j.stop<-matrix(rec.j[stop.dat.ind,], ncol=1)
        stop.dose.mat<-matrix(dose.mat[stop.dat.ind,], ncol=2, byrow=F)
        if(length(stop.dat.ind)!=0){
          if(length(stop.dat.ind)==1){stop.weight.mat<-weight.mat[stop.dat.ind,m+1]
          }else{stop.weight.mat<-weight.mat[1:m,][stop.dat.ind,]}
        }
        if(dim(stop.dose.mat)[1]!=0){
          unique.stop.dose.mat<-uniquecombs(stop.dose.mat)
          if(is.null(dim(unique.stop.dose.mat)[1])==TRUE){unique.stop.dose.mat<-matrix(unique.stop.dose.mat, ncol=2, byrow=T)}
          for(i in 1:dim(unique.stop.dose.mat)[1]){
            index<-which(sapply(1:dim(stop.dose.mat)[1], function(z) identical(stop.dose.mat[z,], unique.stop.dose.mat[i,]))==TRUE) ###
            r.stop[unique.stop.dose.mat[i,1], unique.stop.dose.mat[i,2]]<-ifelse(length(index)==1, ifelse(length(stop.weight.mat)==1, stop.weight.mat, stop.weight.mat[index,m+1]), apply(stop.weight.mat[index, ], 2, sum)[m+1])
            n.stop[unique.stop.dose.mat[i,1], unique.stop.dose.mat[i,2]]<-length(index)
          }
          p.stop<-pbeta(theta, a + r.stop, b + n.stop - r.stop)
          create.stop <- mtc.create(matrices, p.stop, constraint, pconstraint, 
                                    epsilon, admis, rec.i.stop, rec.j.stop, n.stop, contour.select = contour.select, 
                                    hsegments = hsegments, vsegments = vsegments, 
                                    S, non.admissible, reweight, R, P)
          if (all(!create.stop$admissible)) {
            rec.i <- rbind(rec.i, matrix(0, nrow = N/c + 1 - m, ncol = doses))
            rec.j <- rbind(rec.j, matrix(0, nrow = N/c + 1 - m, ncol = doses))
            break
          }
        }
      }else{      if (all(!create$admissible)) {
        rec.i <- rbind(rec.i, matrix(0, nrow = N/c + 1 - m, ncol = doses))
        rec.j <- rbind(rec.j, matrix(0, nrow = N/c + 1 - m, ncol = doses))
        break
      }
      }
      if (!is.null(stop)) {
        stop.dose.mat<-dose.mat[which(weight.mat[1:m,m+1]==0 | weight.mat[1:m,m+1]==1),]
        if(dim(stop.dose.mat)[1]!=0){
          stop.weight.mat<-weight.mat[1:m,][which(weight.mat[1:m,m+1]==0 | weight.mat[1:m,m+1]==1),]
          unique.stop.dose.mat<-uniquecombs(stop.dose.mat)
          if(is.null(dim(unique.stop.dose.mat)[1])==TRUE){unique.stop.dose.mat<-matrix(unique.stop.dose.mat, ncol=2, byrow=T)}
          for(i in 1:dim(unique.stop.dose.mat)[1]){
            index<-which(sapply(1:dim(stop.dose.mat)[1], function(z) identical(stop.dose.mat[z,], unique.stop.dose.mat[i,]))==TRUE) ###
            r.stop[unique.stop.dose.mat[i,1], unique.stop.dose.mat[i,2]]<-ifelse(length(index)==1, stop.weight.mat[index,m+1], apply(stop.weight.mat[index, ], 2, sum)[m+1])
            n.stop[unique.stop.dose.mat[i,1], unique.stop.dose.mat[i,2]]<-length(index)
          }
          p.stop<-pbeta(theta, a + r.stop, b + n.stop - r.stop)
          if (1 - p.stop[1, 1] > stop) {
            rec.i <- rbind(rec.i, matrix(0, nrow = N/c + 1 - m, ncol = doses))
            rec.j <- rbind(rec.j, matrix(0, nrow = N/c + 1 - m, ncol = doses))
            break
          }
        }
      }
      if(n[rec.i[m, 1], rec.j[m, 1]] < min.cohort.size){
        rec.i = rbind(rec.i, rec.i[m,1])
        rec.j = rbind(rec.j, rec.j[m,1])
      } else{
        nxt <- mtc(create$dominant, create$admissible, strategy, 
                   rec.i, rec.j, pi.theta, mat, p, alternate, psmooth = create$mat.lik)
        if(is.na(nxt$rec.i[m+1,1]) | is.na(nxt$rec.j[m+1,1])){
          rec.i<-matrix(rec.i[1:m,],ncol=1)
          rec.j<-matrix(rec.j[1:m,],ncol=1)
          consider.stop.vec<-cbind(consider.stop.vec, c(m,arrival.times[m]))
          pat.index<-which(weight.mat[1:m,m+1]!=0 & (weight.mat[1:m,m+1]!=1 | (weight.mat[1:m,m+1]==1 & arrival.times[1:m]==arrival.times[m+1])))
          fails<-lapply(pat.index, function(x) failure[[x]])
          fail.inds<-fail.ind.fn(fails, length(pat.index))
          arrival.inds<-arrivals.fn(arrival.times[pat.index], Tmax)
          fail.index<-which(fail.inds!=0)
          comb.ind.vec<-replace(arrival.inds, fail.index, fail.inds[fail.index])
          index<-which(arrival.times[(m+1):(N+1)]<max(comb.ind.vec))+m # edited colon
          if(arrival.times[m+1] < max(comb.ind.vec)){arrival.times[index]<-max(comb.ind.vec)}
          if(wait.to.enter==TRUE & m==N & failure[[m]]$y.i==0){arrival.times[m+1]<-arrival.times[m]+Tmax}  
          weight.mat<-weight.mat.fn(times=arrival.times[1:(m+1)], dlt.tox=sapply(1:m, function(z) failure[[z]]$y.i), dlt.times=sapply(1:m, function(z) failure[[z]]$time.i), weight.type = weight.type, Tmax=Tmax) ###-arrival.times[z]
          dose.mat<-cbind(rec.i, rec.j)
          unique.dose.mat<-uniquecombs(dose.mat)
          if(is.null(dim(unique.dose.mat)[1])==TRUE){unique.dose.mat<-matrix(unique.dose.mat, ncol=2, byrow=T)}
          for(i in 1:dim(unique.dose.mat)[1]){
            index<-which(sapply(1:dim(dose.mat)[1], function(z) identical(dose.mat[z,], unique.dose.mat[i,]))==TRUE) ###
            r[unique.dose.mat[i,1], unique.dose.mat[i,2]]<-ifelse(length(index)==1, weight.mat[index,m+1], apply(weight.mat[index, ], 2, sum)[m+1])
          }
          p <- pbeta(theta, a + r, b + n - r)
          if (!is.null(uppertox.constraint)) {
            pconstraint <- pbeta(uppertox.constraint, a + r, b + n - r)
          }
          else {
            pconstraint <- NULL
          }
          create <- mtc.create(matrices, p, constraint, pconstraint, 
                               epsilon, admis, rec.i, rec.j, n, contour.select = contour.select, 
                               hsegments = hsegments, vsegments = vsegments, 
                               S, non.admissible, reweight, R, P)
          mat <- create$mat
          pi.theta <- 1/(a + b + n)
          if (all(!create$admissible)) {
            rec.i <- rbind(rec.i, matrix(0, nrow = N/c + 1 - m, ncol = doses))
            rec.j <- rbind(rec.j, matrix(0, nrow = N/c + 1 - m, ncol = doses))
            break
          }
          nxt<-mtc(create$dominant, create$admissible, strategy, 
                   rec.i, rec.j, pi.theta, mat, p, alternate, psmooth = create$mat.lik)
        }
        rec.i = nxt$rec.i
        rec.j = nxt$rec.j
      }
    }
    consider.stop[[s]]<-consider.stop.vec
    r.sim[[s]] <- r
    n.sim[[s]] <- n
    failure.list[[s]]<-failure
    rec.i.sim[s, , ] <- rec.i[-(m + 1), ]
    rec.j.sim[s, , ] <- rec.j[-(m + 1), ]
    if (any(create$admissible)) {
      create.rpII <- mtc.create(matrices, p, constraint = "none", 
                                pconstraint = pconstraint, epsilon, admis = "closest", 
                                rec.i, rec.j, n, contour.select = contour.select, 
                                hsegments = hsegments, vsegments = vsegments, 
                                S, non.admissible, reweight = reweight, R = R, 
                                P = P)
      rpIIs <- create.rpII$dominant & create.rpII$mat == 0 & n != 0 & create$matupper == 0 & create$matupper2 == 0
      rpII.i <- row(mat)[rpIIs]
      rpII.j <- col(mat)[rpIIs]
      for (i in 1:length(rpII.i)) {
        rec[rpII.i[i], rpII.j[i]] <- rec[rpII.i[i], rpII.j[i]] + 1
      }
      n.rpII[s] <- length(rpII.i)
    }
    else {
      rpII.i <- rpII.j <- numeric(0)
      n.rpII[s] <- 0
    }
    rpII <- rbind(rpII.i, rpII.j)
    rownames(rpII) <- c("rpII.A", "rpII.B")
    rpII.list[[s]] <- rpII
    arrival.list[[s]]<-arrival.times 
    original.arrival.list[[s]]<-arrivals$S
    weights.list[[s]]<-weight.mat
    trial.duration.list[[s]]<-arrival.times[sum(n)+1]
    if (S > 1 & s%%50==0) 
      cat(s, "\n")
  }
  exp <- Reduce("+", n.sim)/sum(Reduce("+", n.sim))
  no.not.treated <- N * S - sum(Reduce("+", n.sim))
  dlts <- sapply(1:S, function(l) {
    sum(r.sim[[l]])/sum(n.sim[[l]])
  })
  results <- list(r.sim = r.sim, n.sim = n.sim, rec.i.sim = rec.i.sim, 
                  rec.j.sim = rec.j.sim, exp = exp, rec = rec, p.rec = rec/sum(rec), 
                  dlts = dlts, mat.list = mat.list, uppermat.list = uppermat.list, 
                  uppermat2.list = uppermat2.list, r.list = r.list, n.list = n.list, 
                  n.rpII = n.rpII, no.not.treated = no.not.treated, pi = pi, 
                  theta = theta, a = a, b = b, pi.theta.list = pi.theta.list, 
                  cdfs = cdfs, means = means, mat.lik.list = mat.lik.list, 
                  h.lik.list = h.lik.list, v.lik.list = v.lik.list, dom.list = dom.list, 
                  admis.list = admis.list, rpII.list = rpII.list,
                  arrival.list = arrival.list, original.arrival.list = original.arrival.list, weights.list = weights.list, trial.duration.list = trial.duration.list, failure.list = failure.list, wait.to.enter = wait.to.enter, consider.stop = consider.stop, min.patient.type = min.patient.type) ###
  if (S > 1) {
    class(results) <- "pipe.sim"
  }
  else {
    class(results) <- "pipe"
  }
  return(results)
}



fail.ind.fn<-function(failure, n){
  sapply(1:n, function(z) failure[[z]]$time.i * (failure[[z]]$y.i==1))
}



arrivals.fn<-function(arrivals, Tmax){
  sapply(1:length(arrivals), function(z) arrivals[z]+Tmax)
}



interarrival.time.fn<-function(n = n, rate = rate, Tmax = Tmax){
  X <- rexp(n = n, rate = rate) # Inter-arrival times   
  S <-c(0, cumsum(X))
  S.startend<-c(S,S[length(S)]+Tmax)
  return(list(X = X, S = S, S.startend = S.startend))
  # X = interarrival times
  # S = original arrival times
  # S.start = same as S with end time included also
}



monotonic.matrices<-function (I, J) 
{
  comb.col <- combinations(2, J, c(0, 1), repeats.allowed = TRUE)
  n.com1 <- dim(comb.col)[1]
  comb.row <- combinations(n.com1, I, repeats.allowed = TRUE)
  n.com2 <- dim(comb.row)[1]
  matrices <- sapply(1:n.com2, function(i) {
    comb.col[comb.row[i, ], ]
  }, simplify = F)
  return(matrices)
}



beta.med<-function (prior.med, prior.ss) 
{
  if (any(prior.med == 0 | prior.med == 1)) 
    stop("`prior.med' must be greater than 0 and less than 1")
  betaprior1 = function(K, x, p) {
    m.lo = 0
    m.hi = 1
    flag = 0
    while (flag == 0) {
      m0 = (m.lo + m.hi)/2
      p0 = pbeta(x, K * m0, K * (1 - m0))
      if (p0 < p) 
        m.hi = m0
      else m.lo = m0
      if (abs(p0 - p) < 1e-04) 
        flag = 1
    }
    return(m0)
  }
  a.med <- b.med <- matrix(NA, nrow = dim(prior.med)[1], ncol = dim(prior.med)[2])
  for (i in 1:dim(prior.med)[1]) {
    for (j in 1:dim(prior.med)[2]) {
      a.med[i, j] <- prior.ss[i, j] * betaprior1(prior.ss[i, 
                                                          j], prior.med[i, j], 0.5)
      b.med[i, j] <- prior.ss[i, j] - a.med[i, j]
    }
  }
  return(list(a = a.med, b = b.med))
}



mtc.create<-function (matrices, p, constraint, pconstraint, epsilon, admis, 
                      rec.i, rec.j, n, contour.select = contour.select, hsegments = hsegments, 
                      vsegments = vsegments, S, non.admissible, reweight, R, P) 
{
  m <- dim(rec.i)[1]
  mtc.num <- unlist(lapply(matrices, function(l) {
    prod((1 - p)[l == 1]) * prod(p[l == 0])
  }))
  mtc.lik <- mtc.num/sum(mtc.num)
  if (reweight) {
    I <- dim(matrices[[1]])[1]
    J <- dim(matrices[[1]])[2]
    Q <- unlist(lapply(matrices, sum))
    newweight <- 1 + R * abs(2/(I * J) * Q - 1)^P
    mtc.lik <- newweight * mtc.lik/sum(newweight * mtc.lik)
  }
  if (contour.select == "median" | S == 1) {
    h.lik <- Reduce("+", lapply(1:length(mtc.lik), function(i) mtc.lik[[i]] * 
                                  hsegments[[i]]))
    v.lik <- Reduce("+", lapply(1:length(mtc.lik), function(i) mtc.lik[[i]] * 
                                  vsegments[[i]]))
    nv <- ncol(h.lik)
    nh <- nrow(v.lik)
    h.med <- apply(h.lik, 2, w.median, x = 0:nh)
    v.med <- apply(v.lik, 1, w.median, x = 0:nv)
    mat.hmed <- do.call(cbind, lapply(h.med, function(i) rep(0:1, 
                                                             c(i, nh - i))))
    mat.vmed <- do.call(rbind, lapply(v.med, function(i) rep(0:1, 
                                                             c(i, nv - i))))
    mtc.hmed <- which(sapply(matrices, function(m) sum(abs(m - 
                                                             mat.hmed)) == 0))
    mtc.vmed <- which(sapply(matrices, function(m) sum(abs(m - 
                                                             mat.vmed)) == 0))
    if (mtc.lik[mtc.hmed] >= mtc.lik[mtc.vmed]) 
      mat.med <- mat.hmed
    else mat.med <- mat.vmed
    mtc.mode <- which.max(mtc.lik)
    mat.mode <- matrices[[mtc.mode]]
  }
  else {
    mtc.mode <- which.max(mtc.lik)
    mat.mode <- matrices[[mtc.mode]]
    h.lik <- v.lik <- mat.med <- NULL
  }
  if (contour.select == "median") 
    mat <- mat.med
  else mat <- mat.mode
  I <- dim(mat)[1]
  J <- dim(mat)[2]
  if (!is.null(pconstraint)) {
    upper.lik <- which.max(unlist(lapply(matrices, function(l) {
      prod((1 - pconstraint)[l == 1]) * prod(pconstraint[l == 
                                                           0])
    })))
    matupper <- matrices[[upper.lik]]
  }
  else {
    matupper <- matrix(0, nrow = I, ncol = J)
  }
  mat.lik <- Reduce("+", sapply(1:length(mtc.lik), function(l) {
    matrices[[l]] * mtc.lik[[l]]
  }, simplify = F))
  if (!is.null(epsilon)) {
    weight.pMTC <- mat.lik
    matupper2 <- weight.pMTC >= epsilon
  }
  else {
    weight.pMTC <- NULL
    matupper2 <- matrix(0, nrow = I, ncol = J)
  }
  if (grepl("neighbouring", constraint)) {
    if (dim(rec.i)[2] > 1) {
      admissible1 <- row(mat) <= max(rec.i[m, 1], 0) + 
        1 & col(mat) <= max(rec.j[m, 1], 0) + 1 & row(mat) >= 
        max(rec.i[m, 1], 0) - 1 & col(mat) >= max(rec.j[m, 
                                                        1], 0) - 1
      admissible2 <- row(mat) <= max(rec.i[m, 2], 0) + 
        1 & col(mat) <= max(rec.j[m, 2], 0) + 1 & row(mat) >= 
        max(rec.i[m, 2], 0) - 1 & col(mat) >= max(rec.j[m, 
                                                        2], 0) - 1
      admissible <- admissible1 | admissible2
    }
    else {
      admissible <- row(mat) <= max(rec.i[m, 1], 0) + 1 & 
        col(mat) <= max(rec.j[m, 1], 0) + 1 & row(mat) >= 
        max(rec.i[m, 1], 0) - 1 & col(mat) >= max(rec.j[m, 
                                                        1], 0) - 1
    }
  }
  else if (grepl("no.dose.skip", constraint)) {
    if (dim(rec.i)[2] > 1) {
      admissible1 <- row(mat) <= max(rec.i[, 1], 0) + 1 & 
        col(mat) <= max(rec.j[, 1], 0) + 1
      admissible2 <- row(mat) <= max(rec.i[, 2], 0) + 1 & 
        col(mat) <= max(rec.j[, 2], 0) + 1
      admissible <- admissible1 | admissible2
    }
    else {
      admissible <- row(mat) <= max(rec.i[, 1], 0) + 1 & 
        col(mat) <= max(rec.j[, 1], 0) + 1
    }
  }
  else {
    admissible <- matrix(TRUE, nrow = I, ncol = J)
  }
  if (grepl("nodiag", constraint)) {
    if (dim(rec.i)[1] >= 1) {
      if (rec.i[m, 1] != I & rec.j[m, 1] != J) {
        admissible[rec.i[m, 1] + 1, rec.j[m, 1] + 1] <- FALSE
      }
    }
  }
  if (!is.null(non.admissible)) {
    admissible <- admissible & !non.admissible
  }
  if (!is.null(pconstraint) | !is.null(epsilon)) {
    admissible <- admissible & matupper == 0 & matupper2 == 
      0
    if (all(!admissible)) {
      test <- abs(rec.i[m, 1] - row(mat)) + abs(rec.j[m, 
                                                      1] - col(mat))
      admissible <- test == min(c(test[matupper == 0 & 
                                         matupper2 == 0], -Inf)) & matupper == 0 & matupper2 == 
        0
    }
  }
  separate <- FALSE
  if (admis == "adjacent") {
    if (I < 2 | J < 2) 
      stop("Admissible doses can only be calculated when both drugs have more than one level")
    admat <- mat
    dominantu <- admat == 1 & (rbind(0, admat[-I, ]) == 0 | 
                                 cbind(0, admat[, -J]) == 0 | rbind(0, cbind(0, admat[, 
                                                                                      -J])[-I, ]) == 0)
    dominantl <- admat == 0 & (rbind(admat[-1, ], 1) == 1 | 
                                 cbind(admat[, -1], 1) == 1 | rbind(cbind(admat[, 
                                                                                -1], 1)[-1, ], 1) == 1)
    dominant <- dominantl | dominantu
    if (!any(dominant & admissible)) {
      separate <- TRUE
    }
  }
  if (admis == "closest" | separate == TRUE) {
    admat <- mat
    admat[admissible == FALSE] = 2
    dominant <- closest(admat)
  }
  return(list(dominant = dominant, admissible = admissible, 
              mat = mat, mat.mode = mat.mode, mat.med = mat.med, matupper = matupper, 
              matupper2 = matupper2, weight.pMTC = weight.pMTC, mat.lik = mat.lik, 
              h.lik = h.lik, v.lik = v.lik))
}



w.median<-function (x, w) 
{
  if (missing(w)) 
    w <- rep(1, length(x))
  ok <- complete.cases(x, w)
  x <- x[ok]
  w <- w[ok]
  ind <- sort.list(x)
  x <- x[ind]
  w <- w[ind]
  ind1 <- min(which(cumsum(w)/sum(w) >= 0.5))
  ind2 <- if ((w[1]/sum(w)) > 0.5) {
    1
  }
  else {
    max(which(cumsum(w)/sum(w) <= 0.5))
  }
  max(x[ind1], x[ind2])
}



closest<-function (mat) 
{
  I <- nrow(mat)
  J <- ncol(mat)
  dominantu <- mat == 1 & rbind(0, mat[-I, ]) %in% c(0, 2) & 
    cbind(0, mat[, -J]) %in% c(0, 2)
  dominantl <- mat == 0 & rbind(mat[-1, ], 1) %in% c(1, 2) & 
    cbind(mat[, -1], 1) %in% c(1, 2)
  dominant <- dominantl | dominantu
  dominant
}



mtc<-function (dominant, admissible, strategy, rec.i, rec.j, pi.theta, 
               mat, p, alternate, psmooth) 
{
  m <- dim(rec.i)[1]
  I <- dim(pi.theta)[1]
  J <- dim(pi.theta)[2]
  k <- I * J
  if (alternate == T) {
    if (sum(dominant & admissible) > 1) {
      if (m == 0) {
        if (any(mat[dominant & admissible] == 0)) 
          dominant[mat == 1] <- FALSE
      }
      else {
        if (mat[rec.i[[m]], rec.j[[m]]] == 1) {
          if (any(mat[dominant & admissible] == 0)) 
            dominant[mat == 1] <- FALSE
        }
        else {
          if (any(mat[dominant & admissible] == 1)) 
            dominant[mat == 0] <- FALSE
        }
      }
    }
  }
  if (strategy == "ss") {
    if(length(pi.theta[dominant & admissible])==0){
      rec.i <- rbind(rec.i, NA)
      rec.j <- rbind(rec.j, NA) 
    }else{
    test <- pi.theta == max(pi.theta[dominant & admissible]) & 
      dominant & admissible
    chosen = ifelse(sum(test) > 1, sample(sum(test), 1), 
                    1)
    rec.i <- rbind(rec.i, row(pi.theta)[test][chosen])
    rec.j <- rbind(rec.j, col(pi.theta)[test][chosen])
    }
  }
  else if (strategy == "ss-random" | strategy == "weighted_mtc") {
    pi.theta[!(dominant & admissible)] = 0
    chosen = sample(k, 1, prob = pi.theta)
    rec.i <- rbind(rec.i, row(pi.theta)[chosen])
    rec.j <- rbind(rec.j, col(pi.theta)[chosen])
  }
  else if (strategy == "p") {
    p[!(dominant & admissible)] = 2
    test <- abs(p - 0.5) == min(abs(p - 0.5)) & dominant & 
      admissible
    chosen = ifelse(sum(test) > 1, sample(sum(test), 1), 
                    1)
    rec.i <- rbind(rec.i, row(pi.theta)[test][chosen])
    rec.j <- rbind(rec.j, col(pi.theta)[test][chosen])
  }
  else if (strategy == "psmooth") {
    psmooth[!(dominant & admissible)] = 2
    test <- abs(psmooth - 0.5) == min(abs(psmooth - 0.5)) & 
      dominant & admissible
    chosen = ifelse(sum(test) > 1, sample(sum(test), 1), 
                    1)
    rec.i <- rbind(rec.i, row(pi.theta)[test][chosen])
    rec.j <- rbind(rec.j, col(pi.theta)[test][chosen])
  }
  else if (strategy == "ssp") {
    if(length(pi.theta[dominant & admissible])==0){
      rec.i <- rbind(rec.i, NA)
      rec.j <- rbind(rec.j, NA) 
    }else{
    test <- pi.theta == max(pi.theta[dominant & admissible]) & 
      dominant & admissible
    p[!(test)] = 2
    test <- abs(p - 0.5) == min(abs(p - 0.5))
    chosen = ifelse(sum(test) > 1, sample(sum(test), 1), 
                    1)
    rec.i <- rbind(rec.i, row(pi.theta)[test][chosen])
    rec.j <- rbind(rec.j, col(pi.theta)[test][chosen])
    }
  }
  return(list(rec.i = rec.i, rec.j = rec.j))
}



failure.time.fn<-function(trueprob, type, Tmax, t.i0){
  if(type=="uniform" | type=="weibull" | type=="pareto"){
    if(type=="uniform"){
      y.i<-rbinom(1, 1, trueprob)
      ifelse(y.i==1, time.i<-runif(1, min=0, max=Tmax), time.i<-Tmax+999)
    }
    if(type=="weibull"){
      time.i<-rweibull(1, shape=4, scale=Tmax/((log(1/(1-trueprob)))^(1/4)))
      ifelse(time.i<=Tmax, y.i<-1, y.i<-0)
    }
    if(type=="pareto"){
      time.i<-rpareto(1, scale=Tmax/10, shape=log(1-trueprob)/log(1/10))
      ifelse(time.i<=Tmax, y.i<-1, y.i<-0)
    }
    return(list(y.i=y.i, time.i=time.i+t.i0))
  }else{stop("type must be either uniform, weibull, or pareto.")}
}



weight.mat.fn<-function(times, dlt.tox, dlt.times, weight.type, Tmax){
  m<-length(times)-1
  current.time<-times[m+1] # have m patients in study, what are outcomes at time when patient m+1 arrives (OK)
  df<-diag(nrow=m+1, ncol=m+1)
  for(i in 1:m){
    df[i,(i+1):dim(df)[2]]<-weight.fn(weight.type=weight.type, dlt.tox.i=dlt.tox[i], dlt.times.i=dlt.times[i], times=times[(i+1):(m+1)], 
                                      t.i0=times[i], Tmax, all.t.i0=times, all.fail.tox=dlt.tox, all.fail.times=dlt.times)
    #dlt.tox[i] = failure status of patient i
    #dlt.times[i] = failure time of patient i
    #times[(i+1):(m+1)] = start times of previous patients after patient i and current patient (m+1)
    #t.i0=times[i] = start time of patient i
    #Tmax = observation window in time units
    #all.t.i0=all start times and current time, regardless of i
    #all.fail.tox = failure statuses of all patients
    #all.fail.times = failure times of all patients
  }
  colnames(df)<-sapply(1:(m+1), function(z) paste0("t=",times[z]))
  df
}



weight.fn<-function(weight.type, dlt.tox.i, dlt.times.i, times, t.i0, Tmax, all.t.i0, all.fail.tox, all.fail.times){
  if(weight.type=="uniform" | weight.type=="adaptive"){
    weight<-rep(NA,length(times))
    if(weight.type=="uniform"){
      for(i in 1:length(weight)){ifelse(dlt.tox.i==1 & dlt.times.i<=times[i], weight[i]<-1,
                                        ifelse(times[i]>=Tmax+t.i0, weight[i]<-0, weight[i]<-1 - ((times[i]-t.i0)/Tmax)))}
    }
    if(weight.type=="adaptive"){
      for(i in 1:length(weight)){ifelse(dlt.tox.i==1 & dlt.times.i<=times[i], weight[i]<-1,
                                        ifelse(times[i]>=Tmax+t.i0, weight[i]<-0, weight[i]<-weight.adaptive.fn(all.t.i0, all.fail.tox, all.fail.times, Tmax, t.i0=t.i0, current.time=times[i])))}
    }
    weight
  }else{stop("weight.type must be either uniform or adaptive")}
}



weight.adaptive.fn<-function(all.t.i0, all.fail.tox, all.fail.times, Tmax, t.i0, current.time){
  # all.t.i0 = all start times for patients 1:m and the current time
  # all.fail.tox = all failure statuses for patients 1:m
  # all.fail.times = all failure times for patients 1:m
  # Tmax = maximum observation window
  # t.i0 = start time for patient i
  # current.time = current time
  index<-which(all.fail.tox==1 & all.fail.times<current.time) ## <= ##
  ifelse(length(index)==0, ftimes<-c(0, Tmax), ftimes<-sort(c(0,all.fail.times[index] - all.t.i0[index],Tmax)))
  z<-length(ftimes)-2
  u<-seq(0,Tmax,length=101)
  weights<-sapply(u[-length(u)], function(x) {kappa<-max(which(ftimes<=x)) - 1
  1 - (1/(z+1))*(kappa + ((x - ftimes[kappa+1]))/(ftimes[kappa+2] - ftimes[kappa+1]))})
  current.std.time<-current.time-t.i0
  time.point<-which.min(abs(current.std.time - u))
  weights<-c(weights,0)
  weights[time.point]
}



###############################################################
#
# R Code for generating results after simulations have been run
#
###############################################################

results.table.fn<-function(obj, exp.rec="exp"){
  if(exp.rec=="exp"){
        tab<- c(100*t(as.matrix(print.tite.pipe.sim(obj, print=FALSE)$exp)),
                                        mean(sapply(1:length(obj$n.sim), function(z) sum(obj$n.sim[[z]]) )),
                                        100*mean((sapply(1:length(obj$n.sim), function(z) sum(obj$r.sim[[z]]) ))/(sapply(1:length(obj$n.sim), function(z) sum(obj$n.sim[[z]]) ))))
  }else{
        stop.store<-length(which(sapply(1:length(obj$n.sim), function(z) obj$n.rpII[[z]])==0 ) )/length(obj$n.sim)
        n.max<-max(sapply(1:length(obj$n.sim), function(z) sum(obj$n.sim[[z]])))
        tab<- c(100*t(as.matrix(print.tite.pipe.sim(obj, print=FALSE)$rec)),
                                        mean(sapply(1:length(obj$n.sim), function(z) obj$n.rpII[[z]] )),
                                        100*stop.store, early.stop.fn(obj, nmax=n.max))
  }
    if(exp.rec=="rec"){digit.index<-c(0,0,0,0,0,1,1,1)}else{digit.index<-c(0,0,0,0,0,1,0)}
  tab<-data.frame(t(sapply(1:length(digit.index), function(z) round(tab[z], digit.index[z]))))
  if(exp.rec == "rec"){
    colnames(tab)<-c("0-15", "15-24", "25-34", "35-44", "45+", "Mean No. MTDCs", "Trials with No MTDCs", "Trials that stopped early (%)")
  }else{
    colnames(tab)<-c("0-15", "15-24", "25-34", "35-44", "45+", "Average Sample Size", "Mean DLTs (%)")
  }
  tab
}


early.stop.fn<-function(obj, nmax){
  L<-length(obj$n.sim)
  round(100*sum(sapply(1:L, function(z) sum(obj$n.sim[[z]])<nmax))/L,1)
}


trial.duration.fn<-function(obj){
  duration<-mean(sapply(1:length(obj$n.sim), function(z) obj$trial.duration.list[[z]]))
  duration
}


print.tite.pipe.sim<-function (x, pi = x$pi, cut.points = c(0, 15, 25, 35, 45.1, 100)/100, 
                               digits = 1, print = TRUE, ...) 
{
  exp <- x$exp
  rec <- x$p.rec
  cuts <- cut(pi, cut.points, right = F)
  exp.table <- sapply(levels(cuts), function(i) {
    sum(exp[cuts == i])
  })
  rec.table <- sapply(levels(cuts), function(i) {
    sum(rec[cuts == i])
  })
  if (print) {
    cat("\n Experimentation percentages by true toxicity: \n")
    print(round(100 * exp.table, digits))
    cat("\n Recommendation percentages by true toxicity: \n")
    print(round(100 * rec.table, digits))
  }
  return(list(exp.table = exp.table, rec.table = rec.table))
}


plot.tite.pipe.sim<-function (x, pi = x$pi, theta = x$theta, plot = "both", ...) 
{
  exp <- x$exp
  rec <- x$p.rec
  I <- dim(x$n.sim[[1]])[1]
  J <- dim(x$n.sim[[1]])[2]
  mat.true <- pi > theta
  x.true <- 1:(I + 1) - 0.5
  y.true <- c(apply(mat.true, 1, function(i) {
    min(which(i == 1), J + 1) - 0.5
  }), 0.5)
  s <- length(x$n.sim)
  if (plot == "exp") {
    df <- data.frame(x = rep(1:I, J), y = rep(1:J, each = I), 
                     z = 100 * c(exp))
    df2 <- data.frame(x = x.true, y = y.true)
    v1 <- ggplot() + geom_tile(aes(x = x, y = y, fill = z), 
                               data = df) + scale_fill_gradient(name = "Experimentation percentages", 
                                                                low = "white", high = "red") + geom_step(aes(x = x, 
                                                                                                             y = y), data = df2, size = 1.5, colour = "green", 
                                                                                                         linetype = 4) + xlab("Drug A level") + ylab("Drug B level")
    print(v1)
  }
  if (plot == "rec") {
    df <- data.frame(x = rep(1:I, J), y = rep(1:J, each = I), 
                     z = 100 * c(rec))
    df2 <- data.frame(x = x.true, y = y.true)
    v1 <- ggplot() + geom_tile(aes(x = x, y = y, fill = z), 
                               data = df) + scale_fill_gradient(name = "Recommendation percentages", 
                                                                low = "white", high = "red") + geom_step(aes(x = x, 
                                                                                                             y = y), data = df2, size = 1.5, colour = "green", 
                                                                                                         linetype = 4) + xlab("Drug A level") + ylab("Drug B level")
    print(v1)
  }
  if (plot == "both") {
    df.exp <- data.frame(x = rep(1:I, J), y = rep(1:J, each = I), 
                         z = 100 * c(exp), type = "Experimentation")
    df.rec <- data.frame(x = rep(1:I, J), y = rep(1:J, each = I), 
                         z = 100 * c(rec), type = "Recommendation")
    df <- rbind(df.exp, df.rec)
    df2 <- data.frame(x = x.true, y = y.true)
    v1 <- ggplot() + geom_tile(aes(x = x, y = y, fill = z), 
                               data = df) + facet_grid(~type) + scale_fill_gradient(name = "Percent", 
                                                                                    low = "white", high = "red") + geom_step(aes(x = x, 
                                                                                                                                 y = y), data = df2, size = 1.5, colour = "green", 
                                                                                                                             linetype = 4) + xlab("Drug A level") + ylab("Drug B level")
    print(v1)
  }
}


#################################
# Setting up a simulation study #
#################################

## The true probability of DLT for a 4x4 grid of dose combinations
## Taken from Mander and Sweeting (2015) [Simulation 2]
pi.true<-matrix(c(4, 8, 12, 16, 10, 14, 18, 22, 16, 20, 24, 28, 22, 26, 30, 34)/100, nrow=4, ncol=4, byrow=F)
# Prior median - assume here prior medians are same as true DLT probabilities
prior.med <- pi.true


############################################################
# Specifications for scenarios:
#
# - Maximum sample size of 40 patients
# - Simulate nsims trials
# - cohort size of one patient
# - TTL = 0.20
# - pi = true DLT probabilities
# - prior median specified, with weight of one patient split over all doses
# - strategy - dose-escalation strategy: "ss" = sample size, "ss-random" = weight randomisation
# - constraint - dose-skipping constraints: "neighbouring" (doses around current) or "no.dose.skip" (doses within one dose of any previously tested)
# - with non-diagonal option of "-nodiag".
# - epsilon - doses not considered admissible if P(dose > MTC)>epsilon
# - admis - "adjacent" or "closest" - adjacent includes doses next to MTC but may be dominated by another safe dose
# - alternate - logical: should design always de-escalate if above the MTD and escalate if below? (subject to other constraints) - is design coherent?
# - upper.tox.constraint - default NULL - number between 0 and 1, so MTC for this number represents upper tox boundary
# - stop - alternative stopping rule - stop trial if P(lowest dose > theta)>stop.
# - non.admissible - logical matrix for doses - determining which ones are permanently barred from trial
# - seed - random seed number
# - Tmax - follow-up time interval length
# - lambda - rate of arrival times (average of one patient arriving every 1/lambda time units, i.e. lambda patients per time unit)
# - first.full.eval - number of patients that have to have completed follow-up before new patients can enter trial
# - wait.to.enter - do all patients have to wait until follow-up has been completed on patients?
# - failure.type - "uniform", "weibull", or "pareto": mechanism to generate failure times
# - weight.type - "uniform", or "adaptive": weight mechanism (see Cheung and Chappell nsims)
#
############################################################

N <- 40
nsims <- 100
cohort.size <- 1
theta <- 0.20
pi <- pi.true
prior.med <- prior.med
prior.ss <- matrix(1/prod(dim(pi)),nrow=dim(pi)[1],ncol=dim(pi)[2])
epsilon <- 0.80


############################
# Current time parameters:
# - Tmax = 1 time unit
# - lambda = 1 (average of one person per time unit recruitment)
# - first.full.eval = min.pats patients to be fully evaluated before new patients can be assigned to treatment
# - wait.to.enter = FALSE for Partial outcome dependent scenario, TRUE for other
# - failure.type = "uniform" - generate binomial outcome, then if DLT occurred, simulate failure time
# - weight.type = "uniform" - (uniform weighting, not DLT time dependent)
# - min.cohort.size = 2 (require minimum of two people per cohort)
#
###############################

Tmax <- 1
lambda <- 1
first.full.eval <- 2
failure.type <- "uniform"
weight.type <- "uniform"
min.cohort.size <- 2


###################
# Run simulations #
###################

# PIPE
trial.1<-tite.pipe.design(N=N, S=nsims, c=cohort.size, theta=theta, pi=pi.true, prior.med=prior.med,
                          prior.ss=prior.ss, strategy="ss", constraint="neighbouring",
                          epsilon=epsilon, admis="closest", alternate=FALSE,
                          Tmax = Tmax, lambda = lambda, first.full.eval = first.full.eval,
                          wait.to.enter = TRUE, failure.type = failure.type, weight.type = weight.type, seed=TRUE, min.cohort.size=min.cohort.size, min.patient.type = "complete.dose")

# TITE-PIPE-C
trial.2<-tite.pipe.design(N=N, S=nsims, c=cohort.size, theta=theta, pi=pi.true, prior.med=prior.med,
                                        prior.ss=prior.ss, strategy="ss", constraint="neighbouring",
                                        epsilon=epsilon, admis="closest", alternate=FALSE,
                                        Tmax = Tmax, lambda = lambda, first.full.eval = first.full.eval,
                                        wait.to.enter = FALSE, failure.type = failure.type, weight.type = weight.type, seed=TRUE, min.cohort.size=min.cohort.size, min.patient.type = "complete.dose")

# TITE-PIPE-O
trial.3<-tite.pipe.design(N=N, S=nsims, c=cohort.size, theta=theta, pi=pi.true, prior.med=prior.med,
                                        prior.ss=prior.ss, strategy="ss", constraint="neighbouring",
                                        epsilon=epsilon, admis="closest", alternate=FALSE,
                                        Tmax = Tmax, lambda = lambda, first.full.eval = first.full.eval,
                                        wait.to.enter = FALSE, failure.type = failure.type, weight.type = weight.type, seed=TRUE, min.cohort.size=min.cohort.size, min.patient.type = "on.dose")

##############################
# How to analyse the results #
##############################

# Analysing set of simulations from "trial.1"
# Experimentation percentages
results.table.fn(obj = trial.1, exp.rec="exp")
# Recommendation percentages
results.table.fn(obj = trial.1, exp.rec="rec")
# Average trial duration
trial.duration.fn(obj = trial.1)
# Get basic experimentation and recommendation with varying toxicity cut-off intervals
print.tite.pipe.sim(trial.1)
# Plot heatmaps for experimentation and recommendation over doses-toxicity surface
plot.tite.pipe.sim(trial.1)


#######
# END #
#######
