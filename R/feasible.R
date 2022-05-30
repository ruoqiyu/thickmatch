feasible<-function(z,X,p,caliper,dat,ncontrol=1,exact=NULL,nearexact=NULL,fine=rep(1,length(z)),penalty=1000,nearexpenalty=100,rank=FALSE,select_num=0,eps=1000){

  #Check input
  stopifnot(is.data.frame(dat))
  stopifnot(is.vector(z))
  stopifnot(is.vector(p))
  if (is.factor(fine)){
    levels(fine)<-1:nlevels(fine)
    fine<-as.integer(fine)
  }
  stopifnot(is.vector(fine))
  fine<-as.numeric(fine)
  stopifnot(all((z==1)|(z==0)))
  stopifnot((ncontrol==round(ncontrol))&(ncontrol>=1))
  stopifnot(length(z)==length(p))
  stopifnot(length(z)==length(fine))
  nobs<-length(z)
  ntreat<-sum(z)
  ncontr<-sum(1-z)
  stopifnot(ncontr>=(ncontrol*ntreat))
  if (is.data.frame(X)){
    X<-as.data.frame(unclass(X))
    X<-data.matrix(X)
  }
  if (is.vector(X)) X<-matrix(X,length(X),1)
  stopifnot(length(z)==(dim(X)[1]))
  stopifnot(length(z)==(dim(dat)[1]))
  stopifnot(caliper>=0)

  if (!is.null(exact)){
    if (is.factor(exact)){
      levels(exact)<-1:nlevels(exact)
      exact<-as.integer(exact)
    }
  }

  #handle missing data
  for(i in c(1:ncol(X))){
    if(any(is.na(X[,i]))){
      X <- cbind(X,is.na(X[,i]))
      colnames(X)[ncol(X)] <- paste(colnames(X)[i],'NA', sep = '')
      X[which(is.na(X[,i])),i] <- mean(X[,i], na.rm = TRUE)
    }
  }

  #do match
  dist<-dmaha(z,X,ncontrol,exact,nearexact,nearexpenalty,rank)
  dist<-DiPs::addcaliper(dist, z, p, c(-caliper,caliper), stdev = TRUE, penalty)
  dist$d<-as.numeric(dist$d>eps)

  if (!requireNamespace("optmatch",quietly=TRUE)){
    stop("Error: package optmatch (>= 0.9-1) not loaded.  To run rcbalance command, you must install optmatch first and agree to the terms of its license.")
  }

  net<-netvr(z,dist,ncontrol,ncontrol,sum(z)*ncontrol,fine,penalty)
  #net<-netfine(z,dist,ncontrol,fine,penalty)
  output<-rcbalance::callrelax(net)

  if (output$feasible!=1){
    stop("Match is infeasible.  Change dist or ncontrol to obtain a feasible match.")
  }else{
    x<-output$x[(2*ntreat+2):(2*ntreat+1+net$tcarcs)]
    #x<-output$x[1:net$tcarcs]
    if (sum(dist$d[which(x==1)])<=sum(z)-select_num) return (1)
    else return (0)
  }
}
