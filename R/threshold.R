threshold<-function(z,X,p,caliper,dat,ncontrol=1,exact=NULL,nearexact=NULL,fine=NULL,penalty=1000,nearexpenalty=100,rank=FALSE,select_num=0,tol=0.1){

  #check input
  stopifnot(is.vector(z))
  stopifnot(is.vector(p))
  stopifnot(length(z)==length(p))
  stopifnot(all((z==1)|(z==0)))
  stopifnot(sum(z)*ncontrol<=sum(1-z))

  if (!is.null(exact)){
    exact<-as.factor(exact)
    nexactlevels<-nlevels(exact)
    levels(exact)<-1:nexactlevels
    tb<-table(z,exact)
    if (!all(tb[2,]*ncontrol<=tb[1,])){
      stop("An exact match for exact is infeasible for every caliper.")
    }
  }

  if (is.null(fine)) fine<-rep(1,length(z))
  # nomissing.ix<-apply(X,1,function(x) sum(is.na(x))==0)
  # z<-z[nomissing.ix]
  # X<-X[nomissing.ix,]
  # p<-p[nomissing.ix]
  # dat<-dat[nomissing.ix,]
  # exact<-exact[nomissing.ix]
  # fine<-fine[nomissing.ix]

  high<-100
  low<-0

  TF0<-feasible(z,X,p,caliper,dat,ncontrol,exact,nearexact,fine,penalty,nearexpenalty,rank,select_num,low)
  if (TF0==1){
    return (list(epsilon=low,interval=c(0,low),interval.length=low))
  }
  TF1<-feasible(z,X,p,caliper,dat,ncontrol,exact,nearexact,fine,penalty,nearexpenalty,rank,select_num,high)
  while (TF1==0){
    low<-high
    high<-2*high
    TF1<-feasible(z,X,p,caliper,dat,ncontrol,exact,nearexact,fine,penalty,nearexpenalty,rank,select_num,high)
  }
  while ((high-low)>tol){
    mid<-(high+low)/2
    TF<-feasible(z,X,p,caliper,dat,ncontrol,exact,nearexact,fine,penalty,nearexpenalty,rank,select_num,mid)
    print(mid)
    print(TF)
    if (TF==0) low<-mid
    else high<-mid
  }
  list(epsilon=high,interval=c(low,high),interval.length=high-low)

}
