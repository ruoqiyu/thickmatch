feasible<-function(dist,z,dat,ncontrol=1,fine=rep(1,length(z)),penalty=1000,select_num=0,eps=1000){
  
  ntreat=sum(z)
  #do match
  dist$d<-as.numeric(dist$d>eps)
  
  if (!requireNamespace("optmatch",quietly=TRUE)){
    stop("Error: package optmatch (>= 0.9-1) not loaded.  To run rcbalance command, you must install optmatch first and agree to the terms of its license.")
  }
  
  net<-netvr(z,dist,ncontrol,ncontrol,ntreat*ncontrol,fine,penalty)
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

threshold<-function(dist,z,dat,ncontrol=1,fine=NULL,penalty=1000,select_num=0,tol=0.1){
  
  #check input
  stopifnot(is.vector(z))
  stopifnot(all((z==1)|(z==0)))
  stopifnot(sum(z)*ncontrol<=sum(1-z))
  
  if (is.null(fine)) fine<-rep(1,length(z))
  # nomissing.ix<-apply(X,1,function(x) sum(is.na(x))==0)
  # z<-z[nomissing.ix]
  # X<-X[nomissing.ix,]
  # p<-p[nomissing.ix]
  # dat<-dat[nomissing.ix,]
  # exact<-exact[nomissing.ix]
  # fine<-fine[nomissing.ix]
  
  high<-1
  low<-0
  
  TF0<-feasible(dist,z,dat,ncontrol,fine,penalty,select_num,low)
  if (TF0==1){
    return (list(epsilon=low,interval=c(0,low),interval.length=low))
  }
  TF1<-feasible(dist,z,dat,ncontrol,fine,penalty,select_num,high)
  while (TF1==0){
    low<-high
    high<-2*high
    TF1<-feasible(dist,z,dat,ncontrol,fine,penalty,select_num,high)
  }
  while ((high-low)>tol){
    mid<-(high+low)/2
    TF<-feasible(dist,z,dat,ncontrol,fine,penalty,select_num,mid)
    print(mid)
    print(TF)
    if (TF==0) low<-mid
    else high<-mid
  }
  list(epsilon=high,interval=c(low,high),interval.length=high-low)
  
}

threshold_match<-function(dist,z,dat,min.control=1,max.control=min.control,total.control=sum(z)*min.control,fine=rep(1,length(z)),finepenalty=1000,eps=NULL,penalty=10000){
  
  ntreat=sum(z)
  #do match
  if (!is.null(eps)) dist$d[which(dist$d>eps)]<-10*penalty+dist$d[which(dist$d>eps)]
  dist$d[which(dist$d==Inf)]=2*max(dist$d[which(dist$d!=Inf)])
  if (!requireNamespace("optmatch",quietly=TRUE)){
    stop("Error: package optmatch (>= 0.9-1) not loaded.  To run rcbalance command, you must install optmatch first and agree to the terms of its license.")
  }
  
  net<-netvr(z,dist,min.control,max.control,total.control,fine,finepenalty)
  
  output<-rcbalance::callrelax(net)
  if (output$feasible!=1){
    warning("Match is infeasible.  Change dist or ncontrol to obtain a feasible match.")
    m<-list(feasible=output$feasible,d=NULL)
  }else{
    x<-output$x[(2*ntreat+2):(2*ntreat+1+net$tcarcs)]
    treated<-net$startn[(2*ntreat+2):(2*ntreat+1+net$tcarcs)]-2
    control<-net$endn[(2*ntreat+2):(2*ntreat+1+net$tcarcs)]-2
    newd<-dist$d
    newd[which(x!=1)]=10*max(newd)
    
    n<-length(z)
    ntreat<-sum(z)
    id1<-(1:n)[z==1]
    id0<-(1:n)[z==0]
    
    if (!is.null(eps)){
      select_index<-which(newd<=eps)
      select_treated<-treated[select_index]
      select_control<-control[select_index]
      distance<-sort(newd)[1:length(select_index)]
      #smatchid<-matrix(c(id1[as.numeric(row.names(smatches))],id0[as.vector((smatches-sum(z)))]),ncol=ncontrol+1)
      smatchid<-matrix(c(id1[select_treated],id0[(select_control-sum(z))]),ncol=2)
      smatchid<-as.vector(t(smatchid))
      sid<-c(id1[select_treated],id0)
    }
    
    match.df<-data.frame('treat'=treated,'x'=x,'control'=control)
    matched.or.not<-plyr::daply(match.df,plyr::.(match.df$treat),
                                function(treat.edges) c(as.numeric(as.character(treat.edges$treat[1])),sum(treat.edges$x)),.drop_o=FALSE)
    if(any(matched.or.not[,2]==0)){
      match.df<-match.df[-which(match.df$treat %in% matched.or.not[which(matched.or.not[,2]==0),1]),]
    }
    #match.df$treat<-as.factor(as.character(match.df$treat))
    #matches<-as.matrix(plyr::daply(match.df, plyr::.(match.df$treat),
    #                               function(treat.edges) treat.edges$control[treat.edges$x==1],.drop_o=FALSE))
    #matchid<-matrix(c(id1[as.numeric(row.names(matches))],id0[as.vector((matches-sum(z)))]),ncol=ncontrol+1)
    #matchid<-as.vector(t(matchid))
    matches<-match.df[match.df$x==1,]
    ixs<-cbind(id1[matches[,1]],id0[matches[,3]-ntreat])
    ixs=as.vector(t(ixs))
    matchid<-ixs[!duplicated(ixs)]
    
    #restid<-setdiff(1:n,id0[(select_control-sum(z))])
    dat1<-dat[matchid,]
    seti=0
    mset=c()
    for (mi in 1:length(matchid)){
      if ((z[matchid])[mi]==1){
        seti=seti+1
      }
      mset=c(mset,seti)
    }
    
    dat1<-cbind(mset,dat1)
    if (!is.null(eps)) m<-list(feasible=output$feasible,d=dat1,zm=z[matchid],dms=dat[smatchid,],number=net$tcarcs)
    else m<-list(feasible=output$feasible,d=dat1,zm=z[matchid],dms=dat[matchid,],number=net$tcarcs)
  }
  if(m[[1]]==0) {
    warning("The match you requested is infeasible, perhaps because the caliper is too small.")
  }else{
    if (!is.null(eps)) list(data=m$d,sdata=m$dms,closest=distance)
    else list(data=m$d,sdata=m$dms,closest=NULL)
  }
}
