threshold_match<-function(z,p,caliper,X,dat,min.control=1,max.control=min.control,total.control=sum(z)*min.control,exact=NULL,fine=rep(1,length(z)),finepenalty=1000,nearexact=NULL,nearexpenalty=100,eps=NULL,penalty=10000,rank=FALSE){

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
  stopifnot(max.control>=min.control)
  stopifnot(length(z)==length(p))
  stopifnot(length(z)==length(fine))
  nobs<-length(z)
  ntreat<-sum(z)
  ncontr<-sum(1-z)
  stopifnot(ncontr>=total.control)
  stopifnot(total.control>=(min.control*ntreat))
  if (is.data.frame(X)){
    X<-as.data.frame(unclass(X))
    X<-data.matrix(X)
  }
  if (is.vector(X)) X<-matrix(X,length(X),1)
  if (is.data.frame(nearexact)) nearexact<-data.matrix(nearexact)
  if (is.factor(nearexact)){
    levels(nearexact)<-1:nlevels(nearexact)
    nearexact<-as.numeric(nearexact)
  }
  if (is.vector(nearexact)) nearexact<-matrix(nearexact,length(nearexact),1)
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
  dist<-dmaha(z,X,min.control,exact,nearexact,nearexpenalty,rank)
  dist<-DiPs::addcaliper(dist, z, p, c(-caliper,caliper), stdev = TRUE, penalty)
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
    if (!is.null(eps)) m<-list(feasible=output$feasible,d=dat1,Xm=X[matchid,],zm=z[matchid],Xs=X[sid,],zs=z[sid],zms=z[smatchid],dms=dat[smatchid,],Xms=X[smatchid,],number=net$tcarcs)
    else m<-list(feasible=output$feasible,d=dat1,Xm=X[matchid,],zm=z[matchid],Xs=X,zs=z,zms=z[matchid],dms=dat[matchid,],Xms=X[matchid,],number=net$tcarcs)
  }
  if(m[[1]]==0) {
    if (!is.null(exact)) warning("The match you requested is infeasible.  Reconsider caliper, ncontrol and exact.")
    else warning("The match you requested is infeasible, perhaps because the caliper is too small.")
  }else{
    balance=DiPs::check(X,m$Xm,z,m$zm)
    sbalance=DiPs::check(m$Xs,m$Xms,m$zs,m$zms)
    if (!is.null(eps)) list(data=m$d,sdata=m$dms,balance=balance,sbalance=sbalance,closest=distance)
    else list(data=m$d,sdata=m$dms,balance=balance,sbalance=sbalance,closest=NULL)
  }
}
