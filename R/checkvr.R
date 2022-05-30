checkvr=function(xf,x,zf,z,mset){
  d=ncol(x)
  nostratum=length(unique(mset))
  res=matrix(0,nrow=d,ncol=5)
  rownames(res)=colnames(x)
  colnames(res)=c('Treated Mean','Control Match Mean','Control All Mean', 'Control Match SMD','Control All SMD')
  for (k in 1:d){
    xk=x[,k]
    var.xtreated=var(xf[zf==1,k]);
    var.xcontrol=var(xf[zf==0,k]);
    combinedsd=sqrt(.5*(var.xtreated+var.xcontrol));
    std.diff.before.matching=(mean(xf[zf==1,k])-mean(xf[zf==0,k]))/combinedsd;
    diff.in.stratum=rep(0,nostratum);
    treated.in.stratum=rep(0,nostratum);
    meancontrol.in.stratum=rep(0,nostratum)
    for(i in 1:nostratum){
      diff.in.stratum[i]=mean(xk[mset==i & z==1])-mean(xk[mset==i & z==0]);
      treated.in.stratum[i]=sum(mset==i & z==1);
      meancontrol.in.stratum[i]=mean(xk[mset==i & z==0])
    }
    std.diff.after.matching=(sum(treated.in.stratum*diff.in.stratum)/sum(treated.in.stratum))/combinedsd;
    mean.control.after=sum(meancontrol.in.stratum*treated.in.stratum)/sum(treated.in.stratum)
    res[k,]=c(mean(xk[z==1]),mean.control.after,mean(xf[zf==0,k]),std.diff.after.matching,std.diff.before.matching)
  }
  res
}
