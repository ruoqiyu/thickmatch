netvr<-function(z,dist,min.control=1,max.control=min.control,total.control=sum(z)*min.control,fine=rep(1,length(z)),penalty=1000){

  #check input
  stopifnot(is.vector(z))
  stopifnot(is.vector(fine))
  fine<-as.numeric(fine)
  stopifnot(all((z==1)|(z==0)))
  stopifnot(max.control>=min.control)
  nobs<-length(z)
  ntreat<-sum(z)
  ncontr<-sum(1-z)
  stopifnot(nobs==length(fine))
  stopifnot(total.control>=ntreat*min.control)
  stopifnot(total.control<=ncontr)

  #create basic treated-vs-control bipartite graph
  fine1<-fine[z==1]
  fine0<-fine[z==0]

  #source and extras
  startn<-rep(1,ntreat+1)
  endn<-2:(ntreat+2)
  cost<-rep(0,ntreat+1)
  ucap<-c(total.control-min.control*ntreat,rep(min.control, ntreat))
  b<-c(total.control,rep(0,ntreat+1)) #supply for source

  #connect extras to treated
  startn<-c(startn,rep(2,ntreat))
  endn<-c(endn,3:(ntreat+2))
  cost<-c(cost,rep(0,ntreat))
  ucap<-c(ucap,rep(max.control-min.control, ntreat))

  #connect treated to controls
  startn<-c(startn,dist$start+2)
  endn<-c(endn,dist$end+2)
  cost<-c(cost,dist$d)
  tcarcs<-length(dist$start) # number of treatment-control arcs
  ucap<-c(ucap,rep(1,tcarcs))
  b<-c(b,rep(0,ncontr)) #supply for control nodes
  #Make costs integer
  cost<-round(100*cost)

  #create a duplicate for each control to make sure each control is only used once
  startn<-c(startn,(ntreat+3):(nobs+2))
  endn<-c(endn,(nobs+3):(nobs+ncontr+2))
  cost<-c(cost,rep(0,ncontr))
  ucap<-c(ucap,rep(1,ncontr))
  b<-c(b,rep(0,ncontr))

  #Add structure to the bipartite graph for near fine balance

  tb<-table(z,fine)
  nc<-as.vector(tb[1,])
  nt<-as.vector(tb[2,])
  nwant<-ceiling(nt*total.control/ntreat) #desired number
  nlow<-pmin(nc,nwant) #available number
  nextra<-nwant-nlow #gap between desired and available
  finelevels<-as.vector(as.numeric(colnames(tb)))

  #Add a node for fine balance category k
  sinks<-NULL
  for (k in 1:length(nlow)){
    if(nlow[k]>0){
      sinkk<-length(b)+1
      sinks<-c(sinks,sinkk)
      who0<-fine0==finelevels[k]
      b<-c(b,0)
      if (sum(who0)>0){
        startn<-c(startn,rep(nobs+2,sum(who0))+which(who0))
        endn<-c(endn,rep(sinkk,sum(who0)))
        ucap<-c(ucap,rep(1,sum(who0)))
        cost<-c(cost,rep(0,sum(who0)))
      }
    }
  }

  #Add a node to take the extras
  sinkex<-length(b)+1
  b<-c(b,0)
  startn<-c(startn,(nobs+3):(nobs+2+ncontr))
  endn<-c(endn,rep(sinkex,ncontr))
  ucap<-c(ucap,rep(1,ncontr))
  cost<-c(cost,rep(0,ncontr))

  #Add a sink
  finalsink<-length(b)+1
  b<-c(b,-total.control) #finalsink absorbs all flow
  #Connect balance nodes to finalsink
  startn<-c(startn,sinks)
  endn<-c(endn,rep(finalsink,length(sinks)))
  ucap<-c(ucap,nlow[nlow>0])
  cost<-c(cost,rep(0,length(sinks)))

  #Connect sinkex to finalsink
  startn<-c(startn,sinkex)
  endn<-c(endn,finalsink)
  ucap<-c(ucap,ntreat*max.control)
  cost<-c(cost,penalty)

  #Make costs integer
  cost<-round(cost)
  net<-list(startn=startn,endn=endn,ucap=ucap,b=b,cost=cost,tcarcs=tcarcs)
  net
}
