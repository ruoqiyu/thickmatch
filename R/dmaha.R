dmaha<-function(z,X,min.control=1,exact=NULL,nearexact=NULL,penalty=1000,rank=FALSE){
  Xmatrix<-function(x){
    if (is.vector(x) || is.factor(x)) x<-matrix(x,nrow=length(z))

    if(is.data.frame(x) || is.character(x)){
      if(!is.data.frame(x)) x <- as.data.frame(x)
      X.chars <- which(plyr::laply(x, function(y) 'character' %in% class(y)))
      if(length(X.chars) > 0){
        for(i in X.chars){
          x[,i] <- factor(x[,i])

        }
      }
      #if some variables are factors convert to dummies
      X.factors <-  which(plyr::laply(x, function(y) 'factor' %in% class(y)))

      #handle missing data
      for(i in which(plyr::laply(x, function(y) any(is.na(y))))){
        if(i %in% X.factors){
          #for factors, make NA a new factor level
          x[,i] <- addNA(x[,i])
        }else{
          #for numeric/logical, impute means and add a new indicator for missingness
          x[[paste(colnames(x)[i],'NA', sep = '')]] <- is.na(x[,i])
          x[which(is.na(x[,i])),i] <- mean(x[,i], na.rm = TRUE)
        }
      }
      for(i in rev(X.factors)){
        dummyXi <- model.matrix(as.formula(
          paste('~',colnames(x)[i], '-1')),data=x)
        x <- cbind(x[,-i], dummyXi)
      }

    }else{
      #handle missing data
      for(i in c(1:ncol(x))){
        if(any(is.na(x[,i]))){
          x <- cbind(x,is.na(X[,i]))
          colnames(x)[ncol(x)] <- paste(colnames(X)[i],'NA', sep = '')
          x[which(is.na(x[,i])),i] <- mean(x[,i], na.rm = TRUE)
        }
      }

    }

    #get rid of columns that do not vary
    varying <- apply(x,2, function(y) length(unique(y)) > 1)
    x <- x[,which(varying),drop = FALSE]

    as.matrix(x)
  }

  X<-Xmatrix(X)
  if (is.vector(X)) X<-as.matrix(X,ncol=1)
  if (is.data.frame(nearexact)) nearexact<-data.matrix(nearexact)
#  if (is.vector(nearexact)) nearexact<-as.matrix(nearexact,ncol=1)
  if (is.factor(nearexact)){
    levels(nearexact)=1:nlevels(nearexact)
    nearexact<-as.numeric(nearexact)
  }
  stopifnot(is.matrix(X))
  stopifnot(is.vector(z))
  stopifnot(all((z==1)|(z==0)))
  stopifnot(length(z)==(dim(X)[1]))
  stopifnot(sum(is.na(X))==0)

  if (!is.null(exact)){
    stopifnot(length(exact)==length(z))
    tb<-table(z,exact)
    if (!all(tb[2,]*min.control<=tb[1,])){
      stop("An exact match for exact is infeasible.")
    }
  }

  n<-dim(X)[1]

  ids<-1:n
  m<-sum(z)

  #X<-data.matrix(Filter(function(x)(length(unique(x))>1), data.frame(X)))
  k<-dim(X)[2]
  if (rank==T){
    for (j in 1:k) X[, j]<-rank(X[, j])
    cv<-stats::cov(X)
    vuntied<-stats::var(1:n)
    rat<-sqrt(vuntied/diag(cv))
    cv<-diag(rat)%*%cv%*%diag(rat)
  }else{
    cv<-stats::cov(X)
  }


  icov<-MASS::ginv(cv)
  #LL<-chol(cv)

  if (is.vector(X)) X<-matrix(X,ncol=1)
  X0<-X[z==0,]
  X1<-X[z==1,]
  if (is.vector(X0)) X0<-matrix(X0,ncol=1)
  if (is.vector(X1)){
    if (sum(z)==1) X1<-matrix(X1,nrow=1)
    else X1<-matrix(X1,ncol=1)
  }

  out <- matrix(NA, m, n - m)
  rownames(out) <- rownames(X)[z == 1]
  colnames(out) <- rownames(X)[z == 0]
  #for (i in 1:m) out[i, ] <- mvnfast::maha(X0,t(as.matrix(X1[i,])),LL,isChol=TRUE)
  for (i in 1:m) out[i, ] <- stats::mahalanobis(X0,t(as.matrix(X1[i,])),icov,inverted = T)
  if (!is.null(exact)){
    dif <- outer(exact[z == 1], exact[z == 0], "!=")
    out[dif] <- Inf
  }

  if (!is.null(nearexact)){
    if (is.vector(nearexact)){
      dif <- outer(nearexact[z == 1], nearexact[z == 0], "!=")
      out <- out + dif * penalty
    }else{
      nnearex=dim(nearexact)[2]
      for (kk in 1:nnearex){
        dif <- outer(nearexact[z == 1,kk], nearexact[z == 0,kk], "!=")
        out <- out + dif * penalty
      }
    }

  }
  d<-t(out)
  dim(d)=c(1,m*(n-m))
  list(dm=out,d=as.vector(d),start=rep(1:m,each=n-m),end=rep((m+1):n,m))
}
