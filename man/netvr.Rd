\name{netvr}
\alias{netvr}
\title{Optimal variable-ratio match from a distance matrix.}
\usage{
netvr(z,dist,min.control=1,max.control=min.control,
total.control=sum(z)*min.control,
fine=rep(1,length(z)),penalty=1000)
}
\arguments{
\item{z}{ A vector whose ith coordinate is 1 for a treated unit and is 0 for a control. }
\item{dist}{ A distance list with the starting node (treated subjec), ending node (control), the extra distance between them.}
\item{min.control}{ A positive integer giving the minimum number of controls to be matched to each treated subject.  If min.control is too large, the match will be infeasible.}
\item{max.control}{ A positive integer giving the maximum number of controls to be matched to each treated subject.}
\item{total.control}{ A positive integer giving the total number of controls to be matched to each treated subject.  If total.control is too large, the match will be infeasible.}
\item{fine}{ A vector of with length(z)=length(fine) giving the nominal levels that are to be nearly-finely balanced.}
\item{penalty}{ A numeric penalty imposed for each violation of fine balance.}
}
\value{ A network for optimal variable-ratio matching.}
\description{
The function creates the network for optimal variable-ratio matching to be passed via callrelax to the Fortran code for Bertsekas and Tseng's (1988) Relax IV.

Of limited interest to most users; function netvr() would typically be called by some other functions.
}
\details{
The network contains a bipartite graph for treated and control subjects plus additional nodes for fine balance categories, plus additional nodes accept needed deviations from fine balance yielding near-fine balance.
}
\references{
Bertsekas, D. P. and Tseng, P. (1988) The relax codes for linear minimum cost network flow problems. Annals of Operations Research, 13, 125-190.
Fortran and C code: http://www.mit.edu/~dimitrib/home.html.  Available in R via the optmatch package.
}
