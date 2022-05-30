\name{dmaha}
\alias{dmaha}
\title{Creates a Mahalanobis distance for matching based on a dense network.}
\usage{
dmaha(z, X, min.control=1,exact=NULL, nearexact=NULL,
penalty=1000, rank=FALSE)
}
\arguments{
\item{z}{ A vector whose ith coordinate is 1 for a treated unit and is 0 for a control. }
\item{X}{ A matrix with length(z) rows giving the covariates.  X should be of full column rank.}
\item{min.control}{ A positive integer giving the minimum number of controls to be matched to each treated subject.  If min.control is too large, the match will be infeasible.}
\item{exact}{ If not NULL, then a vector of length(z)=length(p) giving variable that need to be exactly matched. }
\item{nearexact}{ If not NULL, then a vector of length length(z) giving variable that need to be exactly matched. }
\item{penalty}{ The penalty for a mismatch on nearexact.}
\item{rank}{ If rank=TRUE, a rank-based Mahalanobis distance will be  calculated. Otherwise (with default value FALSE), a traditional Mahalanobis distance will be computed.}
}

\value{
\item{d}{ A distance object for each pair of treated and control.}
\item{start}{ The treated subject for each distance.}
\item{end}{ The control subject for each distance.}
}
\description{
Computes a Mahalanobis distance list,either the traditional version or the rank-based version, for use in dense matching, i.e. the distance for all possible pairs of treated and control.

This function and its use are discussed in Rosenbaum (2010).  The rank-based Mahalanobis distance in described in Chapter 8 of Rosenbaum (2010).

}
\details{
The usual Mahalanobis distance works well for multivariate Normal covariates, but can exhibit odd behavior with typical covariates. Long tails or an outlier in a covariate can yield a large estimated variance, so the usual Mahalanobis distance pays little attention to large differences in this covariate. Rare binary covariates have a small variance, so a mismatch on a rare binary covariate is viewed by the usual Mahalanobis distance as extremely important. If you were matching for binary covariates indicating US state of residence, the usual Mahalanobis distance would regard a mismatch for Wyoming as much worse than a mismatch for California.

The robust Mahalanobis distance uses ranks of covariates rather than the covariates themselves, but the variances of the ranks are not adjusted for ties, so ties do not make a variable more important. Binary covariates are, of course, heavily tied.
}
\examples{
data("nysr")
attach(nysr)
X<-cbind(family.income,family.structure,highest.education.parent.in.household,
         female,race.black,race.hispanic,age.teenager,school.dropout)
dist<-dmaha(intense,X)
head(dist$d)
detach(nysr)

}
\references{

Rosenbaum, P. R. (2010) Design of Observational Studies.  New York: Springer.

}
