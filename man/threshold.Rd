\name{threshold}
\alias{threshold}
\title{Smallest threshold for thick description.}
\usage{
threshold(z,X,p,caliper,dat,ncontrol=1,exact=NULL,nearexact=NULL,fine=NULL,
penalty=1000,nearexpenalty=100,rank=FALSE,select_num=0,tol=0.1)
}
\arguments{
\item{z}{ A vector whose ith coordinate is 1 for a treated unit and is 0 for a control.}
\item{X}{ A matrix with length(z) rows giving the covariates.  X should be of full column rank.}
\item{p}{ A vector of with length(z)=length(p) giving the variable used to define the caliper.  For instance,
p might be the propensity score.}
\item{caliper}{ If the treated-minus-control
difference (in the scale of sd(p)) in p is < -caliper or > caliper, then penalty is added to the distance.}
\item{dat}{ A data frame with length(z) rows.  If the match is feasible, the matched portion of dat is returned with additional columns that define the match.}
\item{ncontrol}{ A positive integer giving the number of controls to be matched to each treated subject.  If ncontrol is too large, the match will be infeasible.}
\item{exact}{ If not NULL, then a vector of length(z)=length(p) giving variable that need to be exactly matched. }
\item{nearexact}{ If not NULL, then a vector of length length(z) giving variable that need to be exactly matched. }
\item{fine}{ A vector of with length(z)=length(fine) giving the nominal levels that are to be nearly-finely balanced.}
\item{penalty}{ A numeric penalty imposed for each violation of fine balance.}
\item{nearexpenalty}{ The penalty for a mismatch on nearexact.}
\item{rank}{ If rank=TRUE, a rank-based Mahalanobis distance will be  calculated. Otherwise (with default value FALSE), a traditional Mahalanobis distance will be computed.}
\item{select_num}{ A positive number giving the required number of matched pairs with distance not exceeding eps.}
\item{tol}{ The tolerance.  The smallest threshold is determined with an error of at most tol.}
}

\value{
If the match is infeasible, a warning is issued.  Otherwise, a list of results is returned.

A match may be infeasible if the caliper on p is too small, or ncontrol is too large, or if exact matching for exact is impossible.

\item{epsilon}{The smallest threshold, with an error of at most tol.  This threshold is a little too large, at most tol too large, but because its error is on the high side, a match with this threshold ensures at least select_num matched pairs with distance not exceeding epsilon.}
\item{interval}{An interval that contains the best threshold. The upper bound of the interval was returned as epsilon above.}
\item{interval.length}{The length of interval. By definition, length.interval<=tol.}
}

\description{
Finds the smallest threshold on such that a treated-control matching with that at least select_num matched pairs having distance not exceeding the threshold exists.
}


\details{
The method uses binary search to find the best threshold. It applies threshold algorithm with function feasible; details see Rosenbaum (2017).

Often, we need a small and feasible threshold, and we prefer to estimate the threshold very precisely.  Making tol smaller makes the number of closest pairs close to select_num.

You MUST install and load the optmatch package to use threshold.
}
\examples{
# To run this example, you must load the optmatch package.
\donttest{
data("nysr")
attach(nysr)
X<-cbind(family.income,family.structure,highest.education.parent.in.household,
female,race.black,race.hispanic,age.teenager,school.dropout)
detach(nysr)
threshold(z=nysr$intense,X=X,p=nysr$plogit,caliper=0.2,dat=nysr,select_num=10,tol=0.00001)
}
}
\references{
Rosenbaum, P.R. (2017) Imposing Minimax and Quantile Constraints on Optimal Matching in Observational Studies, Journal of Computational and Graphical Statistics, 26:1, 66-78, DOI: 10.1080/10618600.2016.1152971.
}
