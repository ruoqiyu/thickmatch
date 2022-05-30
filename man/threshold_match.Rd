\name{threshold_match}
\alias{threshold_match}
\title{Minimum-distance threshold matching.}
\usage{
threshold_match(z,p,caliper,X,dat,min.control=1,
max.control=min.control,total.control=sum(z)*min.control,
exact=NULL,fine=rep(1,length(z)),finepenalty=1000,nearexact=NULL,
nearexpenalty=100,eps=NULL,penalty=10000,rank=FALSE)
}
\arguments{
\item{z}{ A vector whose ith coordinate is 1 for a treated unit and is 0 for a control.}
\item{p}{ A vector of with length(z)=length(p) giving the variable used to define the caliper.  For instance,
p might be the propensity score.}
\item{caliper}{ If the treated-minus-control
difference (in the scale of sd(p)) in p is < -caliper or > caliper, then penalty is added to the distance.}
\item{X}{ A matrix with length(z) rows giving the covariates.  X should be of full column rank.}
\item{dat}{ A data frame with length(z) rows.  If the match is feasible, the matched portion of dat is returned with additional columns that define the match.}
\item{min.control}{ A positive integer giving the minimum number of controls to be matched to each treated subject.  If min.control is too large, the match will be infeasible.}
\item{max.control}{ A positive integer giving the maximum number of controls to be matched to each treated subject.}
\item{total.control}{ A positive integer giving the total number of controls to be matched to each treated subject.  If total.control is too large, the match will be infeasible. Fine balance constraint can be determined based on total.control.}
\item{exact}{ If not NULL, then a vector of length(z)=length(p) giving variable that need to be exactly matched. }
\item{fine}{ A vector of with length(z)=length(fine) giving the nominal levels that are to be nearly-finely balanced.}
\item{finepenalty}{ A numeric penalty imposed for each violation of fine balance.}
\item{nearexact}{ If not NULL, then a vector of length length(z) giving variable that need to be exactly matched. }
\item{nearexpenalty}{ The penalty for a mismatch on nearexact.}
\item{eps}{ The threshold whose feasibility is examined. If eps is NULL, the conventional optimal match with the propensity score caliper, fine balance, exact and near-exact match constraints is returned.}
\item{penalty}{ A numeric penalty imposed for each distance greater than eps.}
\item{rank}{ If rank=TRUE, a rank-based Mahalanobis distance will be  calculated. Otherwise (with default value FALSE), a traditional Mahalanobis distance will be computed.}
}

\value{
If the match is infeasible, a warning is issued.  Otherwise, a list of results is returned.

A match may be infeasible if min.control or total.control is too large, or if exact matching for exact is impossible.

\item{data}{ The matched sample, selected rows of dat.}
\item{sdata}{ The matched closest pairs, selected rows of dat.}
\item{balance}{ Balance table of the matched sample, including 5 columns: treated mean, matched control mean, all control mean, matched SMD and all SMD. }
\item{sbalance}{ Balance table of the matched closest pairs, including 5 columns: treated mean, matched control mean, all control mean, matched SMD and all SMD. }
}

\description{
The program finds an optimal threshold match with a given threhold on distance, plus near-fine balance, exact match and near-exact match constraints. That is, it finds a match that minimizes the penalized Mahalanobis distance.
}


\details{
The match minimizes the total distance between treated subjects and their matched controls subject to a threshold which imposes a penalty on distances above the threshold.

For discussion of the choice of threshold, see Rosenbaum (2017).

You MUST install and load the optmatch package to use threshold_match.
}
\examples{
# To run this example, you must load the optmatch package.
\donttest{
data("nysr")
attach(nysr)
X<-cbind(family.income,family.structure,highest.education.parent.in.household,
female,race.black,race.hispanic,age.teenager,school.dropout)
detach(nysr)

eps=threshold(z=nysr$intense,X=X,p=nysr$plogit,caliper=0.2,
dat=nysr,select_num=10,tol=0.00001)$epsilon
m<-threshold_match(z=nysr$intense,p=nysr$plogit,caliper=0.2,X=X,dat=nysr,min.control=2,eps=eps)
dim(m$sdata)
m$sbalance
m$balance
}
}
\references{
Rosenbaum, P.R. (2017) Imposing Minimax and Quantile Constraints on Optimal Matching in Observational Studies, Journal of Computational and Graphical Statistics, 26:1, 66-78, DOI: 10.1080/10618600.2016.1152971.
}
