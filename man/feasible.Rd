\name{feasible}
\alias{feasible}
\title{Feasilibity of a chosen threshold.}
\usage{
feasible(z,X,p,caliper,dat,ncontrol=1,exact=NULL,
nearexact=NULL,fine=rep(1,length(z)),penalty=1000,nearexpenalty=100,
rank=FALSE,select_num=0,eps=1000)
}
\arguments{
\item{z}{ A vector whose ith coordinate is 1 for a treated unit and is 0 for a control.}
\item{X}{ A matrix with length(z) rows giving the covariates.  X should be of full column rank.}
\item{p}{ A vector of with length(z)=length(p) giving the variable used to define the caliper.  For instance,
p might be the propensity score.}
\item{caliper}{ If the treated-minus-control difference (in the scale of sd(p)) in p is < -caliper or > caliper, then penalty is added to the distance.}
\item{dat}{ A data frame with length(z) rows.  If the match is feasible, the matched portion of dat is returned with additional columns that define the match.}
\item{ncontrol}{ A positive integer giving the number of controls to be matched to each treated subject.  If ncontrol is too large, the match will be infeasible.}
\item{exact}{ If not NULL, then a vector of length(z)=length(p) giving variable that need to be exactly matched. }
\item{nearexact}{ If not NULL, then a vector of length length(z) giving variable that need to be exactly matched. }
\item{fine}{ A vector of with length(z)=length(fine) giving the nominal levels that are to be nearly-finely balanced.}
\item{penalty}{ A numeric penalty imposed for each violation of fine balance.}
\item{nearexpenalty}{ The penalty for a mismatch on nearexact.}
\item{rank}{ If rank=TRUE, a rank-based Mahalanobis distance will be  calculated. Otherwise (with default value FALSE), a traditional Mahalanobis distance will be computed.}
\item{select_num}{ A positive number giving the required number of matched pairs with distance not exceeding eps.}
\item{eps}{ The threshold whose feasibility is examined.}
}

\value{
If the match is infeasible, a warning is issued.  Otherwise, a binary indicator of whether there are at least select_num matched pairs with distance not exceeding eps.

A match may be infeasible if ncontrol is too large, or if exact matching for exact is impossible.
}

\description{
The program determines whether it is possible to find at least select_num matched pairs with distance not exceeding eps.
}


\details{
If there is a feasible matching with at least select_num matched pairs with distance not exceeding eps, then eps is said to be feasible and 1 is returned. Otherwise, eps is said to be infeasible and 0 is returned.

For details, see Rosenbaum (2017).

You MUST install and load the optmatch package to use feasible.
}
\examples{
# To run this example, you must load the optmatch package.
data("nysr")
attach(nysr)
X<-cbind(family.income,family.structure,highest.education.parent.in.household,
         female,race.black,race.hispanic,age.teenager,school.dropout)
detach(nysr)
feasible(z=nysr$intense,X=X,p=nysr$plogit,caliper=0.2,dat=nysr,select_num=10,eps=0.5)
}
\references{
Rosenbaum, P.R. (2017) Imposing Minimax and Quantile Constraints on Optimal Matching in Observational Studies, Journal of Computational and Graphical Statistics, 26:1, 66-78,
DOI: 10.1080/10618600.2016.1152971.
}
