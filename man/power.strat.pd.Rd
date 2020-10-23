\name{power.strat.pd}
\alias{power.strat.pd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate attained powers and construct power distribution for stratified CRT
}
\description{
Use simulations or analytic formula to calculate attained powers and construct pre-randomization power distribution for all allocations that a randomization algorithm can produce in stratified CRT
}
\usage{
power.strat.pd(I, P, K, S, mu0, Tx.effect, Time.effect = NULL, pwr.thrd = NULL, factor.time = TRUE,  family, design, rho = NULL, sigma.a = NULL, sigma.e = NULL, sig.level = 0.05, method = "analytic", plot = FALSE, gen.all = TRUE, n.allocs = NULL, n.sims, seed = NULL, rep = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{I}{
A vector indicating the number of clusters in each size group (i.e., 8 clusters with six small and two large clusters has I = c(6, 2))
}
  \item{P}{
P is a list. Each element within the list consists of a stratum name and a vector specifying the number of clusters transit at each period within that stratum. (i.e., P=list("a"=c(1,2,1),"b"=c(1,1,2)) means that there are two strata named “a” and “b”. In stratum a, one cluster transit at the first and third transition period, and two clusters transit at the second transition period. In stratum b, one cluster transit at the first and second transition period, and two clusters transit at the third transition period)
}
  \item{K}{
For SW,  K is the cluster size per period for each cluster. Then the total cluster size is K× (number of transition steps+1).
For PCRT, K is the cluster size.
Length of K equals sum(I) - total number of clusters. K should be entered as the order in I (i.e., I = c(6, 2) and K = c(6, 6, 6, 6, 6, 6, 12, 12) means the first six small clusters have 6 pts per period and the two large clusters have 12 pts per period).
K can be a non-integer. If the resulting cluster size is also a non-integer, the cluster size is determined by randomly picking one of the two nearest integers with probabilities that ensure the non-integer cluster size is achieved in expectation (e.g., a cluster has K = 6.1 with four transition periods, the total cluster size is 30.5. Therefore, cluster size for this cluster will be either 30 or 31 with probabilities 0.9 an 0.1)
}
  \item{S}{
Stratum ID for each cluster. S is a factor variable and should be entered as the same order as I. (i.e., I = c(6, 2) and S = c(“a”, “b”, “b”, “a”, “a”, “b”, “b”, “a”) means among six small clusters, small cluster #1 is at stratum a, small cluster #2 is at stratum b, etc., and the two large clusters are at stratum a and b respectively. Note that the levels in S should match the name of elements in P
}
  \item{mu0}{
Baseline mean/risk/rate on the outcome’s natural scale (i.e., not transformed by the link function)
}
  \item{Tx.effect}{
Treatment effect (natural scale: risk difference gaussian outcome, odds ratio for binary outcome, risk ratio for poisson outcome)
}
  \item{Time.effect}{
Time effect (linear scale). Default is no time effect. The time effect can be a continuous linear trend (a single number) or categorical time effect
}
  \item{pwr.thrd}{
The threshold used to define “low power”. If NULL no risk will be calculated, risk in the output returns NULL by default
}
  \item{factor.time}{
Logical value. Whether treating time as a categorical variable in analysis model. In simulation, user can treat time differently in data generator and analysis model. (i.e., continuous time effect in data generator, but categorical time effect in the analysis). However, these two must match in analytic calculation. In P-CRT, the value is ignored.
}
  \item{family}{
Indicate the type of outcome. One of “gaussian”, “binomial” or “poisson” (case sensitive)
}
  \item{design}{
Choice of parallel cluster randomization trial (“pcrt”) or stepped-wedge trial (“sw”)
}
  \item{sigma.a}{
Between-cluster SD
}
  \item{sigma.e}{
Individual-level SD (at most, two of rho, sigma.a, and sigma.e should be provided)
}
  \item{sig.level}{
Significance level (Default = 0.05)
}
  \item{method}{
“sim”, “analytic” or “both”. Choice of using simulation, analytic formula or both methods to calculate the power
}
  \item{plot}{
TRUE or FALSE. Whether to print the histogram of pre-randomization power distribution
}
Only used when method = “sim” or “both”
  \item{gen.all}{
Whether generate all possible unique allocations before sampling. If TRUE, all possible unique allocation will be listed.
}
  \item{n.allocs}{
“A” (by default) or a number. The number of randomization order(s) to evaluate. It could be an "A" where all possible unique allocations will be evaluated. If a number is entered, a sample of unique allocations will be evaluated, then the predicted powers for rest allocations will be evaluated
}
  \item{n.sims}{
The number of trials the user wants to simulate to calculation power (should be greater than 1)
}
  \item{seed}{
A number for set.seed () to reproduce simulation
}

Only available when method = “analytic”
  \item{rep}{
The number of times analytical formulae should be repeated when K is non-integer. We recommend a number of around 1,000 to reach a stable estimate. Leave at 1 or NULL when the cluster sizes are integers
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Method = "sim"

\item{coefs}{
A matrix displays the power, estimated coefficient, standard error for estimated coefficient, treatment-vs-time correlation (TTC), control group size (control.size), treatment group size (trt.size), and absolute treatment group imbalance (TGI) for each evaluated allocation.
}
\item{predicted.power}{The predicted power for unevaluated allocations}
\item{PREP}{Average of all attained powers, including both evaluated and predicted powers}
\item{weights}{The weights used for all (evaluated) allocations}
\item{PREP}{Pre-randomization expected power according to CV-based formulae (not for power.strat.pd)}
\item{risk}{The risk of obtaining an attained power lower than the threshold (pwr.thrd)}
\item{allocation}{A matrix of all possible unique allocations}

Method = "analytic"
\item{attained.power}{Attained power for each unique allocation}
\item{PREP}{Average of all attained powers}
\item{allocation}{A matrix of all possible unique allocations or user-specified allocations}
\item{PREP.CV}{Pre-randomization expected power according to CV-based formulae (not for power.strat.pd)}
\item{risk}{The risk of obtaining an attained power lower than the threshold (pwr.thrd)}
}
\references{
Longford NT. Logistic regression with random coefficients. Computational Statistics & Data Analysis 1994; 17: 1–15.

Hussey MA, Hughes JP. Design and analysis of stepped wedge cluster randomized trials. Contemp Clin Trials 2007; 28: 182–191.

}
\author{
Yongdong Ouyang, Liang Xu, Hubert Wong
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
power.strat.pd(I = c(6, 2), P = list("a"= c(1, 1, 1, 1),"b"= c(1, 1, 1, 1)),
               S = factor(c("a", "b", "b", "a", "b", "a", "a", "b")),
               K = c(rep(6, 6),rep(31, 2)), mu0 = 1, Tx.effect = 0.5,
               Time.effect = 0.5*log(0.5), pwr.thrd = NULL, factor.time = FALSE, design = "sw",
               rho = 0.1,n.allocs = "A", plot = FALSE, sig.level = 0.05, family = "poisson",
               method ="sim", seed = 111,
               n.sims = 1000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
