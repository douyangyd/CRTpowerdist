\name{power.ap}
\alias{power.ap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate attained power for specific allocation(s)
}
\description{
Use simulations or analytic formula to calculate attained power for specific allocations
}
\usage{
power.ap (I, P , K , mu0, Tx.effect, Time.effect = NULL, user.allocs,
          factor.time = TRUE, family, design, rho = NULL,
          sigma.e = NULL, sigma.a = NULL, sig.level = 0.05,
          method = "analytic", adj.pwr=FALSE, seed = NULL,
          n.sims = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{I}{
A vector indicating the number of clusters of each group (i.e., 8 clusters with six small and two large cluster has I = c(6, 2))
}
  \item{P}{
For SW-CRT, P is a vector indicating the number of clusters transitioning at each step. Length of P should be the total number of the transition steps (e.g., if four transition steps, two clusters transitioning at each step, then P = c(2, 2, 2, 2))
For P-CRT, P is a vector indicating the number of clusters allocated to the control and intervention groups (e.g., P = c(5, 3) means 5 clusters in control and 3 clusters in intervention groups)
}
  \item{K}{
K is a vector indicating the number of participants in each cluster.
Length of K equals the total number of clusters. K should be  ordered to match I (e.g., I = c(6, 2) and K = c(6, 6, 6, 6, 6, 6, 12, 12) means the first six small clusters have six pts per cluster and the two large clusters have 12 pts per cluster)
}
  \item{mu0}{
Baseline mean/risk/rate on the outcome’s natural scale (i.e., not transformed by the link function)
}
  \item{Tx.effect}{
Treatment effect (natural scale: risk difference Gaussian outcome, odds ratio for binary outcome, risk ratio for Poisson outcome)
}
  \item{Time.effect}{
Time effect (linear predictor scale). Default is no time effect. The time effect can be a continuous linear trend (a single number) or categorical time effect (a vector of J numbers reflecting the time effects for each period excluding baseline)
}
  \item{user.allocs}{
A matrix. Each row represents a specific allocation (a vector indicating when the cluster transit from control to intervention) that the user wants to evaluate for power
}
  \item{factor.time}{
Logical indicator for whether time is treated as a categorical variable in the analysis. In simulation, user can treat time differently in data generator and analysis model. (i.e., categorical time effect in data generator, but continuous time effect in the analysis). However, these two must match in analytic calculation. For the P-CRT, this value is ignored
}
  \item{family}{
Indicate the type of outcome. One of "gaussian", "binomial" or "poisson" (case sensitive)
}
  \item{design}{
Can be "pcrt" for a P-CRT or "sw" for a SW-CRT
}
  \item{sig.level}{
Significance level (Default=0.05). If you have different allocations to evaluate the sig.level can be a vector in which the i-th entry specifies the significance level for the allocation in the i-th row of user.allocs
}
  \item{rho}{
Intra-cluster coefficient (ICC)
}
  \item{sigma.e}{
Between-cluster SD
}
  \item{sigma.a}{
Individual-level SD (at most, two of rho, sigma.a, and sigma.e should be provided)
}
  \item{method}{
"sim", "analytic" or "both". Choice of using simulation, analytic formula or both methods to calculate the power
}
  \item{seed}{
  A number for set.seed () to reproduce simulation
}
Only used with method = "sim" or "both"
  \item{adj.pwr}{
Whether evaluate power adjusted for the P-value cut-off that corresponds to achieving the nominal Type I error rate
}
  \item{n.sims}{
The number of trials the user wants to simulate to calculation power (should be greater than 1)
}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}

\value{
A list with two named components "results" and "inputs", each of which is a list. Components of "results" report the results of the function call while the components of "inputs" reproduce the inputs to the function call. Components associated with a specific method will appear only if that method was requested in the function call (e.g., the component "attained.power.analytic" is generated only if method = "analytic" or "both").

\bold{$result}
\item{attained.power.sim}{Attained power for the allocation(s) evaluated using simulation}
\item{trt.estimate}{The average of the treatment effect coefficients across simulation runs}
\item{SE}{The average of the standard errors for the estimated treatment coefficients across simulation runs}
\item{TGI}{Absolute treatment group imbalance for evaluated allocation(s)}
\item{TTC}{Treatment-vs-time correlation for evaluated allocation(s)}
\item{attained.power.analytic}{Attained power for the allocation(s) evaluated using analytic approximation}
\item{CV}{Coefficient of variation of the cluster sizes used in the analytic approximation}


\bold{$inputs}
\item{I, P, K, ..., etc}{Input values to the function call}
}

\references{
Longford NT. Logistic regression with random coefficients. Computational Statistics & Data Analysis 1994; 17: 1-15.

Hussey MA, Hughes JP. Design and analysis of stepped wedge cluster randomized trials. Contemp Clin Trials 2007; 28: 182-191.

Karla Hemming, Jessica Kasza, Richard Hooper, Andrew Forbes, Monica Taljaard, A tutorial on sample size calculation for multiple-period cluster randomized parallel, cross-over and stepped-wedge trials using the Shiny CRT Calculator, International Journal of Epidemiology, Volume 49, Issue 3, June 2020, Pages 979-995.

Ouyang Y, Karim ME, Gustafson P, Field TS, Wong H. Explaining the variation in the attained power of a stepped-wedge trial with unequal cluster sizes. BMC Med Res Methodol. 2020;20(1):166. Published 2020 Jun 24. doi:10.1186/s12874-020-01036-5

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
## Continous time effect
power.ap(I = c(9, 3), P = c(3, 3, 3, 3), K = c(rep(6, 9), rep(62, 3)),
         user.allocs = rbind(c(1, 1, 1, 2, 2, 2,
                               3, 3, 3, 4, 4, 4)),
         mu0 = 0, Tx.effect = 0.26, Time.effect = 1,
         factor.time = FALSE, design = "sw", rho = 0.01, sigma.e = 1,
         sig.level = 0.05, adj.pwr =FALSE, family = "gaussian",
         method = "sim", seed = 222,
         n.sims = 2000)

## Categorical time effect
power.ap(I = c(9, 3), P = c(3, 3, 3, 3), K = c(rep(30, 9), rep(310, 3)),
         user.allocs = rbind(c(1, 1, 1, 2, 2, 2,
                               3, 3, 3, 4, 4, 4)),
         mu0 = 0, Tx.effect = 0.26, Time.effect = c(1,2,3,4), factor.time = FALSE,
         design = "sw", rho = 0.01, sigma.e = 1, sig.level = 0.05, adj.pwr =FALSE,
         family = "gaussian", method = "sim", seed = 222,
         n.sims = 2000)
}
