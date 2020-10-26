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
power.ap(I, P, K, mu0, Tx.effect, Time.effect = NULL, user.allocs, factor.time = TRUE, family, design, rho = NULL, sigma.e = NULL, sigma.a = NULL, sig.level = 0.05, method = "analytic", adj.pwr = FALSE, seed = NULL, n.sims, rep = NULL)
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
For SW,  K is the cluster size per period for each cluster. Then the total cluster size is K× (number of transition steps+1).
For PCRT, K is the cluster size.
Length of K equals sum(I) - total number of clusters. K should be entered as the order in I (i.e., I = c(6, 2) and K = c(6, 6, 6, 6, 6, 6, 12, 12) means the first six small clusters have 6 pts per period and the two large clusters have 12 pts per period).
K can be a non-integer. If the resulting cluster size is also a non-integer, the cluster size is determined by randomly picking one of the two nearest integers with probabilities that ensure the non-integer cluster size is achieved in expectation (e.g., a cluster has K = 6.1 with four transition periods, the total cluster size is 30.5. Therefore, cluster size for this cluster will be either 30 or 31 with probabilities 0.9 an 0.1)
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
  \item{user.allocs}{
A matrix. Each row represents a specific allocation (a vector indicating when the cluster transit from control to intervention) that the user wants to evaluate for power
}
  \item{factor.time}{
Logical values. Whether treating time as a categorical variable in  the analysis. In simulation, user can treat time differently in data generator and analysis model. (i.e., categorical time effect in data generator, but continuous time effect in the analysis). However, these two must match in analytic calculation. For the P-CRT, this value is ignored.
}
  \item{family}{
Indicate the type of outcome. One of "gaussian", "binomial" or "poisson" (case sensitive)
}
  \item{design}{
Can be: parallel cluster randomized trial ("pcrt") or stepped-wedge trial ("sw")
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
  \item{sig.level}{
Significance level (Default=0.05). If you have different allocations to evaluate the sig.level can be a vector in which the i-th entry specifies the significance level for the allocation in the i-th row of user.allocs
}
  \item{method}{
"sim", "analytic" or "both". Choice of using simulation, analytic formula or both methods to calculate the power
}
  \item{adj.pwr}{
Whether evaluate power adjusted for the P-value cut-off that corresponds to achieving the nominal Type I error rate
}
  \item{seed}{
  A number for set.seed () to reproduce simulation
}

Only used with method = "sim" or "both"
  \item{n.sims}{
The number of trials the user wants to simulate to calculation power (should be greater than 1)
}

Only used with method = "analytic"
  \item{rep}{
The number of times analytical formulae should be repeated when K is non-integer. We recommend a number of around 1,000 to reach a stable estimate. Leave at 1 or NULL when the cluster sizes are integers.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Method = "sim"
A matrix with columns representing the power, control group size, treatment group size, treatment effect estimate, and standard error for the evaluated allocation(s)


Method = "analytic"
\item{attained.power}{Attained power for each unique allocation}
\item{CV}{Coefficient of variation of cluster sizes}
\item{allocation}{A matrix of all possible unique allocations or user-specified allocations}
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
power.ap(I = c(18, 6), P = c(6, 6, 6, 6), K = c(rep(7, 18), rep(15, 6)),
         user.allocs = rbind(c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2,
                               3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4)),
         factor.time = FALSE, mu0 = 0, Tx.effect = 0.26, Time.effect = 1,
         family = "gaussian", design = "sw", sig.level = 0.05, rho = 0.01,
         sigma.a = NULL, sigma.e = 1, adj.pwr =FALSE, method="sim",
         seed = 111, n.sims = 2000)
}
