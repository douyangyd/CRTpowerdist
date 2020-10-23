\name{power.pd}
\alias{power.pd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate attained powers and construct power distribution for CRT
}
\description{
Use simulations or analytic formula to calculate attained powers and construct pre-randomization power distribution for all allocations that a randomization algorithm can produce
}
\usage{
power.pd(I, P, K, mu0, Tx.effect, Time.effect = NULL, pwr.thrd = NULL, factor.time = TRUE, family,  design, rho = NULL, sigma.e = NULL, sigma.a = NULL, sig.level = 0.05, method = "analytic", plot = FALSE, gen.all = TRUE, n.allocs = NULL, n.sims, seed = NULL, rep = NULL)
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
  \item{pwr.thrd}{
The threshold used to define “low power”. If NULL no risk will be calculated, risk in the output returns NA by default
}
  \item{factor.time}{
Logical values. Whether treating time as a categorical variable in  the analysis. In simulation, user can treat time differently in data generator and analysis model. (i.e., categorical time effect in data generator, but continuous time effect in the analysis). However, these two must match in analytic calculation. In P-CRT, the value is ignored
}
  \item{family}{
Indicate the type of outcome. One of “gaussian”, “binomial” or “poisson” (case sensitive)
}
  \item{design}{
Can be: parallel cluster randomized trial (“pcrt”) or stepped-wedge trial (“sw”)
}
  \item{rho}{
Intra-cluster coefficient (ICC)
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

Only used when method = “sim” and “both”
  \item{gen.all}{
Whether generate all possible unique allocations before sampling. If TRUE, all possible unique allocation will be listed
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
The number of times analytic formulae should be repeated when K is non-integer. We recommend a number of around 1,000 to reach a stable estimate. Leave at 1 or NULL when the cluster sizes are integers.
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
\item{PREP}{Average of all attained powers using analytic formulae}
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
power.pd(I = c(9, 3), P = c(3, 3, 3, 3), K = c(rep(6, 9), rep(62, 3)),
         mu0 = 0, Tx.effect = 0.26, Time.effect = 1, pwr.thrd = 0.75, factor.time = FALSE,
         design = "sw", rho = 0.01, sigma.e = 1, plot= TRUE, sig.level = 0.05,
         family = "gaussian", method = "sim", seed = 111, gen.all = TRUE,
         n.allocs = "A",  n.sims = 1000)
}
