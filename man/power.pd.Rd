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
power.pd (I, P, K, mu0, Tx.effect, Time.effect = NULL,
          pwr.thrd = NULL, factor.time = TRUE, design, rho = NULL,
          family, sig.level = 0.05, sigma.e = NULL, sigma.a = NULL,
          method = "analytic", plot = FALSE, gen.all = TRUE, n.allocs,
          n.sims, seed = NULL)
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
Length of K equals the total number of clusters. K should be  ordered to match I (e.g., I = c(6, 2) and K = c(6, 6, 6, 6, 6, 6, 12, 12) means the first six small clusters have 6 pts per period and the two large clusters have 12 pts per period)
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
  \item{pwr.thrd}{
The threshold used to define "low power". If NULL no risk will be calculated, risk in the output returns NA by default
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
"sim", "analytic" or "both". Choice of using simulation, analytic formula or both methods to calculate the power
}
  \item{plot}{
TRUE or FALSE. Whether to print the histogram of pre-randomization power distribution
}

Only used when method = "sim" and "both"
  \item{gen.all}{
Whether generate all possible unique allocations before sampling. If TRUE, all possible unique allocation will be listed
}
  \item{n.allocs}{
NULL (by default) or a number. The number of unique allocations to evaluate via simulation. If NULL, all possible unique allocations will be evaluated. If a number, a random sample of n.allocs unique allocations will be evaluated, and the powers for the remaining allocations will be predicted
}
  \item{n.sims}{
The number of trials to simulate when calculating the attained power (should be greater than 1)
}
  \item{seed}{
A number for set.seed () to enable reproducibility of simulation results
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list with two named components "results" and "inputs", each of which is a list. Components of "results" report the results of the function call while the components of "inputs" reproduce the inputs to the function call. Components associated with a specific method will appear only if that method was requested in the function call (e.g., the component "attained.power.analytic" is generated only if method = "analytic" or "both"). If plot = TRUE, a plot of the PD will be displayed.

\bold{$results}
\item{allocations}{A matrix with each row identifying one allocation used in the evaluation of the PD}
\item{weights}{The weights (proportional to the multiplicities of each allocation) applied to the attained powers when calculating the PD}
\item{attained.power.sim}{A vector of the attained power for each allocation evaluated using simulation (or prediction)}
\item{coefs}{A matrix displaying the average treatment effect coefficient, the average standard error for coefficient, the treatment-vs-time correlation (TTC), the control group size (control.size), the treatment group size (trt.size), and the absolute treatment group imbalance (TGI) for each allocation evaluated by simulation.}
\item{PREP.sim}{PREP calculated as a weighted average of all attained powers evaluated by simulation or prediction}
\item{risk.sim}{The risk of obtaining an attained power lower than the threshold (pwr.thrd) based on the PD obtained via simulation}
\item{attained.power.analytic}{A vector of the attained power for each allocation evaluated using analytic approximation}
\item{PREP.analytic}{PREP calculated as a weighted average of all attained powers using analytic approximation}
\item{CV}{Coefficient of variation of the cluster sizes used in the analytic approximation}
\item{PREP.CV}{PREP calculated according to CV-based formulae}
\item{risk.analytic}{The risk of obtaining an attained power lower than the threshold (pwr.thrd) based on the PD obtained via analytic approximation}

\bold{$inputs}
\item{I, P, K, ..., etc}{Input values to the function call}
}

\references{
Longford NT. Logistic regression with random coefficients. Computational Statistics & Data Analysis 1994; 17: 1-15.

Hussey MA, Hughes JP. Design and analysis of stepped wedge cluster randomized trials. Contemp Clin Trials 2007; 28: 182-191.

Karla Hemming, Jessica Kasza, Richard Hooper, Andrew Forbes, Monica Taljaard, A tutorial on sample size calculation for multiple-period cluster randomized parallel, cross-over and stepped-wedge trials using the Shiny CRT Calculator, International Journal of Epidemiology, Volume 49, Issue 3, June 2020, Pages 979-995.

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
