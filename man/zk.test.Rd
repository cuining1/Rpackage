\name{zk.test}
\alias{zk.test}
\title{
The analogue ZK of the Kolmogorov-Smirnov statistic
}
\description{
The new statistics ZK appear similar to the Kolmogorov-Smirnov statistic, but it's generally much more powerful,see Jin Zhang(2002).
}
\usage{
zk.test(x, y, para = NULL, N = 1000)
}
\arguments{
  \item{x}{
a numeric vector of data values, the number of which must be greater than 7. Missing values are allowed.
}
  \item{y}{
When tested as a single sample, a numeric vector of data values, the number of which must be greater than 7. Missing values are allowed.

When tested as two-sample, a character string indicating what type of test was performed.
Distributions "unif", "exponential", "normal", "lognormal", "gamma","t" and "weibull" are recognised.Here,the names of "exponential", "normal" and "lognormal" are simplified as "exp" ,"norm" and "lognorm" respectively.
}
  \item{para}{
A named list giving the parameters of the distribution specified and this can be omitted.
}
  \item{N}{
The number of replications in Monte Carlo simulation.The default value is 1000
}
}
\details{
The Zk test is an EDF omnibus test for the composite hypothesis of distribution.
The test statistic is
\deqn{Z_{k}=\max \limits_{1 \leq i \leq n}\{(i-0.5)\ln\frac{i-0.5}{nF_{0}(X_{(i)})}+(n-i+0.5)\ln\frac{n-i+0.5}{n[1-F_0(X_{(i)})]}\},
}
where \eqn{F_{0}(x)} is a hypothesized distribution function to be tested.Here, \eqn{F_{0}(X_{(i)}) = \Phi(x)}, \eqn{\Phi} is the cumulative distribution function of the specificed distribution.The p-value is computed  by Monte Carlo simulation.
}
\value{
A list with class \dQuote{htest} containing the following components:
  \item{statistic}{the value of the ZK statistic.}
  \item{p.value }{the p-value for the test.}
  \item{method}{the character string \dQuote{ZK test for given distribution}.}
  \item{data.name}{a character string giving the name(s) of the data.}
}
\references{
Jin Zhang: Goodness-of-Fit Tests Based on the Likelihood Ratio.Journal of the Royal Statistical Society,64,281-294.

Jin Zhang,Yuehua Wu: Likelihood-ratio tests for normality.Computational Statistics & Data Analysis,49,709-721.

Jin Zhang: Powerful Two-Sample Tests Based on the Likelihood Ratio. Technometrics, 48:1, 95-103.
}
\author{
Ning Cui
}
\note{
The Zk test is the recommended EDF test by Jin Zhang.
}

\seealso{
\code{\link{ks.test}} for performing a one- or two-sample Kolmogorov-Smirnov test.
\code{\link{zc.test}},\code{\link{za.test}} for performing a powerful goodness-of-fit test based on the likelihood ratio.
}

\examples{
x<-rbeta(50,shape1 = 0.6,shape2 = 0.8)
y<-rnorm(50)
zk.test(x,y)
zk.test(x,"unif")
zk.test(x,"norm")
zk.test(x,"unif",para = list(min=1,max=2))
zk.test(x,"exp",para = list(rate=1))
zk.test(x,"norm",para = list(mean=1,sd=2))
zk.test(x,"lognorm",para = list(mean=1,sd=2))
zk.test(x,"weibull",para = list(shape=1,scale=2))
zk.test(x,"gamma",para = list(shape=2,scale=1))
zk.test(x,"t",para = list(df=3))
}
\keyword{ htest }
