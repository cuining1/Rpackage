zk.test <- function(x, y, para = NULL, N = 1000) {
  DNAME <- deparse(substitute(x))
  x <- x[!is.na(x)]
  n1 <- length(x)
  if (n1 < 8)
    stop("sample size must be greater than 7")
  if (is.numeric(y)) {
    DNAME <- paste(DNAME, "and", deparse(substitute(y)))
    METHOD <- "Two-sample Zk test"
    y <- y[!is.na(y)]
    n2 <- length(y)
    R <- ceiling(rank(c(x, y)))
    n <- n1 + n2
    P <- (1:n - 0.5)/n
    w <- n * (P * log(P) + (1 - P) * log(1 - P))
    g <- function(m, r, M) {
      d <- sort(r)
      D <- c(1, d, M + 1)
      p <- rep(0:m, D[2:(m + 2)] - D[1:(m + 1)])
      p[d] <- p[d] - 0.5
      p <- p/m
      m * (p * log(p + 1e-10) + (1 - p) * log(1 - p + 1e-10))
    }
    ZK <- max(g(n1, R[1:n1], n) + g(n2, R[(n1 + 1):n], n) - w)
    S <- 0
    for (j in 1:N) {
      R <- sample(n)
      zk <- max(g(n1, R[1:n1], n) + g(n2, R[(n1 + 1):n], n) - w)
      S <- S + (zk > ZK)
    }
    p.value <- S/N
  }
  if (is.character(y)) {
    METHOD <- "One-sample Zk test"
    x <- sort(x)
    distr <- y
    Distr <- c("unif", "exp", "norm", "lognorm", "gamma", "weibull",
               "t")
    if (!is.element(distr, Distr))
      stop("unsupported distribution")
    if (missing(para))
      para <- NULL
    if (distr == "unif") {
      if (is.null(para)) {
        # Moment estimation
        a <- mean(x) - sqrt(3) * sd(x)
        b <- mean(x) + sqrt(3) * sd(x)
      } else {
        a <- para$min
        b <- para$max
      }
      p <- punif(x, min = a, max = b, lower.tail = TRUE, log.p = FALSE)
      for (i in 1:n1) {
        if (p[i] == 0) {
          p[i] = p[i] + 1e-10
        }
        if (p[i] == 1) {
          p[i] = p[i] - 1e-10
        }
      }
    }
    if (distr == "exp") {
      if (any(x <= 0)) {
        stop("sample must be greater than 0")
      }
      if (is.null(para)) {
        a <- 1/mean(x)  # Maximum likelihood estimation
      } else {
        a <- para$rate
      }
      p <- pexp(x, rate = a, lower.tail = TRUE, log.p = FALSE)
      for (i in 1:n1) {
        if (p[i] == 0) {
          p[i] = p[i] + 1e-10
        }
        if (p[i] == 1) {
          p[i] = p[i] - 1e-10
        }
      }
    }
    if (distr == "norm") {
      if (is.null(para)) {
        p <- pnorm((x - mean(x))/sd(x), mean = 0, sd = 1, lower.tail = TRUE,
                   log.p = FALSE)
      } else {
        a <- para$mean
        b <- para$sd
        p <- pnorm(x, mean = a, sd = b, lower.tail = TRUE, log.p = FALSE)
      }
    }
    if (distr == "lognorm") {
      if (any(x <= 0)) {
        stop("sample must be greater than 0")
      }
      if (is.null(para)) {
        p <- pnorm((log(x) - mean(log(x)))/sd(log(x)), mean = 0,
                   sd = 1, lower.tail = TRUE, log.p = FALSE)
      } else {
        a <- para$mean
        b <- para$sd
        p <- pnorm(log(x), mean = a, sd = b, lower.tail = TRUE,
                   log.p = FALSE)
      }
    }
    if (distr == "gamma") {
      if (any(x <= 0)) {
        stop("sample must be greater than 0")
      }
      if (is.null(para)) {
        # Maximum likelihood estimation
        H <- log(mean(x)/((cumprod(x)[n1])^(1/n1)))
        a <- 0
        if (H > 0 & H <= 0.5772) {
          a <- (0.500876 + 0.1648852 * H - 0.054427 * H^2)/H
        }
        if (H > 0.6772 & H <= 17) {
          a <- (8.898919 + 9.05995 * H + 0.9775373 * H^2)/(H *
                                                             (17.79728 + 11.968477 * H + H^2))
        }
        b <- mean(x)/a
      } else {
        a <- para$shape
        b <- para$scale
      }
      p <- pgamma(x, shape = a, scale = b, lower.tail = TRUE, log.p = FALSE)
    }
    if (distr == "weibull") {
      if (any(x <= 0)) {
        stop("sample must be greater than 0")
      }
      if (is.null(para)) {
        # Estimation based on extreme value distribution
        a <- 1.283/sd(log(x))
        b <- exp(mean(log(x)) + 0.5772 * 0.7797 * sd(log(x)))
      } else {
        a <- para$shape
        b <- para$scale
      }
      p <- pweibull(x, shape = a, scale = b, lower.tail = TRUE, log.p = FALSE)
    }
    if (distr == "t") {
      if (is.null(para)) {
        s <- fitdistr(x, densfun = "t")
        a <- unname(s$estimate)[3]
      } else {
        a <- para$df
      }
      p <- pt(x, df = a, lower.tail = TRUE, log.p = FALSE)
    }

    logp1 <- log((seq(1:n1) - 0.5)/(n1 * p))
    logp2 <- log((n1 - seq(1:n1) + 0.5)/(n1 * (1 - p)))
    h <- (seq(1:n1) - 0.5) * logp1 + (n1 - seq(1:n1) + 0.5) * logp2
    ZK <- max(h)

    # Monte Carlo simulation
    S <- 0
    if (distr == "unif") {
      for (j in 1:N) {
        R <- runif(n1, a, b)
        R <- sort(R)
        p <- punif(R, min = a, max = b, lower.tail = TRUE, log.p = FALSE)
        for (i in 1:n1) {
          if (p[i] == 0) {
            p[i] = p[i] + 1e-10
          }
          if (p[i] == 1) {
            p[i] = p[i] - 1e-10
          }
        }
        logp1 <- log((seq(1:n1) - 0.5)/(n1 * p))
        logp2 <- log((n1 - seq(1:n1) + 0.5)/(n1 * (1 - p)))
        h <- (seq(1:n1) - 0.5) * logp1 + (n1 - seq(1:n1) + 0.5) *
          logp2
        zk <- max(h)
        S <- S + (zk > ZK)
      }
    }
    if (distr == "exp") {
      for (j in 1:N) {
        R <- rexp(n1, a)
        R <- sort(R)
        p <- pexp(R, rate = a, lower.tail = TRUE, log.p = FALSE)
        logp1 <- log((seq(1:n1) - 0.5)/(n1 * p))
        logp2 <- log((n1 - seq(1:n1) + 0.5)/(n1 * (1 - p)))
        h <- (seq(1:n1) - 0.5) * logp1 + (n1 - seq(1:n1) + 0.5) *
          logp2
        zk <- max(h)
        S <- S + (zk > ZK)
      }
    }
    if (distr == "norm" || distr == "lognorm") {
      if (is.null(para)) {
        for (j in 1:N) {
          R <- rnorm(n1)
          R <- sort(R)
          p <- pnorm((R - mean(R))/sd(R), mean = 0, sd = 1, lower.tail = TRUE,
                     log.p = FALSE)
          logp1 <- log((seq(1:n1) - 0.5)/(n1 * p))
          logp2 <- log((n1 - seq(1:n1) + 0.5)/(n1 * (1 - p)))
          h <- (seq(1:n1) - 0.5) * logp1 + (n1 - seq(1:n1) + 0.5) *
            logp2
          zk <- max(h)
          S <- S + (zk > ZK)
        }
      } else {
        for (j in 1:N) {
          R <- rnorm(n1, a, b)
          R <- sort(R)
          p <- pnorm((R - mean(R))/sd(R), mean = 0, sd = 1, lower.tail = TRUE,
                     log.p = FALSE)
          logp1 <- log((seq(1:n1) - 0.5)/(n1 * p))
          logp2 <- log((n1 - seq(1:n1) + 0.5)/(n1 * (1 - p)))
          h <- (seq(1:n1) - 0.5) * logp1 + (n1 - seq(1:n1) + 0.5) *
            logp2
          zk <- max(h)
          S <- S + (zk > ZK)
        }
      }
    }
    if (distr == "gamma") {
      for (j in 1:N) {
        R <- rgamma(n1, shape = a, scale = b)
        R <- sort(R)
        p <- pgamma(R, shape = a, scale = b, lower.tail = TRUE,
                    log.p = FALSE)
        logp1 <- log((seq(1:n1) - 0.5)/(n1 * p))
        logp2 <- log((n1 - seq(1:n1) + 0.5)/(n1 * (1 - p)))
        h <- (seq(1:n1) - 0.5) * logp1 + (n1 - seq(1:n1) + 0.5) *
          logp2
        zk <- max(h)
        S <- S + (zk > ZK)
      }
    }
    if (distr == "weibull") {
      for (j in 1:N) {
        R <- rweibull(n1, shape = a, scale = b)
        R <- sort(R)
        p <- pweibull(R, shape = a, scale = b, lower.tail = TRUE,
                      log.p = FALSE)
        logp1 <- log((seq(1:n1) - 0.5)/(n1 * p))
        logp2 <- log((n1 - seq(1:n1) + 0.5)/(n1 * (1 - p)))
        h <- (seq(1:n1) - 0.5) * logp1 + (n1 - seq(1:n1) + 0.5) *
          logp2
        zk <- max(h)
        S <- S + (zk > ZK)
      }
    }
    if (distr == "t") {
      for (j in 1:N) {
        R <- rt(n1, df = a)
        R <- sort(R)
        p <- pt(R, df = a, lower.tail = TRUE, log.p = FALSE)
        logp1 <- log((seq(1:n1) - 0.5)/(n1 * p))
        logp2 <- log((n1 - seq(1:n1) + 0.5)/(n1 * (1 - p)))
        h <- (seq(1:n1) - 0.5) * logp1 + (n1 - seq(1:n1) + 0.5) *
          logp2
        zk <- max(h)
        S <- S + (zk > ZK)
      }
    }
    p.value <- S/N
  }
  names(ZK) <- "Zk"
  RVAL <- list(statistic = ZK, p.value = p.value, alternative = NULL,
               method = METHOD, data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
}
