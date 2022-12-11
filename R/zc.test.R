zc.test <- function(x, y, para = NULL, N = 1000) {
  DNAME <- deparse(substitute(x))
  x <- x[!is.na(x)]
  n1 <- length(x)
  if (n1 < 8)
    stop("sample size must be greater than 7")
  if (is.numeric(y)) {
    DNAME <- paste(DNAME, "and", deparse(substitute(y)))
    METHOD <- "Two-sample Zc test"
    y <- y[!is.na(y)]
    n2 <- length(y)
    R <- rank(c(x, y))
    n <- n1 + n2
    S <- 0
    g <- function(m, r, M) sum(log(m/(1:m - 0.5) - 1) * log(M/(r -
                                                                 0.5) - 1))
    ZC <- (g(n1, sort(R[1:n1]), n) + g(n2, sort(R[(n1 + 1):n]), n))/n
    for (j in 1:N) {
      R <- sample(n)
      zc <- (g(n1, sort(R[1:n1]), n) + g(n2, sort(R[(n1 + 1):n]),
                                         n))/n
      S <- S + (zc < ZC)
    }
    p.value <- S/N
  }
  if (is.character(y)) {
    METHOD <- "One-sample Zc test"
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

    q <- (n1 - 0.5)/(seq(1:n1) - 0.75) - 1
    h <- log((p^(-1) - 1)/q)
    ZC <- sum(h^2)

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
        q <- (n1 - 0.5)/(seq(1:n1) - 0.75) - 1
        h <- log((p^(-1) - 1)/q)
        zc <- sum(h^2)
        S <- S + (zc > ZC)
      }
    }
    if (distr == "exp") {
      for (j in 1:N) {
        R <- rexp(n1, a)
        R <- sort(R)
        p <- pexp(R, rate = a, lower.tail = TRUE, log.p = FALSE)
        q <- (n1 - 0.5)/(seq(1:n1) - 0.75) - 1
        h <- log((p^(-1) - 1)/q)
        zc <- sum(h^2)
        S <- S + (zc > ZC)
      }
    }
    if (distr == "norm" || distr == "lognorm") {
      if (is.null(para)) {
        for (j in 1:N) {
          R <- rnorm(n1)
          R <- sort(R)
          p <- pnorm((R - mean(R))/sd(R), mean = 0, sd = 1, lower.tail = TRUE,
                     log.p = FALSE)
          q <- (n1 - 0.5)/(seq(1:n1) - 0.75) - 1
          h <- log((p^(-1) - 1)/q)
          zc <- sum(h^2)
          S <- S + (zc > ZC)
        }
      } else {
        for (j in 1:N) {
          R <- rnorm(n1, a, b)
          R <- sort(R)
          p <- pnorm(R, mean = a, sd = b, lower.tail = TRUE, log.p = FALSE)
          q <- (n1 - 0.5)/(seq(1:n1) - 0.75) - 1
          h <- log((p^(-1) - 1)/q)
          zc <- sum(h^2)
          S <- S + (zc > ZC)
        }
      }
    }
    if (distr == "gamma") {
      for (j in 1:N) {
        R <- rgamma(n1, shape = a, scale = b)
        R <- sort(R)
        p <- pgamma(R, shape = a, scale = b, lower.tail = TRUE,
                    log.p = FALSE)
        q <- (n1 - 0.5)/(seq(1:n1) - 0.75) - 1
        h <- log((p^(-1) - 1)/q)
        zc <- sum(h^2)
        S <- S + (zc > ZC)
      }
    }
    if (distr == "weibull") {
      for (j in 1:N) {
        R <- rweibull(n1, shape = a, scale = b)
        R <- sort(R)
        p <- pweibull(R, shape = a, scale = b, lower.tail = TRUE,
                      log.p = FALSE)
        q <- (n1 - 0.5)/(seq(1:n1) - 0.75) - 1
        h <- log((p^(-1) - 1)/q)
        zc <- sum(h^2)
        S <- S + (zc > ZC)
      }
    }
    if (distr == "t") {
      for (j in 1:N) {
        R <- rt(n1, df = a)
        R <- sort(R)
        p <- pt(R, df = a, lower.tail = TRUE, log.p = FALSE)
        q <- (n1 - 0.5)/(seq(1:n1) - 0.75) - 1
        h <- log((p^(-1) - 1)/q)
        zc <- sum(h^2)
        S <- S + (zc > ZC)
      }
    }
    p.value <- S/N
  }
  names(ZC) <- "Zc"
  RVAL <- list(statistic = ZC, p.value = p.value, alternative = NULL,
               method = METHOD, data.name = DNAME)
  class(RVAL) <- "htest"
  return(RVAL)
}

