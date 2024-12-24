bVAR <- R6::R6Class(
  "Bayes VAR",
  inherit = VAR,
  public = list(
    # prior and post parameters
    prior.type = NULL,
    alpha.prior = NULL,
    Sigma.prior = NULL,
    alpha.post = NULL,
    Sigma.post = NULL,
    B.NSR = NULL,
    Sigma.NSR = NULL,
    beta.NSR = NULL,
    get.result.tmp = FALSE,
    initialize = function(data = NA, p.lag = NA, prior.type = "info", priors = NA, cons = 1) {
      # need to check data and p.lag input
      super$initialize(data, p.lag)
      self$prior.type <- prior.type
      if (!is.na(priors)) {
        self$alpha.prior <- list(
          "mean" = priors[[1]],
          "cov" = priors[[2]]
        )
        self$Sigma.prior <- list(
          "mean" = priors[[3]],
          "nu" = priors[[4]]
        )
        self$prior.type <- "user"
        cat("* Bayesian VAR priors set by user.\n")
      } else if (prior.type == "flat") {
        self$alpha.prior <- list(
          "mean" = as.vector(t(cbind(rep(0, self$n.var), diag(n.var), matrix(0, self$n.var, self$n.var * (self$p.lag - 1))))),
          "cov" = 0.1 * diag(self$n.var * (self$n.var * self$p.lag + 1))
        )
        self$Sigma.prior <- list(
          "mean" = diag(self$n.var),
          "nu" = 0
        )
        cat("* Bayesian VAR priors set as cheap guess.\n")
      } else if (prior.type == "info") {
        if (det(t(self$X) %*% self$X) == 0) stop("X'X not invertible, cannot use OLS estimates as prior!")
        self$beta <- solve(t(self$X) %*% self$X, t(self$X) %*% self$Y)
        self$U <- self$Y - self$X %*% self$beta
        self$Sigma <- (t(self$U) %*% self$U) / (self$T.est)
        self$alpha.prior <- list(
          "mean" = as.vector(self$beta),
          "cov" = 0.1 * diag(self$n.var * (self$n.var * self$p.lag + 1))
        )
        self$Sigma.prior <- list(
          "mean" = self$Sigma,
          "nu" = 0
        )
        cat("* Bayesian VAR priors set as OLS estimates.\n")
      } else {
        stop("please follow the instructions to set priors")
      }
    },
    est = function(Y = self$Y, X = self$X, method = "gibbs", burn.in = 100, draw = 1000, thin = 1, post.save = NA, post.method = 0) {
      if (method == "gibbs") {
        if (is.na(post.save)) post.save <- draw / 2
        priors <- c(self$alpha.prior, self$Sigma.prior)
        cat("* Gibbs sampler start...\n")
        tic("Gibbs sampler time usage")
        out <- gibbs_sampler(Y, X, priors, burn.in, draw, thin, post.save, post.method)
        toc()
        ## save posterior parameters
        self$alpha.post <- list("mean" = out$alpha_mean_post, "cov" = out$alpha_cov_post)
        self$Sigma.post <- list("mean" = out$Sigma_mean_post, "nu" = out$nu_post)
        ## temporally save parameters estimated, if not using sign restrictions
        self$beta <- matrix(rowMeans(out$alpha), self$n.var * self$p.lag + 1, self$n.var)
        self$Sigma <- apply(out$Sigma, c(1, 2), mean)
        cat("* posterior parameters estimated using ", method, " method.\n")
        ## return loglik for potential analysis
        return(out$loglik)
      } else if (method == "conjugate") {
        # define quantities following Kilian (2017)
        # V <- 0.1 * diag(self$n.var * self$p.lag + 1)
        # A.star <- t(matrix(self$alpha.prior$mean, self$n.var * self$p.lag + 1, self$n.var))
        # Sigma.mu <- self$Sigma
        # S.star <- self$Sigma.prior$mean
        # conjugate estimation of posteriors, instead of Gibbs sampler, to save time
        conj.out <- private$conjugate.post()
        self$alpha.post <- list("mean" = conj.out$A.bar, "cov" = conj.out$Sigma.alpha.bar)
        self$Sigma.post <- list("mean" = conj.out$S, "nu" = conj.out$tau)
        cat("* posterior parameters estimated using ", method, " method.\n")
      } else {
        stop("please set appropriate estimation method!")
      }
    },
    identify.seq = function(SR = NA, EBR = NA, NSR = NA, draw = 5000, M = 1000, IV = NA) {
      restrictions.out <- private$get.restrictions(SR, EBR, NSR)
      restrictions <- restrictions.out$restrictions
      max.irf <- restrictions.out$max.irf
      is.nsr <- restrictions.out$is.nsr
      res.num <- length(restrictions)
      if (!self$get.result.tmp) {
        cat("* ", res.num, " restrictions get.\n* First impose SR and EBR on drawn samples:\n")
        mystr <- paste("*", as.character(draw), "draws saved with sign and elasticity restrictions satisfied, total time usage")
        tic(msg = mystr)
        sign.out <- impose_SR_and_EBR(self$alpha.post, self$Sigma.post, restrictions, self$Y, self$X, draw, self$p.lag, max.irf)
        toc()
        # save drawn structural parameters and IRFs
        self$beta.set <- sign.out$beta_saved
        self$B.set <- sign.out$B_saved
        self$Sigma.set <- sign.out$Sigma_saved
        self$get.result.tmp <- TRUE
      }
      if (is.nsr) {
        cat("* checking NSR...\n")
        mystr <- "* Narrative restrictions imposed, time usage"
        tic(msg = mystr)
        NSR.out <- impose_NSR(restrictions, self$Y, self$X, self$B.set, self$beta.set, self$Sigma.set, self$p.lag, max.irf, M)
        toc()
        self$beta.NSR <- NSR.out$beta_NSR
        self$B.NSR <- NSR.out$B_NSR
        self$Sigma.NSR <- NSR.out$Sigma_NSR
      }
      cat("* identification done.\n")
    },
    # identify.all = function(SR = NA, EBR = NA, NSR = NA, draw = 5000, M = 1000) {
    #   restrictions.out <- private$get.restrictions(SR, EBR, NSR)
    #   restrictions <- restrictions.out$restrictions
    #   max.irf <- restrictions.out$max.irf
    #   res.num <- length(restrictions)
    #   cat("* ", res.num, " restrictions get.\n* impose all restrictions on drawn samples:\n")
    #   mystr <- paste("*", as.character(draw), "draws saved with all restrictions satisfied, total time usage")
    #   tic(msg = mystr)
    #   sign.out <- impose_all_restrictions(self$alpha.post, self$Sigma.post, restrictions, self$Y, self$X, draw, self$p.lag, max.irf, M)
    #   toc()
    #   # save drawn structural parameters and IRFs
    #   self$beta.draw <- sign.out$beta_saved
    #   self$B.draw <- sign.out$B_saved
    #   self$Sigma.draw <- sign.out$Sigma_saved
    #   cat("identification done.\n")
    # },
    tool = function(tool = NA_character_, hor = NA_integer_, prob = NA_real_, start = NA_integer_, end = NA_integer_, which.var = NA_character_) {
      if (tool == "IRF") {
        # self$IRF.compute(hor = hor)
        IRF.out <- self$IRF.compute.batch(hor = hor, prob = prob)
        # IRF.out$base <- self$IRF
        return(IRF.out)
      } else if (tool == "FECD") {
        # self$FEVD.compute(hor = hor)
        FEVD.out <- self$FEVD.compute.batch(hor = hor, prob = prob)
        # FEVD.out$base <- self$FEVD
        return(FEVD.out)
      } else if (tool == "HDC") {
        # HDC.base <- self$HDC.compute(start, end, which.var)
        HDC.out <- self$HDC.compute.batch(start, end, which.var, prob)
        # HDC.boot$base <- HDC.base
        return(HDC.out)
      } else {
        stop("* Please specify correct VAR tool name.")
      }
    },
    IRF.compute.NSR = function(hor = NA_integer_, prob = NA_real_) {
      msg <- paste("* IRF of NSR computed and", prob, "HDP get, time usage")
      HDP <- c((1 - prob) / 2, 1 - (1 - prob) / 2)
      tic(msg)
      # self$Psi.set <- IRF_compute_batch(self$beta.NSR, self$B.NSR, hor)$Psi
      IRF.set <- IRF_compute_batch(self$beta.NSR, self$B.NSR, hor)$IRF
      IRF.med <- apply(IRF.set, c(1, 2), median)
      IRF.avg <- apply(IRF.set, c(1, 2), mean)
      IRF.ub <- apply(IRF.set, c(1, 2), function(x) quantile(x, probs = HDP[2]))
      IRF.lb <- apply(IRF.set, c(1, 2), function(x) quantile(x, probs = HDP[1]))
      toc()
      IRF.list <- list("avg" = IRF.avg, "med" = IRF.med, "ub" = IRF.ub, "lb" = IRF.lb)
      IRF.out <- lapply(IRF.list, mat2cube)
      return(IRF.out)
    }
  ),
  private = list(
    conjugate.post = function() {
      X <- self$X
      K <- dim(X)[2]
      T <- dim(X)[1]
      N0 <- matrix(0, K, K)
      B_0 <- self$beta
      S_0 <- self$Sigma
      v0 <- 0
      NT <- N0 + t(X) %*% X
      BbarT <- solve(NT) %*% (N0 %*% B_0 + (t(X) %*% X) %*% B_0)
      b_post <- as.vector(BbarT)
      vT <- T + v0
      ST <- (v0 / vT) * S_0 + (T / vT) * S_0 + (1 / vT) * t(B_0 - B_0) %*% N0 %*% solve(NT) %*% (t(X) %*% X) %*% (B_0 - B_0)
      S_post <- vT * ST
      v_post <- vT
      invS_post <- solve(S_post)
      Sigma_draw <- solve(rWishart(1, v_post, invS_post)[, , 1])
      V_B_post <- Sigma_draw %x% solve(NT)
      return(
        list(
          "A.bar" = b_post,
          "Sigma.alpha.bar" = V_B_post,
          "S" = S_post,
          "tau" = v_post
        )
      )
    },
    get.restrictions = function(SR = NA, EBR = NA, NSR = NA) {
      restrictions <- list()
      max.irf <- 0
      flag <- 1
      is.NSR <- FALSE
      if (any(!is.na(SR))) {
        for (i in seq_along(SR)) {
          len.var <- length(SR[[i]][[2]])
          shock <- SR[[i]][[1]]
          h <- SR[[i]][[3]]
          if (h > max.irf) max.irf <- h
          sign <- SR[[i]][[4]]
          for (j in 1:len.var) {
            var <- SR[[i]][[2]][j]
            restrictions[[flag]] <- list(
              "type" = "SR",
              "shock" = which(self$var.names == shock),
              "var" = which(self$var.names == var),
              "h" = h,
              "sign" = sign
            )
            flag <- flag + 1
          }
        }
      }
      if (any(!is.na(EBR))) {
        for (i in seq_along(EBR)) {
          shock <- EBR[[i]][[1]]
          h <- EBR[[i]][[3]]
          if (h > max.irf) max.irf <- h
          mb <- EBR[[i]][[4]]
          lb <- EBR[[i]][[5]]
          if (is.na(mb)) mb <- Inf
          if (is.na(lb)) lb <- -Inf
          var.1 <- EBR[[i]][[2]][1]
          var.2 <- EBR[[i]][[2]][2]
          restrictions[[flag]] <- list(
            "type" = "EBR",
            "shock" = which(self$var.names == shock),
            "var.1" = which(self$var.names == var.1),
            "var.2" = which(self$var.names == var.2),
            "h" = h,
            "max.bound" = mb,
            "low.bound" = lb
          )
          flag <- flag + 1
        }
      }
      if (any(!is.na(NSR))) {
        is.NSR <- TRUE
        for (i in seq_along(NSR)) {
          shock <- NSR[[i]][[1]]
          type <- NSR[[i]][[2]]
          start <- which(self$time == NSR[[i]][[3]]) - self$p.lag
          end <- which(self$time == NSR[[i]][[4]]) - self$p.lag
          if ((end - start) > max.irf) max.irf <- end - start
          if (length(NSR[[i]]) == 5) {
            sign <- NSR[[i]][[5]]
            restrictions[[flag]] <- list(
              "type" = "NSR",
              "NSR.type" = "sign",
              "shock" = which(self$var.names == shock),
              "period" = start:end,
              "sign" = sign
            )
            flag <- flag + 1
          } else {
            var <- NSR[[i]][[5]]
            sign <- NSR[[i]][[6]]
            intensity <- NSR[[i]][[7]]
            restrictions[[flag]] <- list(
              "type" = "NSR",
              "NSR.type" = "contribution",
              "shock" = which(self$var.names == shock),
              "var" = which(self$var.names == var),
              "period" = start:end,
              "sign" = sign,
              "intensity" = intensity
            )
            flag <- flag + 1
          }
        }
      }
      return(list(
        "restrictions" = restrictions,
        "max.irf" = max.irf,
        "is.nsr" = is.NSR
      ))
    }
  )
)
