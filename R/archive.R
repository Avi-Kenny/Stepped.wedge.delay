
###########################################.
##### Old methods from run_analysis() #####
###########################################.

# 1. Two-stage methods
if (analysis$type %in% c("2S LM", "2S GEE", "2S LMM", "2S HL")) {

  if (analysis$type=="2S LM") {

    # Run linear model
    if (data_type=="normal") {
      model <- lm(
        y ~ factor(j) + factor(l),
        data = data$data
      )
    } else if (data_type=="binomial") {
      model <- glm(
        y ~ factor(j) + factor(l),
        data = data$data,
        family = binomial(link="log")
      )
    }

    # Extract coefficients and SEs
    coeff_names <- names(model$coefficients)
    theta_l_hat <- as.numeric(model$coefficients)
    sigma_l_hat <- vcov(model)

  }

  if (analysis$type=="2S GEE") {

    if (data_type=="normal") {
      model <- geeglm(
        y ~ factor(j) + factor(l),
        data = data$data,
        id = i,
        family = "gaussian",
        corstr = analysis$params$corr
      )
    } else if (data_type=="binomial") {
      model <- geeglm(
        y ~ factor(j) + factor(l),
        data = data$data,
        id = i,
        family = binomial(link="log"),
        corstr = corstr
      )
    }

    # Extract estimates and covariance matrix
    coeff_names <- names(model$coefficients)
    theta_l_hat <- as.numeric(model$coefficients)
    sigma_l_hat <- model$geese$vbeta

  }

  if (analysis$type=="2S LMM") {

    if (analysis$params$REML && data_type=="binomial") {
      stop("REML is not well-defined for binomial GLMMs")
    }

    if (data_type=="normal") {
      model <- lmer(
        y ~ factor(j) + factor(l) + (1|i),
        data = data$data,
        REML = analysis$params$REML
      )
    } else if (data_type=="binomial") {
      model <- glmer(
        y ~ factor(j) + factor(l) + (1|i),
        data = data$data,
        family = binomial(link="log")
      )
    }

    # Extract estimates and covariance matrix
    coeff_names <- names(summary(model)$coefficients[,1])
    theta_l_hat <- as.numeric(summary(model)$coefficients[,1])
    sigma_l_hat <- vcov(model)

  }

  if (analysis$type=="2S HL") {

    # Run linear model
    if (data_type=="normal") {

      model <- glmmTMB(
        y ~ factor(j) + factor(l) + (1|i),
        data = data$data
      )

    } else if (data_type=="binomial") {

      model <- glmmTMB(
        y ~ factor(j) + factor(l) + (1|i),
        data = data$data,
        family = binomial(link="log"),
        REML = TRUE,
        start = list(beta=rep(-1,2*L$n_time_points-1)) # !!!!! Avoids a gradient error; perhaps implement this in the 2S LMM as well
      )
      # as.numeric(diag(sigma_l_hat))/as.numeric(diag(sigma_l_hat2)) # !!!!! Directly compare SE estimates with REML=FALSE

    }

    # Extract coefficients and SEs
    coeff_names <- dimnames(summary(model)$coefficients$cond)[[1]]
    theta_l_hat <- as.numeric(summary(model)$coefficients$cond[,1])
    sigma_l_hat <- unname(vcov(model)[[1]])

  }

  # Truncate theta_l_hat and sigma_l_hat
  indices <- c(1:length(coeff_names))[str_sub(coeff_names,1,9)=="factor(l)"]
  coeff_names <- coeff_names[indices]
  theta_l_hat <- theta_l_hat[indices]
  sigma_l_hat <- sigma_l_hat[indices,indices]
  sigma_l_hat <- as.matrix(sigma_l_hat)

  # Generate ML estimates of theta and d
  nll <- function(par) {
    return (
      neg_log_lik(
        theta = par[1],
        d = par[2],
        J = L$n_time_points,
        theta_l_hat = theta_l_hat,
        sigma_l_hat = sigma_l_hat
      )
    )
  }
  opt <- optim(par=c(theta=-0.1, d=1), fn=nll)
  theta_hat <- opt$par[["theta"]]
  d_hat <- opt$par[["d"]]
  h <- optimHess(par=c(theta_hat, d_hat), fn=nll)

  # Use hessian to estimate SEs and extract standard errors
  # Note: tryCatch necessary because sometimes estimated variances are <0
  se_theta_hat <- NA
  tryCatch(
    expr = {
      se_theta_hat <- sqrt(solve(h)[1,1])
    },
    error = function(cond) {}, # !!!!! catch/process error
    warning = function(cond) {} # !!!!! catch/process error
  )

  return (list(
    theta_hat = theta_hat,
    se_theta_hat = se_theta_hat
  ))

}

# 2. Last data point
if (analysis$type=="Last") {

  # !!!!! This is possibly incomplete

  # Run GLMM
  if (data_type=="normal") {
    model <- lmer(
      y ~ factor(j) + factor(l) + (1|i),
      data = data$data
    )
  } else if (data_type=="binomial") {
    model <- glmer(
      y ~ factor(j) + factor(l) + (1|i),
      data = data$data,
      family = binomial(link="log")
    )
  }

  # Extract estimates and covariance matrix
  coeff_names <- names(summary(model)$coefficients[,1])
  theta_l_hat <- as.numeric(summary(model)$coefficients[,1])
  sigma_l_hat <- vcov(model)

  # Truncate theta_l_hat and sigma_l_hat
  indices <- c(1:length(coeff_names))[str_sub(coeff_names,1,9)=="factor(l)"]
  coeff_names <- coeff_names[indices]
  theta_l_hat <- theta_l_hat[indices]
  sigma_l_hat <- sigma_l_hat[indices,indices]
  sigma_l_hat <- as.matrix(sigma_l_hat)

  # Use last theta_l_hat as estimate
  theta_hat <- theta_l_hat[length(theta_l_hat)]
  se_theta_hat <- sqrt(sigma_l_hat[nrow(sigma_l_hat),nrow(sigma_l_hat)])

  return (list(
    theta_hat = theta_hat,
    se_theta_hat = se_theta_hat
  ))

}

# 3. Washout period
if (analysis$type=="WASH") {

  discard <- c(1:analysis$params$length)

  # Filter data based on washout period
  data_filtered <- data$data %>% filter(!(l %in% discard))

  # Run GLMM
  if (data_type=="normal") {
    model <- lmer(
      y ~ factor(j) + x_ij + (1|i),
      data = data_filtered
    )
  } else if (data_type=="binomial") {
    model <- glmer(
      y ~ factor(j) + x_ij + (1|i),
      data = data_filtered,
      family = binomial(link="log")
    )
  }

  # Extract estimate and SE
  theta_hat <- summary(model)$coefficients["x_ij",1]
  se_theta_hat <- summary(model)$coefficients["x_ij",2]

  return (list(
    theta_hat = theta_hat,
    se_theta_hat = se_theta_hat
  ))

}

# 4. Smoothing spline (LTE)
if (analysis$type=="SS") {

  J <- L$n_time_points

  if (analysis$params$t==1) {
    model <- gamm(
      y ~ factor(j) + s(l, k=J, fx=FALSE, bs="cr", m=2, pc=0),
      random = list(i=~1),
      data = data$data
    )
  } else if (analysis$params$t==2) {
    model <- gamm(
      y ~ s(j, k=J, fx=FALSE, bs="cr", m=2, pc=0) +
        s(l, k=J, fx=FALSE, bs="cr", m=2, pc=0),
      random = list(i=~1),
      data = data$data
    )
  }

  theta_hat <- predict(model$gam, newdata=list(j=1, l=J-1), type = "terms")[2]
  se_theta_hat <- summary(model$gam)$se[[length(summary(model$gam)$se)]]

  return (list(
    theta_hat = theta_hat,
    se_theta_hat = se_theta_hat
  ))

}

# 5. Failed attempt to use mono.con with gamm()
if (analysis$type=="MSS") {

  J <- L$n_time_points

  # !!!!! Using mono.con()
  # https://stats.stackexchange.com/questions/197509/how-to-smooth-data-and-force-monotonicity
  # https://r.789695.n4.nabble.com/Use-pcls-in-quot-mgcv-quot-package-to-achieve-constrained-cubic-spline-td4660966.html
  # !!!!! The second link helps fit a point constraint; try this out
  {
    # Modification of example using gamm() object instead
    # !!!!! This works
    x <- sort(runif(100)*4-1)
    f <- exp(4*x)/(1+exp(4*x))
    y <- f+rnorm(100)*0.1
    plot(x,y)
    dat <- data.frame(x=x,y=y)
    f.ug2 <- gamm(y~s(x,k=10,bs="cr")) # Show regular spline fit (and save fitted object); small values of k = more smoothing
    lines(x,fitted(f.ug2$gam))
    sm <- smoothCon(s(x,k=10,bs="cr"),dat,knots=NULL)[[1]]
    FF <- mono.con(sm$xp, up=TRUE) # Monotonicity constraints
    M <- list(
      X = sm$X, # Design matrix
      C = matrix(0,0,0), # No equality constraints
      sp = f.ug2$gam$sp, # Initial guess for param estimate #1
      p = sm$xp, # Initial guess for param estimate #1
      y = y, # Outcome
      w = y*0 + 1, # Weights; vector of 1s of length(y)
      Ain = FF$A, # Inequality constraint matrix to enforce monotonicity
      bin = FF$b, # Inequality constraint vector to enforce monotonicity
      S = sm$S, # Smoothness penalty matrix for cubic spline
      off = 0 # Location of offsets in the penalty matrix
    )
    p <- pcls(M) # Fit spine using penalized constrained least squares
    fv <- Predict.matrix(sm,data.frame(x=x))%*%p # using the monotone spline sm, predict values for the vector x
    lines(x,fv,col=2)

    # !!!!! Adapt this to SW data
    # !!!!! Code runs but gives incorrect results
    # !!!!! One option is to use scam() without a random effect and try to manually produce results
    f.ug2 <- gamm(
      y ~ factor(j) + s(l, k=J, fx=FALSE, bs="cr", m=2, pc=0), # !!!!! try J-1 or 2
      random = list(i=~1),
      data = data$data
    ) # Show regular spline fit (and save fitted object); small values of k = more smoothing

    # plot(data$data$j,data$data$y)
    # lines(data$data$j,fitted(f.ug2$gam))
    sm <- smoothCon(
      s(l, k=J, fx=FALSE, bs="cr", m=2, pc=0), # !!!!! try J-1 or 2
      data = data$data,
      knots = NULL
    )[[1]]
    FF <- mono.con(sm$xp, up=TRUE) # Monotonicity constraints !!!!! switch up to FALSE ?????
    M <- list(
      X = sm$X, # Design matrix
      C = matrix(0,0,0), # No equality constraints
      sp = f.ug2$gam$sp, # Initial guess for param estimate #1
      p = sm$xp, # Initial guess for param estimate #1
      y = data$data$y, # Outcome
      w = data$data$y*0 + 1, # Weights; vector of 1s of length(y)
      Ain = FF$A, # Inequality constraint matrix to enforce monotonicity
      bin = FF$b, # Inequality constraint vector to enforce monotonicity
      S = sm$S, # Smoothness penalty matrix for cubic spline
      off = 0 # Location of offsets in the penalty matrix
    )
    p <- pcls(M) # Fit spline using penalized constrained least squares
    data_predict <- data.frame(
      i = rep(1,6),
      j = c(1:6),
      # k = rep(1,6),
      l = c(0:5)
      # l = rep(0,6)
    )
    fv <- Predict.matrix(sm,data_predict)%*%p # using the monotone spline sm, predict values for the vector x
    lines(x,fv,col=2)

  }

  # # !!!!! Same as above but use scam() instead of gamm()
  # {
  #   model <- gamm(
  #     y ~ s(j, k=J, fx=FALSE, bs="cr", m=2, pc=0) +
  #       s(l, k=J, fx=FALSE, bs="cr", m=2, pc=0),
  #     random = list(i=~1),
  #     data = data$data
  #   )
  #   theta_l_hat <- sapply(c(1:(J-1)), function(l) {
  #     predict(model$gam, newdata=list(j=1, l=l), type = "terms")[2]
  #   })
  # }

  model <- scam(
    y ~ s(j, k=J, fx=FALSE, bs="cr", m=2, pc=0) +
      s(l, k=J, fx=FALSE, bs="cr", m=2, pc=0),
    # random = list(i=~1),
    data = data$data
  )
  theta_l_hat <- sapply(c(1:(J-1)), function(l) {
    predict(model, newdata=list(j=1, l=l), type = "terms")[2]
    # predict(model$gam, newdata=list(j=1, l=l), type = "terms")[2]
  })


  # !!!!! No random effect; two smooth terms
  model <- scam(
    y ~ s(j, k=J, fx=FALSE, bs="cr", m=2, pc=0) +
      s(l, k=J, m=2, bs="mpd"),
    data = data$data
  )

  l_first <- predict(model, newdata=list(j=1, l=0), type = "terms")[[2]]
  l_last <- predict(model, newdata=list(j=1, l=J-1), type = "terms")[[2]]
  theta_hat <- l_last - l_first
  se_theta_hat <- summary(model)$se[[length(summary(model)$se)]]

  # !!!!! Plot estimates
  if (FALSE) {
    # d1 <- sapply(seq(0,6,0.1), function(x) {ifelse(x>0,1,0) * (1-exp(-x/1.4))})
    d1 <- sapply(seq(0,6,0.1), function(x) {(-1/16)*x^2 + (1/2)*x})
    ggplot(
      data.frame(
        x = c(seq(0,6,0.1),c(1:6)),
        y = c(log(0.5)*d1, theta_l_hat),
        fn = c(rep("Parabola",61), rep("Estimates",6))
      ),
      aes(x=x, y=y, color=fn)
    ) +
      geom_point() +
      labs(x="Time (steps)", y="% effect achieved", title="Delay models")
  }

  return (list(
    theta_hat = theta_hat,
    se_theta_hat = se_theta_hat
  ))

  # !!!!! Testing: END !!!!!

}

# 6. This version of SPL calculated the LTE
{
  # if (analysis$type=="SPL") {
  #
  #   # Dynamically build formula
  #   # !!!!! monotonic spline currently only works with lm/glm
  #   if (analysis$params$mono) {
  #     formula <- "y ~ factor(j) + t0"
  #   } else {
  #     formula <- "y ~ factor(j) + (1|i) + t0"
  #   }
  #
  #   # Add spline covariates to dataset
  #   k <- analysis$params$knots
  #   df <- data$data
  #   df$t0 <- df$l
  #   if (length(k)>1) {
  #     for (i in 1:(length(k)-1)) {
  #       df[paste0("t",i)] <- pmax(0,df$l-k[i]) - pmax(0,df$l-k[length(k)])
  #       formula <- paste0(formula," + t",i)
  #     }
  #   }
  #
  #   # !!!!! monotonic spline currently only works with lm/glm
  #   if (analysis$params$mono) {
  #
  #     # Run LM with spline terms
  #     if (data_type=="normal") {
  #       model <- lm(
  #         formula = formula,
  #         data = df
  #       )
  #     } else if (data_type=="binomial") {
  #       model <- glm(
  #         formula = formula,
  #         data = df,
  #         family = binomial(link="log")
  #       )
  #     }
  #
  #     # !!!!! Currently hard-coded for num_times=7
  #     res <- restriktor(
  #       object = model,
  #       constraints = rbind(
  #         c(0,0,0,0,0,0,0,-1,0,0,0,0,0),
  #         c(0,0,0,0,0,0,0,-1,-1,0,0,0,0),
  #         c(0,0,0,0,0,0,0,-1,-1,-1,0,0,0),
  #         c(0,0,0,0,0,0,0,-1,-1,-1,-1,0,0),
  #         c(0,0,0,0,0,0,0,-1,-1,-1,-1,-1,0),
  #         c(0,0,0,0,0,0,0,-1,-1,-1,-1,-1,-1)
  #       ),
  #       rhs = c(0,0,0,0,0,0)
  #     )
  #     s <- summary(res)
  #
  #     # Extract coefficients and SEs
  #     coeff_names <- names(summary(res)$coefficients[,1])
  #     coeffs <- as.numeric(summary(res)$coefficients[,1])
  #     sigma_hat <- summary(res)$V
  #
  #   } else {
  #
  #     # Run GLMM with spline terms
  #     if (data_type=="normal") {
  #       model <- lmer(
  #         formula = formula,
  #         data = df
  #       )
  #     } else if (data_type=="binomial") {
  #       model <- glmer(
  #         formula = formula,
  #         data = df,
  #         family = binomial(link="log")
  #       )
  #     }
  #
  #     # Extract estimates and covariance matrix
  #     coeff_names <- names(summary(model)$coefficients[,1])
  #     coeffs <- as.numeric(summary(model)$coefficients[,1])
  #     sigma_hat <- vcov(model)
  #
  #   }
  #
  #   # Truncate s_hat and sigma_s_hat
  #   indices <- c((length(coeff_names)-length(k)+1):length(coeff_names))
  #   coeff_names <- coeff_names[indices]
  #   coeffs <- coeffs[indices]
  #   sigma_hat <- sigma_hat[indices,indices]
  #   sigma_hat <- as.matrix(sigma_hat)
  #
  #   # Calculate estimators
  #   if (length(k)>1) {
  #     kdiffs <- k[length(k)] - c(0, k[1:(length(k)-1)])
  #     theta_hat <- (kdiffs %*% coeffs)[1,1]
  #     se_theta_hat <- sqrt((kdiffs %*% sigma_hat %*% kdiffs)[1,1])
  #   } else {
  #     theta_hat <- k * coeffs
  #     se_theta_hat <- k * sqrt(sigma_hat[1,1])
  #   }
  #
  #   return (list(
  #     theta_hat = theta_hat,
  #     se_theta_hat = se_theta_hat
  #   ))
  #
  # }
}

# 7. Using gam() and gamm()
if (FALSE) {

  # Random effect
  model <- gamm(
    y ~ factor(j) + s(l, k=J, fx=FALSE, bs="cr", m=2, pc=0),
    random = list(i=~1),
    data = data$data
  )
  predict(model$gam, newdata=list(j=1, l=J-1), type = "terms")[2]
  # Naive
  model <- gam(
    y ~ factor(j) + s(l, k=J, fx=FALSE, bs="cr", m=2, pc=0)
    + s(i, bs="re"), # Random intercept
    data = data$data
  )
  predict(model, newdata=list(i=1, j=1, l=J-1), type = "terms")[2]
}

# 8. Linear spline with hard-coded basis
if (analysis$type=="SPL") {

  # Dynamically build formula
  # !!!!! monotonic spline currently only works with lm/glm
  if (analysis$params$mono) {
    formula <- "y ~ factor(j) + t0"
  } else {
    formula <- "y ~ factor(j) + (1|i) + t0"
  }

  # Add spline covariates to dataset
  k <- analysis$params$knots
  df <- data$data
  df$t0 <- df$l
  if (length(k)>1) {
    for (i in 1:(length(k)-1)) {
      df[paste0("t",i)] <- pmax(0,df$l-k[i]) - pmax(0,df$l-k[length(k)])
      formula <- paste0(formula," + t",i)
    }
  }

  # !!!!! monotonic spline currently only works with lm/glm
  if (analysis$params$mono) {

    # Run LM with spline terms
    if (data_type=="normal") {
      model <- lm(
        formula = formula,
        data = df
      )
    } else if (data_type=="binomial") {
      model <- glm(
        formula = formula,
        data = df,
        family = binomial(link="log")
      )
    }

    # !!!!! Currently hard-coded for num_times=7
    res <- restriktor(
      object = model,
      constraints = rbind(
        c(0,0,0,0,0,0,0,-1,0,0,0,0,0),
        c(0,0,0,0,0,0,0,-1,-1,0,0,0,0),
        c(0,0,0,0,0,0,0,-1,-1,-1,0,0,0),
        c(0,0,0,0,0,0,0,-1,-1,-1,-1,0,0),
        c(0,0,0,0,0,0,0,-1,-1,-1,-1,-1,0),
        c(0,0,0,0,0,0,0,-1,-1,-1,-1,-1,-1)
      ),
      rhs = c(0,0,0,0,0,0)
    )
    s <- summary(res)

    # Extract coefficients and SEs
    coeff_names <- names(summary(res)$coefficients[,1])
    coeffs <- as.numeric(summary(res)$coefficients[,1])
    sigma_hat <- summary(res)$V

  } else {

    # Run GLMM with spline terms
    if (data_type=="normal") {
      model <- lmer(
        formula = formula,
        data = df
      )
    } else if (data_type=="binomial") {
      model <- glmer(
        formula = formula,
        data = df,
        family = binomial(link="log")
      )
    }

    # Extract estimates and covariance matrix
    coeff_names <- names(summary(model)$coefficients[,1])
    coeffs <- as.numeric(summary(model)$coefficients[,1])
    sigma_hat <- vcov(model)

  }

  # Truncate s_hat and sigma_s_hat
  indices <- c((length(coeff_names)-length(k)+1):length(coeff_names))
  coeff_names <- coeff_names[indices]
  coeffs <- coeffs[indices]
  sigma_hat <- sigma_hat[indices,indices]
  sigma_hat <- as.matrix(sigma_hat)

  # Calculate estimators
  if (length(k)>1) {
    kdiffs <- k[length(k)] - c(0, k[1:(length(k)-1)])
    theta_hat <- (kdiffs %*% coeffs)[1,1]
    se_theta_hat <- sqrt((kdiffs %*% sigma_hat %*% kdiffs)[1,1])
  } else {
    theta_hat <- k * coeffs
    se_theta_hat <- k * sqrt(sigma_hat[1,1])
  }

  return (list(
    theta_hat = theta_hat,
    se_theta_hat = se_theta_hat
  ))

}

# 9. Smoothing spline using gamlss()
if (analysis$type=="SS ATE 2") {

  # !!!!! Testing to see if gamlss gives the same estimates as mgcv

  J <- L$n_time_points

  model <- gamlss(
    y ~ factor(j) + random(factor(i)) + pb(l),
    data = data$data
  )

  pred <- lpred(model, what="mu", type="terms", terms="pb(l)", se.fit=TRUE)
  spline_vals <- data.frame(
    "spl_fit" = as.numeric(pred$fit),
    "spl_se" = as.numeric(pred$se.fit),
    "l" = data$data$l
  )
  spline_vals %<>% group_by(l) %>% summarize(
    spl_fit = mean(spl_fit),
    spl_se = mean(spl_se)
  )
  subtract <- (spline_vals %>% filter(l==0))$spl_fit
  spline_vals %<>% mutate(
    spl_fit = spl_fit - subtract
  )
  theta_l_hat <- (spline_vals %>% filter(l!=0))$spl_fit
  se_theta_hats <- (spline_vals %>% filter(l!=0))$spl_se

  # Trapezoid sum: sum of first len-1 values plus half the last value
  # !!!!! Covariances ares currently being ignored; need entire covariance mtx
  len <- length(theta_l_hat)
  A <- matrix(c(rep(1/len,len-1),(1/len)/2), nrow=1)
  ate_hat <- (A %*% theta_l_hat)[1,1]
  se_ate_hat <- sqrt(A %*% diag(se_theta_hats^2) %*% t(A))[1,1]

  return (list(
    ate_hat = ate_hat,
    se_ate_hat = se_ate_hat
  ))

}

# 10. Monotone smoothing spline (can't get entire covariance matrix)
if (analysis$type=="MSS") {

  J <- L$n_time_points

  # !!!!! Main difference is pbm() instead of pb()
  model <- gamlss(
    y ~ factor(j) + random(factor(i)) + pbm(l, mono="down"),
    data = data$data
  )

  pred <- lpred(model, what="mu", type="terms",
                terms="pbm(l, mono = \"down\")", se.fit=TRUE)
  spline_vals <- data.frame(
    "spl_fit" = as.numeric(pred$fit),
    "spl_se" = as.numeric(pred$se.fit),
    "l" = data$data$l
  )
  spline_vals %<>% group_by(l) %>% summarize(
    spl_fit = mean(spl_fit),
    spl_se = mean(spl_se)
  )
  subtract <- (spline_vals %>% filter(l==0))$spl_fit
  spline_vals %<>% mutate(
    spl_fit = spl_fit - subtract
  )
  theta_l_hat <- (spline_vals %>% filter(l!=0))$spl_fit
  se_theta_hats <- (spline_vals %>% filter(l!=0))$spl_se

  # Trapezoid sum: sum of first len-1 values plus half the last value
  # !!!!! Covariances ares currently being ignored; need entire covariance mtx
  len <- length(theta_l_hat)
  A <- matrix(c(rep(1/len,len-1),(1/len)/2), nrow=1)
  ate_hat <- (A %*% theta_l_hat)[1,1]
  se_ate_hat <- sqrt(A %*% diag(se_theta_hats^2) %*% t(A))[1,1]

  return (list(
    ate_hat = ate_hat,
    se_ate_hat = se_ate_hat
  ))

}

# 11. ETI MCMC using JAGS
if (analysis$type=="ETI-MCMC") {

  # !!!!! Only coded for Normal data with J=7

  # !!!!! Transform data
  data_jags <- data$data
  data_jags %<>% dummy_cols(select_columns="j", remove_first_dummy=TRUE)
  data_jags %<>% dummy_cols(select_columns="l", remove_first_dummy=TRUE)

  jags_code <- quote("
      model {
        for (n in 1:N) {
          y[n] ~ dnorm(beta0 + beta_j_2*j_2[n] + beta_j_3*j_3[n] +
          beta_j_4*j_4[n] + beta_j_5*j_5[n] + beta_j_6*j_6[n] + beta_j_7*j_7[n]
          + beta_l_1*l_1[n] + beta_l_2*l_2[n] + beta_l_3*l_3[n] +
          beta_l_4*l_4[n] + beta_l_5*l_5[n] + beta_l_6*l_6[n] + alpha[i[n]],
          1/(sigma^2))
        }
        for (n in 1:I) {
          alpha[n] ~ dnorm(0, 1/(tau^2))
        }
        beta_l_6 ~ dnorm(0, 1.0E-4)
        beta_l_5 ~ dnorm(0, 1.0E-4)
        beta_l_4 ~ dnorm(0, 1.0E-4)
        beta_l_3 ~ dnorm(0, 1.0E-4)
        beta_l_2 ~ dnorm(0, 1.0E-4)
        beta_l_1 ~ dnorm(0, 1.0E-4)
        beta_j_7 ~ dnorm(0, 1.0E-4)
        beta_j_6 ~ dnorm(0, 1.0E-4)
        beta_j_5 ~ dnorm(0, 1.0E-4)
        beta_j_4 ~ dnorm(0, 1.0E-4)
        beta_j_3 ~ dnorm(0, 1.0E-4)
        beta_j_2 ~ dnorm(0, 1.0E-4)
        beta0 ~ dnorm(0, 1.0E-4)
        tau <- 1/sqrt(tau_prec)
        tau_prec ~ dgamma(1.0E-3, 1.0E-3)
        sigma <- 1/sqrt(sigma_prec)
        sigma_prec ~ dgamma(1.0E-3, 1.0E-3)
      }
    ")
  jm <- jags.model(
    file = textConnection(jags_code),
    data = list(
      I = length(unique(data_jags$i)),
      N = nrow(data_jags),
      y = data_jags$y,
      i = data_jags$i,
      j_2 = data_jags$j_2,
      j_3 = data_jags$j_3,
      j_4 = data_jags$j_4,
      j_5 = data_jags$j_5,
      j_6 = data_jags$j_6,
      j_7 = data_jags$j_7,
      l_1 = data_jags$l_1,
      l_2 = data_jags$l_2,
      l_3 = data_jags$l_3,
      l_4 = data_jags$l_4,
      l_5 = data_jags$l_5,
      l_6 = data_jags$l_6
    ),
    quiet = FALSE,
    n.chains = 1,
    n.adapt = 1000
  )
  output <- coda.samples(
    model = jm,
    variable.names = c("beta0", "beta_j_2", "beta_j_3", "beta_j_4", "beta_j_5", "beta_j_6", "beta_j_7", "beta_l_1", "beta_l_2", "beta_l_3", "beta_l_4", "beta_l_5", "beta_l_6", "sigma", "tau"),
    n.iter = 1000,
    thin = 1
  )

  # Extract means
  theta_l_hat <- c()
  for (i in 1:6) {
    theta_l_hat[i] <- summary(output)$statistics[paste0("beta_l_",i),"Mean"]
  }

  # Construct covariance matrix of l terms
  sigma_l_hat <- matrix(NA, nrow=6, ncol=6)
  n_samp <- length(output[[1]][,1])
  for (i in 1:6) {
    for (j in 1:6) {
      sigma_l_hat[i,j] <- cov(
        output[[1]][1:n_samp,paste0("beta_l_",i)],
        output[[1]][1:n_samp,paste0("beta_l_",j)]
      )
    }
  }

  return (res(theta_l_hat,sigma_l_hat))

}

# 12. Lasso (not fully fleshed out)
if (method$method=="Lasso") {

  J <- L$n_time_points

  # # If n_extra_time_points>0, recode l terms
  # if (data$params$n_extra_time_points>0) {
  #   data$data %<>% mutate(
  #     l = ifelse(j>J, J-1, l)
  #   )
  # }

  # # If effect_reached>0, recode l terms
  # if (method$effect_reached>0) {
  #   data$data %<>% mutate(
  #     l = ifelse(l>method$effect_reached, method$effect_reached, l)
  #   )
  # }

  data$data %<>% mutate(
    s_1 = l,
    s_2 = pmax(0,l-1),
    s_3 = pmax(0,l-2),
    s_4 = pmax(0,l-3),
    s_5 = pmax(0,l-4),
    s_6 = pmax(0,l-5)
  )


  # Run GLMM
  if (data_type=="normal") {
    model <- glmmLasso(
      fix = y ~ factor(j) + s_1 + s_2 + s_3 + s_4 + s_5 + s_6,
      rnd = list(i=~1),
      data = data$data,
      lambda = 1 # !!!!!
    )





  }
  # } else if (data_type=="binomial") {
  #   model <- glmer(
  #     y ~ factor(j) + factor(l) + (1|i),
  #     data = data$data,
  #     family = "binomial"
  #   )
  # }

  # coeff_names <- names(summary(model)$coefficients[,1])
  # theta_l_hat <- as.numeric(summary(model)$coefficients[,1])
  # sigma_l_hat <- vcov(model)
  # indices <- c(1:length(coeff_names))[str_sub(coeff_names,1,9)=="factor(l)"]
  # coeff_names <- coeff_names[indices]
  # theta_l_hat <- theta_l_hat[indices]
  # # theta_l_hat # !!!!!
  # sigma_l_hat <- sigma_l_hat[indices,indices]
  # sigma_l_hat <- as.matrix(sigma_l_hat)
  #
  # return (res(theta_l_hat,sigma_l_hat,method$effect_reached))




  # Calculate theta_l_hat vector and sigma_l_hat matrix
  B = rbind(
    c(1,0,0,0,0,0),
    c(2,1,0,0,0,0),
    c(3,2,1,0,0,0),
    c(4,3,2,1,0,0),
    c(5,4,3,2,1,0),
    c(6,5,4,3,2,1)
  )
  theta_l_hat <- B %*% beta_s_hat
  sigma_l_hat <- B %*% sigma_s_hat %*% t(B)

  return (res(theta_l_hat,sigma_l_hat,method$effect_reached))



}




#########################.
##### Trapezoid sum #####
#########################.

# Trapezoid sum: sum of first len-1 values plus half the last value
len <- length(theta_l_hat)
A <- matrix(c(rep(1/len,len-1),(1/len)/2), nrow=1)
ate_hat <- (A %*% theta_l_hat)[1,1]
se_ate_hat <- sqrt(A %*% sigma_l_hat %*% t(A))[1,1]



#############################################################.
##### MAIN: Reproduce table 3.1 (Granston dissertation) #####
#############################################################.

if (run_tables3) {

  # Set levels
  sim %<>% set_levels(
    n_clusters = c(12,24,48),
    n_time_points = 7,
    n_ind_per_cluster = 100,
    theta = log(0.5),
    tau = 0,
    sigma = 0.3,
    data_type = c("normal", "binomial"),
    analysis = c("2S LM", "IG LM"), # !!!!! update
    delay_model = list(
      "Exp (d=0)" = list(type="exp", params=list(d=0)),
      "Exp (d=0.5)" = list(type="exp", params=list(d=0.5)),
      "Exp (d=1.4)" = list(type="exp", params=list(d=1.4))
    )
  )

  # Run simulation and save output
  sim %<>% run("one_simulation", sim_uids=.tid)
  saveRDS(sim, file=paste0("../simba.out/sim_",.tid,".simba"))

  # Output results
  if (run_results) {
    sim <- readRDS("../simba.out/sim_tab3.1.simba")
    print(summary(
      sim_obj = sim,
      bias = list(name="bias_theta", truth="theta", estimate="theta_hat")
    ))
  }

}



#############################################################.
##### MAIN: Reproduce table 3.2 (Granston dissertation) #####
#############################################################.

if (run_tables3) {

  # Set levels
  sim %<>% set_levels(
    n_clusters = 24,
    n_time_points = c(5,7,9),
    n_ind_per_cluster = 100,
    theta = log(0.5),
    tau = 0,
    sigma = 0.3,
    data_type = c("normal", "binomial"),
    analysis = c("2S LM", "IG LM"), # !!!!! update
    delay_model = list(
      "Exp (d=0)" = list(type="exp", params=list(d=0)),
      "Exp (d=0.5)" = list(type="exp", params=list(d=0.5)),
      "Exp (d=1.4)" = list(type="exp", params=list(d=1.4))
    )
  )

  # Run simulation and save output
  sim %<>% run("one_simulation", sim_uids=.tid)
  saveRDS(sim, file=paste0("../simba.out/sim_",.tid,".simba"))

  # Output results
  if (run_results) {
    sim <- readRDS("../simba.out/sim_tab3.2.simba")
    print(summary(
      sim_obj = sim,
      bias = list(name="bias_theta", truth="theta", estimate="theta_hat")
    ))
  }

}



#############################################################.
##### MAIN: Reproduce table 3.3 (Granston dissertation) #####
#############################################################.

if (run_tables3) {

  # Set levels
  sim %<>% set_levels(
    n_clusters = 24,
    n_time_points = 7,
    n_ind_per_cluster = c(20,50,100),
    theta = log(0.5),
    tau = 0,
    sigma = 0.3,
    data_type = c("normal", "binomial"),
    analysis = c("2S LM", "IG LM"), # !!!!! update
    delay_model = list(
      "Exp (d=0)" = list(type="exp", params=list(d=0)),
      "Exp (d=0.5)" = list(type="exp", params=list(d=0.5)),
      "Exp (d=1.4)" = list(type="exp", params=list(d=1.4))
    )
  )

  # Run simulation and save output
  sim %<>% run("one_simulation", sim_uids=.tid)
  saveRDS(sim, file=paste0("../simba.out/sim_",.tid,".simba"))

  # Output results
  if (run_results) {
    sim <- readRDS("../simba.out/sim_tab3.3.simba")
    print(summary(
      sim_obj = sim,
      bias = list(name="bias_theta", truth="theta", estimate="theta_hat")
    ))
  }

}



#############################################################.
##### MAIN: Reproduce table 3.4 (Granston dissertation) #####
#############################################################.

if (run_tables3) {

  # Set levels
  sim %<>% set_levels(
    n_clusters = c(12,24,48),
    n_time_points = 7,
    n_ind_per_cluster = 100,
    theta = log(0.5),
    tau = 0,
    sigma = 0.3,
    data_type = c("normal", "binomial"),
    analysis = list("2S LM"=list(type="2S LM")),
    delay_model = list(
      "Exp (d=0.5)" = list(type="exp", params=list(d=0.5)),
      "Exp (d=1.4)" = list(type="exp", params=list(d=1.4))
    )
  )

  # Run simulation and save output
  sim %<>% run("one_simulation", sim_uids=.tid)
  saveRDS(sim, file=paste0("../simba.out/sim_",.tid,".simba"))

  # Output results
  if (run_results) {
    sim <- readRDS("../simba.out/sim_tab3.4.simba")
    print(summary(
      sim_obj = sim,
      bias = list(name="bias_d", truth="d", estimate="d_hat")
    ))
  }

}



#############################################################.
##### MAIN: Reproduce table 3.5 (Granston dissertation) #####
#############################################################.

if (run_tables3) {

  # Set levels
  sim %<>% set_levels(
    n_clusters = 24,
    n_time_points = c(5,7,9),
    n_ind_per_cluster = 100,
    theta = log(0.5),
    tau = 0,
    sigma = 0.3,
    data_type = c("normal", "binomial"),
    analysis = list("2S LM"=list(type="2S LM")),
    delay_model = list(
      "Exp (d=0.5)" = list(type="exp", params=list(d=0.5)),
      "Exp (d=1.4)" = list(type="exp", params=list(d=1.4))
    )
  )

  # Run simulation and save output
  sim %<>% run("one_simulation", sim_uids=.tid)
  saveRDS(sim, file=paste0("../simba.out/sim_",.tid,".simba"))

  # Output results
  if (run_results) {
    sim <- readRDS("../simba.out/sim_tab3.5.simba")
    print(summary(
      sim_obj = sim,
      bias = list(name="bias_d", truth="d", estimate="d_hat")
    ))
  }

}



#############################################################.
##### MAIN: Reproduce table 3.6 (Granston dissertation) #####
#############################################################.

if (run_tables3) {

  # Set levels
  sim %<>% set_levels(
    n_clusters = 24,
    n_time_points = 7,
    n_ind_per_cluster = c(20,50,100),
    theta = log(0.5),
    tau = 0,
    sigma = 0.3,
    data_type = c("normal", "binomial"),
    analysis = list("2S LM"=list(type="2S LM")),
    delay_model = list(
      "Exp (d=0.5)" = list(type="exp", params=list(d=0.5)),
      "Exp (d=1.4)" = list(type="exp", params=list(d=1.4))
    )
  )

  # Run simulation and save output
  sim %<>% run("one_simulation", sim_uids=.tid)
  saveRDS(sim, file=paste0("../simba.out/sim_",.tid,".simba"))

  # Output results
  if (run_results) {
    sim <- readRDS("../simba.out/sim_tab3.6.simba")
    print(summary(
      sim_obj = sim,
      bias = list(name="bias_d", truth="d", estimate="d_hat")
    ))
  }

}



#########################################################.
##### COVERAGE: Investigate coverage issue (no GEE) #####
#########################################################.

if (run_cov_nogee) {

  # Set levels
  sim %<>% set_levels(
    n_clusters = c(12,48),
    n_time_points = 7,
    n_ind_per_cluster = seq(20, 100, 20),
    theta = log(0.5),
    tau = c(0,0.25),
    sigma = 3,
    data_type = c("normal", "binomial"),
    analysis = list(
      "2S LM" = list(type="2S LM"),
      "2S LMM REML" = list(type="2S LMM", params=list(REML=TRUE))
    ),
    delay_model = list(
      "Exp (d=0)" = list(type="exp", params=list(d=0)),
      "Exp (d=1.4)" = list(type="exp", params=list(d=1.4))
    )
  )

  # Run simulation and save output
  sim %<>% run("one_simulation", sim_uids=.tid)
  saveRDS(sim, file=paste0("../simba.out/sim_",.tid,".simba"))

  # Output results
  if (run_results) {
    sim <- readRDS("../simba.out/sim_cov_nogee.simba")
    print(summary(
      sim_obj = sim,
      coverage = list(
        name = "cov_theta",
        truth = "theta",
        estimate = "theta_hat",
        se = "se_theta_hat",
        na.rm = TRUE
      )
    ))
  }

}



######################################################.
##### COVERAGE: Investigate coverage issue (GEE) #####
######################################################.

if (run_cov_gee) {

  # Set levels
  sim %<>% set_levels(
    n_clusters = 24,
    n_time_points = 7,
    n_ind_per_cluster = 50,
    theta = log(0.5),
    tau = c(0,0.25),
    sigma = 0.3,
    data_type = c("normal", "binomial"),
    analysis = list(
      "2S GEE EXC" = list(type="2S GEE", params=list(corr="exchangeable")),
      "2S GEE IND" = list(type="2S GEE", params=list(corr="independence"))
    ),
    delay_model = list("Exp"=list(type="exp", params=list(d=1.4)))
  )

  # Run simulation and save output
  sim %<>% run("one_simulation", sim_uids=.tid)
  saveRDS(sim, file=paste0("../simba.out/sim_",.tid,".simba"))

  # Output results
  if (run_results) {
    sim <- readRDS("../simba.out/sim_cov_gee.simba")
    print(summary(
      sim_obj = sim,
      coverage = list(
        name = "cov_theta",
        truth = "theta",
        estimate = "theta_hat",
        se = "se_theta_hat",
        na.rm = TRUE
      )
    ))
  }

}



##############################################################.
##### COVERAGE: Investigate coverage issue (REML vs. ML) #####
##############################################################.

if (run_cov_reml) {

  # Set levels
  sim %<>% set_levels(
    n_clusters = c(12,48),
    n_time_points = 7,
    n_ind_per_cluster = c(20,100),
    theta = log(0.5),
    tau = c(0,0.25),
    sigma = 0.3,
    data_type = "normal",
    analysis = list(
      "2S LMM REML" = list(type="2S LMM", params=list(REML=TRUE)),
      "2S LMM ML" = list(type="2S LMM", params=list(REML=FALSE))
    ),
    delay_model = list("Exp"=list(type="exp", params=list(d=1.4)))
  )

  # Run simulation and save output
  sim %<>% run("one_simulation", sim_uids=.tid)
  saveRDS(sim, file=paste0("../simba.out/sim_",.tid,".simba"))

  # Output results
  if (run_results) {
    sim <- readRDS("../simba.out/sim_cov_reml.simba")
    print(summary(
      sim_obj = sim,
      coverage = list(
        name = "cov_theta",
        truth = "theta",
        estimate = "theta_hat",
        se = "se_theta_hat",
        na.rm = TRUE
      )
    ))
  }

}



###############################################################.
##### COVERAGE: Investigate coverage issue (H-likelihood) #####
###############################################################.

if (run_cov_hlik) {

  # Set levels
  sim %<>% set_levels(
    n_clusters = 12,
    n_time_points = 7,
    n_ind_per_cluster = c(20,100),
    theta = log(0.5),
    tau = c(0,0.25),
    sigma = 3,
    data_type = "binomial",
    analysis = list(
      "2S LM" = list(type="2S LM"),
      "2S HL" = list(type="2S HL")
    ),
    delay_model = list("Exp"=list(type="exp", params=list(d=1.4)))
  )

  # Run simulation and save output
  sim %<>% run("one_simulation", sim_uids=.tid)
  saveRDS(sim, file=paste0("../simba.out/sim_",.tid,".simba"))

  # Output results
  if (run_results) {
    sim <- readRDS("../simba.out/sim_cov_hlik.simba")
    print(summary(
      sim_obj = sim,
      coverage = list(
        name = "cov_theta",
        truth = "theta",
        estimate = "theta_hat",
        se = "se_theta_hat",
        na.rm = TRUE
      )
    ))
  }

}



#########################.
##### TESTING: PMLE #####
#########################.

if (run_testing_pmle) {

  # First, test the case where we have only one observation per cluster

  # Generate dataset
  data <- generate_dataset(
    alpha = log(0.1),
    tau = 0,
    theta = log(0.5),
    n_clusters = 48,
    n_time_points = 7,
    n_ind_per_cluster = 50,
    data_type = "normal",
    sigma = 3,
    delay_model = list(type="exp", params=list(d=1.4))
  )

  # Estimator with mgcv (smooth for Tx effect)
  model <- gamm(
    y ~ factor(j) + s(l, k=7, fx=FALSE, bs="cr", m=2, pc=0),
    random = list(i=~1),
    data = data$data
  )
  est <- predict(model$gam, newdata=list(j=1, l=6), type = "terms")[2]
  se <- summary(model$gam)$se[[length(summary(model$gam)$se)]]

  # Estimator with mgcv (smooths for Tx effect and time)
  model <- gamm(
    y ~ s(j, k=7, fx=FALSE, bs="cr", m=2, pc=0) +
      s(l, k=7, fx=FALSE, bs="cr", m=2, pc=0),
    random = list(i=~1),
    data = data$data
  )
  est <- predict(model$gam, newdata=list(j=1, l=6), type = "terms")[2]
  se <- summary(model$gam)$se[[length(summary(model$gam)$se)]]

}



##################################################.
##### MISC: Check MLE calculation (Granston) #####
##################################################.

if (run_misc) {

  J <- 5
  theta_hat_l <- matrix(1:J,ncol=1)
  mu_d <- matrix(0.3*c(3:(J+2)),ncol=1)
  A <- matrix(runif(J^2)*2-1, ncol=J)
  sigma_inv <- t(A) %*% A

  neg_log_lik <- function(theta) {
    return(
      (1/2) * t(theta_hat_l-(theta*mu_d)) %*%
        sigma_inv %*%
        (theta_hat_l-(theta*mu_d))
    )
  }

  optim(
    par = 6,
    fn = function(x) {
      return ( neg_log_lik(x) )
    },
    method = "BFGS"
  )

  mle_ank <- (
    (t(theta_hat_l) %*% sigma_inv %*% mu_d) +
      t((t(theta_hat_l) %*% sigma_inv %*% mu_d))) /
    (2* t(mu_d) %*% sigma_inv %*% mu_d)
  mle_tsg <- (t(theta_hat_l) %*% sigma_inv %*% mu_d) /
    (t(mu_d) %*% sigma_inv %*% mu_d)

  print(mle_ank)
  print(mle_tsg)

}



####################################################.
##### MISC: Graph of linear spline alternative #####
####################################################.

if (run_misc) {

  ggplot(data.frame(x=c(0,6)), aes(x=x)) +
    stat_function(fun = function(x) {
      return ( 0.7*x + (0.3/5 - 0.7)*pmax(0,x-1) )
    }) +
    geom_point(aes(x=0, y=0)) +
    geom_point(aes(x=6, y=1), colour="green") +
    geom_point(aes(x=1, y=0.7), colour="purple") +
    labs(
      title = "Linear spline model for R_il",
      y = "R_il (% of treatment effect achieved)",
      x = "Time since implementation (l_i)"
    )

}



#############################################################.
##### ARCHIVE: Old code related to two-stage approaches #####
#############################################################.

if (FALSE) {

  # Negative log lik corresponding to two-stage dissertation approach
  sim %<>% add_method(
    "neg_log_lik",
    function(theta, d, J, theta_l_hat, sigma_l_hat) {

      l_times <- 1:(J-1)
      mu_d <- 1-exp(-l_times/d)
      log_lik <- -0.5 * t(theta_l_hat - theta*mu_d) %*%
        solve(sigma_l_hat) %*% (theta_l_hat - theta*mu_d)

      return (-1*log_lik)

    }
  )

  # Negative log lik corresponding to two-stage spline approach
  sim %<>% add_method(
    "neg_log_lik_spl",
    function(theta, p_x, p_y, J, theta_l_hat, sigma_l_hat) {

      # !!!!! g_x hard-coded for now
      g_x <- J

      l_times <- 1:(J-1)
      mu_spl <- sapply(l_times, function(l) {
        I1 <- ifelse(0<l & l<=p_x, 1, 0)
        I2 <- ifelse(p_x<l & l<=g_x, 1, 0)
        I3 <- ifelse(g_x<l, 1, 0)
        (p_y/p_x)*l*I1 + ((1-p_y)*l+g_x+p_y-p_x-1)/(g_x-p_x)*I2 + I3
      })

      log_lik <- -0.5 * t(theta_l_hat - theta*mu_spl) %*%
        solve(sigma_l_hat) %*% (theta_l_hat - theta*mu_spl)

      return (-1*log_lik)

    }
  )



  # Method corresponding to "semiparametric stochastic..." paper

  # Compute the closed-form estimator
  Y <- data$data$y
  I <- data$params$n_clusters
  J <- data$params$n_time_points
  K <- data$params$n_ind_per_cluster
  int_times <- data$data %>%
    group_by(i,j) %>% summarize(l=max(l))
  s <- matrix(int_times$l, nrow=I, byrow=TRUE)

  for (i in 1:I) {
    for (j in 1:J) {
      N_ij <- matrix(0, nrow=K, ncol=J-1)
      N_ij[,s[i,j]] <- 1
      T_ij <- matrix(0, nrow=K, ncol=J-1)
      if (j!=J) { T_ij[,j] <- 1 }
      if (i==1 && j==1) {
        N <- N_ij
        mtx_T <- T_ij
      } else {
        N <- rbind(N,N_ij)
        mtx_T <- rbind(mtx_T,T_ij)
      }
    }
  }

  Q <- matrix(0, nrow=J-1, ncol=J-3)
  for (j in 1:(J-3)) {
    Q[j,j] <- 1
    Q[j+1,j] <- -2
    Q[j+2,j] <- 1
  }
  R <- matrix(0, nrow=J-3, ncol=J-3)
  for (j in 1:(J-3)) {
    if (j>1) { R[j-1,j] <- 1/6 }
    R[j,j] <- 2/3
    if (j<J-3) { R[j+1,j] <- 1/6 }
  }

  # Construct variance component estimators
  V_s <- tau * (t(B_s) %*% B_s) + V
  X_s <- cbind(X, N %*% mtx_T)
  beta_hat_s <- 999
  V_s_inv <- solve(V_s)
  l_R <- (-1/2) * (
    log(det(V_s)) + log(det( t(X_s) %*% V_s_inv * X_s )) +
      t(Y - (X_s %*% beta_hat_s)) %*% V_s_inv %*% (Y - (X_s %*% beta_hat_s))
  )

  K_star <- Q %*% solve(R) %*% t(Q)
  simga2_e <- 0.01 # !!!!!
  sigma2_gamma <- 0.01 # !!!!!
  var_sum <- simga2_e + sigma2_gamma
  W <- diag(rep(1/var_sum,I*J*K))
  V <- diag(rep(var_sum,I*J*K))
  X <- cbind(matrix(1, nrow=I*J*K, ncol=1), mtx_T)
  W_f <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  lambda <- 20 # !!!!!
  f_hat <- solve(t(N) %*% W_f %*% N + lambda * K_star) %*% t(N) %*% W_f %*% Y

  # Plot delay model vs f_hat
  theta <- log(0.5)
  d2 <- sapply(seq(0,6,0.1), function(x) {
    theta * ifelse(x>0,1,0) * (1-exp(-x/1))
  })
  # Plot functions
  ggplot(
    data.frame(
      x = c(c(0:(J-1)), seq(0,6,0.1)),
      y = c(0,f_hat,d2),
      fn = c(rep("f_hat",J),rep("Exp (d=1)",61))
    ),
    aes(x=x, y=y, color=fn)
  ) +
    geom_line() +
    labs(x="Time (steps)", y="Intervention effect")



  print(paste("Number of available cores:", parallel::detectCores()))
  print(paste("SLURM_ARRAY_JOB_ID:", Sys.getenv("SLURM_ARRAY_JOB_ID")))
  print(paste("SLURM_CPUS_ON_NODE:", Sys.getenv("SLURM_CPUS_ON_NODE")))
  print(paste("SLURM_NODELIST:", Sys.getenv("SLURM_NODELIST")))
  print(paste("SLURM_NNODES:", Sys.getenv("SLURM_NNODES")))
  print(paste("SLURM_NTASKS:", Sys.getenv("SLURM_NTASKS")))

  summary(
    sim_obj = sim,
    sd = list(
      list(name="sd_theta_hat", x="theta_hat"),
      list(name="sd_d_hat", x="d_hat")
    ),
    bias = list(
      list(name="bias_theta", truth="theta", estimate="theta_hat"),
      list(name="bias_d", truth="d", estimate="d_hat")
    ),
    coverage = list(
      list(name="cov_theta", truth="theta",
           estimate="theta_hat", se="se_theta_hat"),
      list(name="cov_d", truth="d", estimate="d_hat", se="se_d_hat")
    )
  )

  # Plots
  plot_sw_design(data_1)
  plot_outcome(data_1, type="no error")
  plot_outcome(data_1, type="realized")

  # Binomial GLM
  model_binomial_gee1 <- geeglm(
    y ~ factor(j) + factor(x_ij),
    data = data$data,
    id = i,
    family = binomial(link = "log"),
    corstr = "exchangeable"
  )
  summary(model_binomial_gee1)
  system.time(
    model_binomial_gee2 <- geeglm(
      y ~ factor(j) + factor(l),
      data = data$data,
      id = i,
      family = binomial(link = "log"),
      corstr = "exchangeable"
    )
  )
  summary(model_gee2)

  # Only estimate theta (step function)
  model_normal_gee1 <- geeglm(
    y ~ factor(j) + factor(x_ij),
    data = data$data,
    id = i,
    family = "gaussian",
    corstr = "exchangeable"
  )
  summary(model_normal_gee1)

}



#############################################.
##### ARCHIVE: log likelihood of spline #####
#############################################.

#' Return the log likelihood for the spline model
#'
#' @param sigma_v !!!!! TO DO
#' @return !!!!! TO DO

log_lik_spline <- function(
  sigma_v, sigma_e, alpha,
  beta_1, beta_2, beta_3, beta_4, beta_5,
  theta, p_x, p_y, g_x, data
) {

  I <- data$params$n_clusters
  J <- data$params$n_time_points
  K <- data$params$n_ind_per_cluster
  beta <- c(beta_1, beta_2, beta_3, beta_4, beta_5)

  # !!!!! Can speed this up by calculating this in advance
  df <- data$data %>% mutate(
    r_ij = r_ij(p_x, p_y, g_x, c_i, j),
    s = y - alpha - beta[j] - theta*r_ij,
    t = s^2
  )

  df_sum <- df %>% group_by(i) %>% summarize(
    s = sum(s),
    t = sum(t)
  )

  s_i <- df_sum$s
  t_i <- df_sum$t

  t1 <- (I/2) *
    log((2*pi*(sigma_e^2)*(sigma_v^2))/((sigma_e^2)+(J*K*(sigma_v^2))))
  t2 <- (sigma_v^2)/(2*(sigma_e^2)*((sigma_e^2)+J*K*(sigma_v^2)))
  t3 <- sum(s_i^2-t_i)

  return(t1 + t2*t3)

}



# Helper function to calculate the spline log likelihood (above)
r_ij <- function(p_x, p_y, g_x, c_i, j) {

  I1 <- ifelse(c_i<j & j<=c_i+p_x, 1, 0)
  I2 <- ifelse(c_i+p_x<j & j<=c_i+g_x, 1, 0)
  I3 <- ifelse(c_i+g_x<j, 1, 0)
  (p_y/p_x)*(j-c_i)*I1 + (((c_i-j+1)*(p_y-1))/(g_x-p_x)+1)*I2 + I3

}



#####################################.
##### MAIN: Process sim results #####
#####################################.

if (run_process_results) {

  # Read in simulation object
  # sim <- readRDS("../simba.out/sim_main_1026.simba")

  # Generate true ATE value
  sim$results %<>% mutate(
    ate = ifelse(delay_model=="EXP (d=0)", round(theta*1,4),
                 ifelse(delay_model=="EXP (d=1.4)", round(theta*0.77,4),
                        ifelse(delay_model=="SPL (k=2,4 s=0.1,0.4)", round(theta*0.57,4),
                               999)))
  )

  # Summarize data
  summ <- summary(
    sim_obj = sim,
    mean = list(all=TRUE, na.rm=TRUE),
    quantile = list(
      list(name="q025_ate", x="ate_hat", prob=0.025, na.rm=TRUE),
      list(name="q975_ate", x="ate_hat", prob=0.975, na.rm=TRUE)
    ),
    coverage = list(
      list(
        name = "cov_ate",
        truth = "ate",
        estimate = "ate_hat",
        se = "se_ate_hat",
        na.rm = TRUE
      ),
      list(
        name = "beta",
        truth = 0,
        estimate = "ate_hat",
        se = "se_ate_hat",
        na.rm = TRUE
      )
    )
  )

  # Transform summary data (1)
  summ %<>% mutate(
    method = factor(method, levels=c("HH","ETI","SS")),
    power = 1 - beta
  )

  # Transform summary data (2)
  summ$theta_log <- rep(NA, nrow(summ))
  summ$theta_log <- ifelse(round(as.numeric(summ$theta),1)==-0.7,
                           "log(0.5)", summ$theta_log)
  summ$theta_log <- ifelse(round(as.numeric(summ$theta),1)==-0.5,
                           "log(0.6)", summ$theta_log)
  summ$theta_log <- ifelse(round(as.numeric(summ$theta),1)==-0.4,
                           "log(0.7)", summ$theta_log)
  summ$theta_log <- ifelse(round(as.numeric(summ$theta),1)==-0.2,
                           "log(0.8)", summ$theta_log)
  summ$theta_log <- ifelse(round(as.numeric(summ$theta),1)==-0.1,
                           "log(0.9)", summ$theta_log)
  summ$theta_log <- ifelse(round(as.numeric(summ$theta),1)==0,
                           "log(1.0)", summ$theta_log)
  summ$theta_log <- as.factor(summ$theta_log)

}



