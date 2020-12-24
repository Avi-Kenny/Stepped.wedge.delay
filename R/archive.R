
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



###################.
##### Section #####
###################.
