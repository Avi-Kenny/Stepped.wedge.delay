
# First parameter: beta should be roughly uniform over (-20,0)
n <- 100000
alpha = rexp(n, rate=10)
# beta = -exp(alpha)
beta = -exp(3-13*alpha)
ggplot(
  data.frame(
    values = c(alpha,beta),
    which = rep(c("alpha","beta"), each=n)
  ), aes(x=values)
) +
  # geom_density() +
  geom_histogram(bins=100) +
  facet_wrap(~which, scales="free")

# Subsequent parameters: beta should be roughly uniform over (-10,10)
# alpha = rbeta(n, shape1=10, shape2=0.85)
n <- 100000
alpha1 = rexp(n, rate=13)
alpha2 = rexp(n, rate=13)
alpha3 = rexp(n, rate=13)
beta1 = - exp(3-13*alpha1)
beta2 = exp(3-13*alpha1) - exp(3-13*alpha2)
beta3 = exp(3-13*alpha2) - exp(3-13*alpha3)
ggplot(
  data.frame(
    # values = c(beta1),
    # which = rep(c("beta1"), each=n)
    values = c(beta1,beta2,beta3),
    which = rep(c("beta1","beta2","beta3"), each=n)
  ), aes(x=values)
) +
  # geom_density() +
  geom_histogram(bins=100) +
  facet_wrap(~which, scales="free")

