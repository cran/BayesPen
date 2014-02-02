metropLambda<-function (tau2, lambda, shape1 = 1.2, shape2 = 1.2, max = 200, ncp = 0)
{
    lambda2 <- lambda^2
    l2_new <- rgamma(rate = sum(tau2)/2, shape = length(tau2),
        n = 1)
    l_new <- sqrt(l2_new)
    logP_old <- sum(dexp(x = tau2, log = TRUE, rate = (lambda2/2))) +
        dbeta(x = lambda/max, log = TRUE, shape1 = shape1, shape2 = shape2) -
        dgamma(shape = sum(tau2)/2, rate = length(tau2), x = (2/lambda2),
            log = TRUE)
    logP_new <- sum(dexp(x = tau2, log = TRUE, rate = (l2_new/2))) +
        dbeta(x = l_new/max, log = TRUE, shape1 = shape1, shape2 = shape2) -
        dgamma(shape = sum(tau2)/2, rate = length(tau2), x = (2/l2_new),
            log = TRUE)
    accept <- (logP_new - logP_old) > log(runif(1))
    if (accept) {
        lambda <- l_new
    }
    return(lambda)
}