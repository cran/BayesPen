dLambda<-function (rate, shape, lambda) {
    tmp <- dgamma(x = I(lambda^2), rate = rate, shape = shape) * 2 * lambda
    return(tmp)
}
