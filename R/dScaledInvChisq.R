dScaledInvChisq<-function (x, df, S){
    tmp <- dchisq(S/x, df = df)/(x^2)
    return(tmp)
}