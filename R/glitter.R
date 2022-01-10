
#given a data frame with following columns (required) - note the column names -
#alt
#total
#tcn
#lcn
#alpha = 0.05 - CI band for expected somatic and germline VAFs.
#plot = logical, plot the excpected somatic and germline VAFs and called VAFs along with CIs
#and purity estimate

#Determine called VAF by where it lies in the CI bound region of expected somatic and germline VAFs


#compute the following -
# ccf - cancer cell fraction
# multiplicity
# clonal_status
# ccf bounds
# expected germline VAF, CI
# expected somatic VAF, CI
# germline flag


source("germlineFilter_main.R")
library(plyr)

glitter<-function(mat, purity, alpha=.05, plot=TRUE, call_status=FALSE){
  all<-list()
  if(is.na(purity)){stop("purity can't be missing")}
  for(i in 1:nrow(mat)){

    ccf = get.ccf(mat$alt[i], mat$total[i], mat$tcn[i],mat$lcn[i] ,purity)
    evaf = get.expectedVAF(ccf$multiplicity, mat$tcn[i], purity, mat$alt[i], mat$total[i], alpha)
    all[[i]]<-c(ccf, evaf)
  }

  mat.out<-ldply(all, data.frame)

  if(plot){
    n = nrow(mat.out)
    uppers = as.numeric(sapply(mat.out$vafsCI, function(x) get.confBound(x)$xx.u))

    lowers = as.numeric(sapply(mat.out$vafsCI, function(x) get.confBound(x)$xx.l))

    upperg = as.numeric(sapply(mat.out$vafgCI, function(x) get.confBound(x)$xx.u))

    lowerg = as.numeric(sapply(mat.out$vafgCI, function(x) get.confBound(x)$xx.l))


    plot(1:n,mat.out$vafg, type="o",ylim=c(0,1), bty="l", xlab="mutations index", ylab="vaf", cex.axis =1.2, cex.lab=1.2)
    lines(1:n,upperg, type="l",lty=2)
    lines(1:n,lowerg, type="l", lty=2)

    lines(1:n,mat.out$vafs, type="o", col="cadetblue")
    lines(1:n,uppers, type="l",lty=2, col="cadetblue")
    lines(1:n,lowers, type="l", lty=2, col="cadetblue")
    pch.call_status=8
    gfilter1 = gfilter2 = rep(NA, n); gfilter1[which(mat.out$flag.germline)] = which(mat.out$flag.germline); gfilter2[which(mat.out$flag.germline2)]=which(mat.out$flag.germline2)
    if(call_status){pch.call_status = rep(3,n); pch.call_status[which(mat$call_status == "FP")] = 13
    points(gfilter1, rep(0.01,n), pch=4, col="green"); points(gfilter2, rep(0.05,n), pch=4, col="darkgreen") }

    lines(1:n, mat$alt/mat$total, type="o", lwd=2, col="indianred", pch=pch.call_status)

    if(call_status){
      legend("topright", c("expected germline", "expected somatic", "observed VAF (TP)", "observed VAF (FP)", paste0((1-alpha)*100, "% CI")), col=c("black", "cadetblue", "indianred", "indianred", "grey"), lty=c(1,1,1,1,2), bty="n", pch=c(1,1,3,13,NA), lwd=2)}

    if(!(call_status)){
      legend("topright", c("expected germline", "expected somatic", "observed VAF",paste0((1-alpha)*100, "% CI")), col=c("black", "cadetblue", "indianred", "grey"), lty=c(1,1,1,2), bty="n", pch=c(1,1,8,NA), lwd=2)}
  }

  return(mat.out)

}
