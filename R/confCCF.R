#Confidence interval for CCF at a default alpha of 0.05 or 95% CI
confCCF=function(alt, ref, tcn, minor, purity, multiplicity, alpha=0.05){

  if(is.na(alt) | is.na(ref) | is.na(tcn) | is.na(minor) | is.na(purity) | is.na(multiplicity)){message("input value is missing"); return(list(lower=NA, upper=NA))}

  depth=alt+ref
  ci=binconf(alt, depth, alpha, method=c("exact"))
  lowerci=ci[,2]
  upperci=ci[,3]

  lower=computeCCF(vaf=lowerci, tcn, minor, purity, multiplicity)$ccf
  upper=computeCCF(vaf=upperci, tcn, minor, purity, multiplicity)$ccf

  return(list(lower=lower, upper=upper))

}
