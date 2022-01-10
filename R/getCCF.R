#Wrapper function to get necessary things from CCF

getCCF<-function(alt, total, tcn,lcn, purity){

  if(is.na(tcn) | is.na(lcn) | is.na(alt) | is.na(total) | is.na(purity)){message("input value is missing"); return(list(ccf=NA, multiplicity=NA, clonal_status=NA_character_, ccfCI = NA_character_))}
  fit=computeCCF(vaf=alt/total, tcn=tcn, minor=lcn, purity=purity)
  #cancer cell fraction for the mutation
  ccf=round(fit$ccf,2) #assumes mut on minor copy
  #mutant copy number (this is the v you want)
  multiplicity=fit$multiplicity

  #calculate confidence interval for CCF
  conf=confCCF(alt=alt, ref=total-alt, tcn=tcn, minor=lcn, purity=purity, multiplicity=multiplicity)
  lower.ccf=round(conf$lower,2)
  upper.ccf=round(conf$upper,2)
  #assign clonal status for each mutation: if the lower bound of 95% CI for the observed CCF is >= 75%, call it clonal mutation.
  clonal_status=ifelse(lower.ccf >= 0.75,"clonal","subclonal")
  clonal_status[lower.ccf < 0.75 & ccf >= 0.80]="likely clonal"
  ccfCI = paste0("[",lower.ccf, "-", upper.ccf,"]")
  return(list(ccf=ccf, multiplicity=multiplicity, clonal_status=clonal_status, ccfCI = ccfCI))

}

