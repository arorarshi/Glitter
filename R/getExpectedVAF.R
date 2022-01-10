#Calculate expected germline and somatic VAFs along with their CI given tcn and multiplicity

getExpectedVAF<-function(multiplicity, tcn, purity, alt, total, alpha = 0.05) {

  if(is.na(multiplicity) | is.na(tcn) | is.na(purity) | is.na(alt) | is.na(total)){message("input value missing"); return(list(vafg=NA, vafgCI=NA_character_,vafs=NA,vafsCI=NA_character_, flag.germline=NA, flag.germline2=NA))}

  vafg = ((purity*multiplicity)+ (1-purity))/((purity*tcn)+2*(1-purity))
  vafs = (purity*multiplicity)/((purity*tcn)+2*(1-purity))

  z = abs(qnorm(alpha/2))

  sigmag=sqrt(vafg*(1-vafg)/total)
  upperg=round(vafg + (z*sigmag),2); lowerg=round(vafg - (z*sigmag),2)

  sigmas=sqrt(vafs*(1-vafs)/total)
  uppers=round(vafs + (z*sigmas),2); lowers= round(vafs - (z*sigmas),2)
  flag.germline = ifelse(round((alt/total),2) >= lowerg, "germline", "pass")
  #flag.germline = ifelse( ((lowerg <= round((alt/total),2)) ), TRUE, FALSE)
  #flag.germline = ifelse( ((lowerg <= round((alt/total),2)) & (uppers < round((alt/total),2))), TRUE, FALSE)

  #flag.germline2 = ifelse( ( ((lowers <= (alt/total)) & ((alt/total) >= uppers)) | ((lowers >= (alt/total)) & ((alt/total) <= uppers))), TRUE, FALSE)
  #flag.germline2 = ifelse( (alt/total) > uppers, TRUE, FALSE)

  return(list(vafg = round(vafg,2), vafgCI = paste0("[", lowerg, "-", upperg, "]"), vafs = round(vafs,2), vafsCI = paste0("[", lowers, "-", uppers, "]"), flag.germline=flag.germline))
}

