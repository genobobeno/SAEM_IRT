SFileString<-function(l,gen,r=NA) {
  f<-paste0(sapply(names(l),function(x) (paste0(substr(x,1,1),l[[x]]))),collapse="_")
  ifelse(gen,
         ifelse(!is.na(r),paste0("Gen_",f,"_",r),paste0("Gen_",f)),
         ifelse(!is.na(r),paste0("Fit_",f,"_",r),paste0("Fit_",f)))
}