###########################################
# This is a file that prints help information, or function options
#    1. ShowCutOptions()
#    2. Info on functions available through Analysis.R

ShowCutOptions <- function() {
  Year<-unique(IDF[,"Year"])
  School<-unique(IDF[,"School"])
  #Class<-unique(IDF[,"Class"])
  print("Years that can be cut:")
  print(Year[order(Year)])
  print("Schools that can be cut:")
  print(School[order(School)])
}
