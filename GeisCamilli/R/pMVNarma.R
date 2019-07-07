pMVNarma<-function(x,pList,thHat,VAR) {
  mvrnormArma(length(pList[[x]]), thHat[pList[[x]],], VAR)
}