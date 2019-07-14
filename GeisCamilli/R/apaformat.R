apaformat <- function(val,digits=2) { 
  sub("^(-?)0\\.", "\\1\\.", sprintf("%.2f", round(val,digits=digits))) 
}
