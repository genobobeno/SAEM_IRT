apaformat <- function(val,digits=2) { 
  sub("^(-?)0\\.", "\\1\\.", sprintf(paste0("%.",digits,"f"), round(val,digits=digits))) 
}
