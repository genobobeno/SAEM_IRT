QueryUser<-function(message,choices=NA,type="character",defaultchoice=NA) {
  stopifnot(is.character(message),type %in% c("character","numeric"))
  if (interactive()) {
    cat(message,"\n")
    if (!is.na(choices)) {
      cat(paste("Your options: \n  ",paste(paste("(",1:length(choices),") ",choices,sep=""),collapse="\n   "),"\n "))
      x <- readline("        Enter choice: ")
      if (x == "") {
        return(defaultchoice)
      } else if (tolower(type) != "character") {
        if (grepl("[0-9]",x)) {
          eval(parse(text=paste0("x<-as.",tolower(type),"(x)")))
        } else {
          print("ERROR: Query was expecting a number, not a character. Please start again.")
          return(NA)
        }
      } else {
        if (!grepl("[A-Z|a-z]",x)) {
          if (x=="") {
            if (is.numeric(defaultchoice)) {
              return(choices[defaultchoice]) 
            } else { 
              return(defaultchoice) 
            }
          } else if (grepl("[0-9]",x)) {
            if (as.numeric(x)%in%1:length(choices)) {
              return(choices[as.numeric(x)])
            } else {
              cat("\n\n That entry does not exist.")
              return(NA)
            }
          } else {
            cat("\n\n Whatever you entered in makes no sense.")
            if (is.numeric(defaultchoice)) {
              return(choices[defaultchoice]) 
            } else { 
              return(defaultchoice) 
            }
          }
        } else {
          cat("\n\n Whatever you entered in makes no sense.")
          if (is.numeric(defaultchoice)) {
            return(choices[defaultchoice]) 
          } else { 
            return(defaultchoice) 
          }
        }
      }
    } else {
      x <- readline("Enter a value: ")
      if (x == "") {
        return(defaultchoice)
      } else if (tolower(type) != "character") {
        if (grepl("[0-9]",x)) {
          eval(parse(text=paste0("x<-as.",tolower(type),"(x)")))
        } else {
          print("ERROR: Query was expecting a number, not a character. Please start again.")
          return(NA)
        }
      } else {
        if (!grepl("[A-Z|a-z|,]",x)) {
          print("ERROR: Query was expecting a character, not a number. Please start again.")
          return(NA)
        }
      }
    }
    if (eval(parse(text=paste("is.",tolower(type),"(x)",sep="")))) {
      return(x)
    } else {
      cat("\n User Input does not match respective Type [default:character]")
    }
  } else {
    cat("R is not running in 'interactive()' mode: returning the defaultchoice")
    return(defaultchoice)
  }
}