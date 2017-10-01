ParamCheck <-
function(plist,dplist) {
  Names = names(dplist)
  for (i in Names) {
    if (!i %in% names(plist)) {
      if (tolower(i) %in% tolower(names(plist))) {
        names(plist)[names(plist)==tolower(i)] <- names(dplist)[names(dplist) %in% i]
      } else {
        plist[i] <- dplist[i]
      }
    } else if (!is.na(plist[i])) {
      if (class(plist[i])!=class(dplist[i])) {
        cat(paste("\nPlease check parameter '",i,
                    "' \n:: Your choice of parameter : ",plist[i]," : has class() ",class(plist[i]),
                  " \n:: The default value is : ",dplist[i]," : and has a class() ",class(dplist[i]),sep=))
        QueryUser("How would you like to remedy this?",c("Change the value"))
      }
    }
    
    return(plist)
  }
