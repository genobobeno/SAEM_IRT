ParamCheck <-
function(name,plist,dplist) {
    if (!(tolower(name) %in% tolower(names(dplist)) && !is.na(plist[name]))) {
      plist[name] = dplist[name]
    }
    names(plist)[names(plist)==name] <- names(dplist)[which(tolower(names(dplist)) %in% tolower(name))]
    return(plist)
  }
