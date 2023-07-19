#### ### ### ### ##
# plyr: round_any #
#### ### ### ### ##

round2any <- function(x, accuracy, f = round) {
  #UseMethod("round2any")
  if(class(x)[[1]] %in% "numeric"){
    round2any.numeric(x, accuracy, f)
  }else if(class(x)[[1]] %in% "integer"){
    round2any.integer(x, accuracy, f)
  }else if(class(x)[[1]] %in% "double"){
    round2any.double(x, accuracy, f)
  }else if(class(x)[[1]] %in% "POSIXct"){
    round2any.POSIXct(x, accuracy, f)
  }else if(class(x)[[1]] %in% "difftime"){
    round2any.difftime(x, accuracy, f)
  }else{
    stop(paste0("Class ", class(x) ," not defined."))
  }
}

round2any.numeric <- function(x, accuracy, f = round) {
  f(x / accuracy) * accuracy
}

round2any.integer <- function(x, accuracy, f = round) {
  f(x / accuracy) * accuracy
}

round2any.double <- function(x, accuracy, f = round) {
  f(x / accuracy) * accuracy
}

round2any.difftime <- function(x, accuracy, f = round) {
  f(x / accuracy) * accuracy
}

round2any.POSIXct <- function(x, accuracy, f = round) {
  tz <- format(x[1], "%Z")
  xr <- round2any(as.numeric(x), accuracy, f)
  as.POSIXct(xr, origin="1970-01-01 00:00.00 UTC", tz=tz)
}
