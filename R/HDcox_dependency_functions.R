#### ### ### ### ##
# plyr: round_any #
#### ### ### ### ##

round_any <- function(x, accuracy, f = round) {
  UseMethod("round_any")
}

round_any.numeric <- function(x, accuracy, f = round) {
  f(x / accuracy) * accuracy
}

round_any.double <- function(x, accuracy, f = round) {
  f(x / accuracy) * accuracy
}

round_any.POSIXct <- function(x, accuracy, f = round) {
  tz <- format(x[1], "%Z")
  xr <- round_any(as.numeric(x), accuracy, f)
  as.POSIXct(xr, origin="1970-01-01 00:00.00 UTC", tz=tz)
}
