percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}


comprss <- function(tx) { 
  div <- findInterval(as.numeric(gsub("\\,", "", tx)), 
                      c(0, 1e3, 1e6, 1e9, 1e12) )  # modify this if negative numbers are possible
  paste(round( as.numeric(gsub("\\,","",tx))/10^(3*(div-1)), 2), 
        c("","K","M","B","T")[div] )}
