sec2ext <- function(x) {
  time <- cbind(c(x %/% 86400, 
                  (x %% 86400) %/% 3600,
                  ((x %% 86400) %% 3600) %/% 60, 
                  round(((x %% 86400) %% 3600) %% 60, 2)),
                c("day", "hour", "min", "sec"))
  paste(
    apply(time[cumsum(time[, 1] > 0) > 0, , drop=FALSE], 1,
          function(x) {
            paste(x[1], " ", x[2], ifelse(x[1] > 1, "s", ""), sep="")
          }), collapse=", ")
}