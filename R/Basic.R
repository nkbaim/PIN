#' @export
mean_na <- function(x) {
  
  if (is.nan(mean(x, na.rm = T))) {
    return(NA)
  } else {
    return(mean(x, na.rm = T))
  }
  
}

#' @export
sd_na <- function(x) {
  
  if (is.nan(sd(x, na.rm = T))) {
    return(NA)
  } else {
    return(sd(x, na.rm = T))
  }
  
}

#' @export
cv_na <- function(x) {
  
  return(sd(x, na.rm = T) / mean(x, na.rm = T))

}

#' @export
median_na <- function(x) {
  
  if (is.nan(median(x, na.rm = T))) {
    return(NA)
  } else {
    return(median(x, na.rm = T))
  }
  
}

#' @export
sum_na <- function(x) {
  
  if (is.nan(sum(x, na.rm = T))) {
    return(NA)
  } else {
    return(sum(x, na.rm = T))
  }
  
}

#' @export
summarize_data <- function(input_dt) {
  
  message(
    "This OpenSWATH result contains ", length(unique(input_dt$PeptideIon)), " peptide ions; ", length(unique(input_dt$ProteinName)), " proteins."
  )
  message("Among them, ", length(which(complete.cases(input_dt))), " are completed rows (no NAs).")
  
}
