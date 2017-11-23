

#' @export
write_tsv <- function(input_table) {
  
  write.table(
    input_table,
    file = paste0(deparse(substitute(input_table)), ".tsv"),
    col.name = T,
    row.name = F,
    sep = "\t",
    quote = F
  )
  
}
