#' @importFrom data.table :=
#' @importFrom data.table as.data.table


#' @export

generate_lib_info <- function(input_sptxt = "cons_fixed.sptxt") {

FullName <- system(paste0("grep \"FullNa\" ", input_sptxt, "| awk '{print $2}'" ), intern = T, ignore.stderr = F)
NMC <- system(paste0(" grep \"^Comm\" ", input_sptxt, " | sed 's|.* NMC=\\([0-9]*\\) .*|\\1|g' " ), intern = T, ignore.stderr = F)
NTT <- system(paste0(" grep \"^Comm\" ", input_sptxt, " | sed 's|.* NTT=\\([0-9]*\\) .*|\\1|g' " ), intern = T, ignore.stderr = F)
Protein <- system(paste0(" grep \"^Comm\" ", input_sptxt, " | sed 's|.* Protein=\\([^[:space:]]*\\) .*|\\1|g' " ), intern = T, ignore.stderr = F)

raw <- as.data.table(cbind(FullName, NTT, NMC, Protein))

data <- unique(raw)
names(data) <- c("FullName", "NTT", "NMC", "Prot")

data[, PeptideIon := paste0(sapply(strsplit(data$FullName, "\\."), "[[", 2), "_", sapply(strsplit(data$FullName, "\\/"), "[[", 2)) ]

write.table(data, file="lib_peptide.tsv", col.names=T, row.names=F, sep="\t", quote=F)

return(data)

}
