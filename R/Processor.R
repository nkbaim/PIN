#' @importFrom data.table dcast copy melt
#' @importFrom outliers chisq.out.test


#' @export
keep_only_proteotypic <- function(input_dt) {

  return(input_dt[which(grepl("^1/", input_dt$ProteinName)), ])

}


#' @export
long2wide <- function(input_dt, level="PeptideIon") {

  if (!level %in% c("PeptideIon","Protein")) {
    stop("Please select a valid type. Options:  \"PeptideIon(default)\", \"Protein\"")
  }

  if(level=="PeptideIon") {
    output_dt <- dcast(input_dt, ProteinName + PeptideIon ~ filename, value.var = "Intensity")
  } else if(level=="Protein") {
    output_dt <- dcast(input_dt, ProteinName ~ run_id, value.var="Quant_sum")
  }

  return(output_dt)

}


#' @export
wide2long <- function(input_dt) {

  output_dt <- melt(input_dt, id.vars=c("PeptideIon", "ProteinName", "FullName", "NTT", "NMC"), 
                    measure.vars=names(input_dt)[which(grepl("mean_", names(input_dt)))], 
                    variable.name="run_id", 
                    value.name="Intensity", na.rm=T)

  return(output_dt)

}


#' @export
keep_only_semi <- function(input_dt) {

  semi_cons_peptides_long <- copy(input_dt)
  semi_cons_peptides_long[which(semi_cons_peptides_long$NTT == 2), ]$Intensity <- 0

  return(semi_cons_peptides_long)

}


#' @export
pept2prot <- function(input_pept_long) {

  output_prot_long <- input_pept_long[, .(Quant_sum = sum_na(Intensity), NumPeptides = length(which(Intensity>0)) ), by=.(ProteinName, run_id)]

  return(output_prot_long)

}


#' @export
normalize_data <- function(input_dt, input_index, normalization="mediancenter") {

  output_dt <- copy(input_dt)
  
  output_dt[, (input_index) := log2(output_dt[, input_index, with=F]) ]
    
  run_median <- sapply(output_dt[, input_index, with=F], median_na)
    
  for(i in 1:length(input_index)) {
    output_dt[, names(output_dt)[input_index[i]]] <- output_dt[, input_index[i], with=F] - run_median[i] + median(run_median)
    output_dt[, names(output_dt)[input_index[i]]] <- 2^output_dt[, input_index[i], with=F]
  }

  return(output_dt)

}


#' @export
calculate_stats <- function(input_dt, input_index) {

  index_semi <- which(input_dt$NTT<2)
  index_missed <- which(input_dt$NMC>0)
  index_normal <- which(input_dt$NTT==2 & input_dt$NMC==0)

  numPept <- apply(input_dt[, input_index, with=F], 2, function(x) length(which(x>0)))

  numSemi <- apply(input_dt[index_semi, input_index, with=F], 2, function(x) length(which(x>0)))
  numMissed <- apply(input_dt[index_missed, input_index, with=F], 2, function(x) length(which(x>0)))
  numNormal <- apply(input_dt[index_normal, input_index, with=F], 2, function(x) length(which(x>0)))

  TIC <- apply(input_dt[, input_index, with=F], 2, sum_na )

  proteomic_measure <- data.frame(cbind( names(input_dt)[input_index], numPept, numSemi, numMissed, numNormal, TIC ))
  names(proteomic_measure) <- c("filename", "numPept", "numSemi", "numMissed", "numNormal", "TIC")

  return(proteomic_measure)

}


#' @export
merge_replicates <- function(input_dt, input_anno) {

  output_dt <- copy(input_dt)

  list_samples <- unique(input_anno$sample_id) 
  
  message("It starts to merge replicates...")
  
  for(i in 1:length(list_samples)) { 
    cat("Processing ", i, " sample: ", list_samples[i], "\n")
    output_dt[, paste0("mean_", list_samples[i]) := apply(.SD, 1, mean_na), .SDcols=input_anno[input_anno$sample_id %in% list_samples[i], ]$filename]
  }
  
  message("Done with merging replicates...")

  return(output_dt)

}


#' @export
generate_iPIS_matrix <- function(input_all_dt, input_semi_dt) {

  index_temp <- which(grepl("^mean", names(input_all_dt)))

  m_iPIS <- copy(input_all_dt)

#for(i in 1:length(index_temp)) {
#m_iPIS[, index_temp[i], with=F] <- 1 - input_semi_dt[, index_temp[i], with=F] / input_all_dt[, index_temp[i], with=F]
#}

  m_iPIS[, index_temp] <- 1 - input_semi_dt[, index_temp, with=F] / input_all_dt[, index_temp, with=F]

  names(m_iPIS)[index_temp] <- sapply(strsplit(names(m_iPIS)[index_temp], "mean_"), "[[", 2)

  m_iPIS <- m_iPIS[order(apply(m_iPIS[, index_temp, with=F], 1, function(x) length(which(is.na(x))))), ]

  return(m_iPIS)

}


#' @export
calculate_PIN <- function(input_iPIS) {

  index_temp <- seq(2, dim(input_iPIS)[2], 1)

  m_PIN <- data.frame( sample_id = names(input_iPIS)[index_temp], 
                            PIN = apply(input_iPIS[,index_temp, with=F], 2, mean_na), 
                            pval = 1.0 )

  m_PIN <- m_PIN[order(m_PIN$PIN), ]

  for(i in 1: (dim(m_PIN)[1]-1) ) {

    if(grepl("lowest", chisq.out.test(m_PIN$PIN[i:dim(m_PIN)[1]], opposite=F)$alternative)) {
      m_PIN[i, ]$pval <- chisq.out.test(m_PIN$PIN[i:dim(m_PIN)[1]], opposite=F)$p.value
    } else {
      m_PIN[i, ]$pval <- chisq.out.test(m_PIN$PIN[i:dim(m_PIN)[1]], opposite=T)$p.value
    }

    if(i > 1) {
      if(m_PIN[i, ]$pval < m_PIN[i-1, ]$pval) {
        m_PIN[i, ]$pval <- m_PIN[i-1, ]$pval + m_PIN[i, ]$pval * 0.05
      }
    }

  }

  return(m_PIN)

}

