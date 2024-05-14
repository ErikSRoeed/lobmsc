#---- VCF_UTILS ----------------------------------------------------------------
#>                                                                            <#
#>    Erik Sandertun Roeed                                                    <#
#>                                                                            <#
#>    Contents                                                                <#
#>    function: import_slim_vcf_as_genind()                                   <#
#>    function: export_genind_as_vcf()                                        <#
#>    function: collect_vcfs_as_pooled_genind()                               <#
#>                                                                            <#
#------------------------------------------------------------------------------#

import_slim_vcf_as_genind <- function(
    vcf_filepath,
    reference_snp_sequence,
    snp_positions,
    position_ids
)
{
  reference_alleles <- reference_snp_sequence |> stringr::str_split_1("")
  
  vcf <- read.vcfR(vcf_filepath) # Loci with ancestral reference missing!
  vcf_n_loci <- nrow(vcf@fix)
  vcf_locus_positions <- vcf@fix[seq_len(vcf_n_loci), "POS"]
  
  new_fixed_matrix <- matrix(
    dimnames = dimnames(vcf@fix),
    ncol = ncol(vcf@fix),
    nrow = length(snp_positions)
  )
  
  new_genotype_matrix <- matrix(
    dimnames = dimnames(vcf@gt),
    ncol = ncol(vcf@gt),
    nrow = length(snp_positions)
  )
  
  new_fixed_reference_row <- setNames(
    object = c("1", rep(NA, 7)),
    nm = colnames(vcf@fix)
  )
  
  new_genotype_reference_row <- setNames(
    object = c("GT", rep("0|0", ncol(vcf@gt) - 1)),
    nm = colnames(vcf@gt)
  )
  
  for (position in snp_positions)
  {
    if (position %in% vcf_locus_positions)
    {
      which_vcf_row <- which(vcf@fix[, "POS"] == position)
      new_fixed_matrix[position, ] <- vcf@fix[which_vcf_row, ]
      new_genotype_matrix[position, ] <- vcf@gt[which_vcf_row, ]
      
      next
    }
    
    new_fixed_reference_row["POS"] <- position
    new_fixed_reference_row["REF"] <- reference_alleles[position]
    
    new_fixed_matrix[position, ] <- new_fixed_reference_row
    new_genotype_matrix[position, ] <- new_genotype_reference_row
  }
  
  new_fixed_matrix[, "ID"] <- position_ids
  
  vcf_including_ancestral_loci <- vcf
  vcf_including_ancestral_loci@fix <- new_fixed_matrix
  vcf_including_ancestral_loci@gt <- new_genotype_matrix
  
  vcf_as_genind <- vcf_including_ancestral_loci |> 
    vcfR2genind(return.alleles = TRUE)
  
  return(vcf_as_genind)
}



export_genind_as_vcf <- function(
    genind,
    drop_na_by = "individuals",
    reference_snp_sequence,
    snp_positions,
    vcf_filepath
)
{
  VCF_META_HEADERS <- c(
    "##fileformat=VCF"
  )
  VCF_COLUMN_NAMES <- c(
    "#CHROM",
    "POS",
    "ID",
    "REF",
    "ALT",
    "QUAL",
    "FILTER",
    "INFO",
    "FORMAT"
  )
  BLANK_VALUE <- "."
  DEFAULT_CHROM <- "1"
  FORMAT_VALUES <- "GT"
  REF_ALLELE_NAME <- "0"
  
  genind <- genind |> drop_genind_NAs(drop = drop_na_by)
  
  reference_alleles <- reference_snp_sequence |> stringr::str_split_1("")
  
  number_of_cols <- length(VCF_COLUMN_NAMES)
  number_of_loci <- length(genind$loc.n.all)
  locus_numbers <- seq_len(number_of_loci)
  
  alleles_at_each_locus <- lapply(
    locus_numbers,
    function(locus_number)
    {
      ref_allele_at_locus <- reference_alleles[locus_number]
      all_alleles_at_locus <- genind$all.names[[locus_number]]
      which_allele_is_ref <- which(all_alleles_at_locus == ref_allele_at_locus)
      
      alt_count <- length(all_alleles_at_locus) - length(which_allele_is_ref)
      alt_allele_names <- seq_len(alt_count)|> as.character()
      
      all_allele_names_at_locus <- alt_allele_names |> 
        append(REF_ALLELE_NAME, after = which_allele_is_ref - 1)
      
      named_alleles_at_locus <- setNames(
        object = all_alleles_at_locus,
        nm = all_allele_names_at_locus
      )
      
      return(named_alleles_at_locus)
    }
  )
  
  vcf_template_matrix <- matrix(
    nrow = number_of_loci,
    ncol = number_of_cols,
    dimnames = list(
      locus_numbers,
      VCF_COLUMN_NAMES
    )
  )
  
  vcf_table <- as.data.frame(vcf_template_matrix) |> 
    mutate(`#CHROM` = rep(DEFAULT_CHROM, number_of_loci)) |>
    mutate(POS = snp_positions) |> 
    mutate(ID = names(genind$all.names)) |> 
    mutate(REF = reference_alleles) |> 
    mutate(
      ALT = lapply(
        alleles_at_each_locus,
        function(locus_alleles)
        {
          which_allele_is_ref <- which(names(locus_alleles) == REF_ALLELE_NAME)
          alt_alleles = locus_alleles[-which_allele_is_ref]
          comma_separated_alt_alleles = paste(alt_alleles, collapse = ",")
          return(comma_separated_alt_alleles)
        }
      ) |> unlist()
    ) |> 
    mutate(QUAL = rep(BLANK_VALUE, number_of_loci)) |> 
    mutate(FILTER = rep(BLANK_VALUE, number_of_loci)) |> 
    mutate(INFO = rep(BLANK_VALUE, number_of_loci)) |> 
    mutate(FORMAT = rep(FORMAT_VALUES, number_of_loci))
  
  genotypes_to_export <- adegenet::genind2df(
    x = genind,
    usepop = FALSE,
    sep = "|"
  )
  
  for (locus_number in locus_numbers)
  {
    locus_allele_numbers <- setNames(
      object = names(alleles_at_each_locus[[locus_number]]),
      nm = alleles_at_each_locus[[locus_number]]
    )
    
    locus_genotypes <- genotypes_to_export[, locus_number]
    
    number_genotypes_at_locus <- sapply(
      locus_genotypes,
      function(individual_genotype_at_locus)
      {
        if(is.na(individual_genotype_at_locus))
        {
          return(NA)
        }
        
        individual_number_genotype_at_locus <- individual_genotype_at_locus |>
          stringr::str_replace_all(pattern = locus_allele_numbers)
        
        return(individual_number_genotype_at_locus)
      }
    )
    
    genotypes_to_export[, locus_number] <- number_genotypes_at_locus
  }
  
  genotypes_to_export_with_loci_as_rows <- t(genotypes_to_export)
  vcf_table <- cbind(vcf_table, genotypes_to_export_with_loci_as_rows)
  
  write.table(
    x = vcf_table,
    file = vcf_filepath,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
  vcf_contents_without_header <- readLines(vcf_filepath)
  with_header <- c(VCF_META_HEADERS, vcf_contents_without_header)
  writeLines(text = with_header, con = vcf_filepath)
}



collect_vcfs_as_pooled_genind <- function(
    vcf_filepaths,
    reference_snp_sequence,
    snp_positions,
    position_ids
)
{
  pooled_genind <- vcf_filepaths |> 
    lapply(
      function(path)
      {
        import_slim_vcf_as_genind(
          vcf_filepath = path,
          reference_snp_sequence = reference_snp_sequence,
          snp_positions = snp_positions,
          position_ids = position_ids
        )
      }
    ) |> 
    repool()
  
  return(pooled_genind)
}
