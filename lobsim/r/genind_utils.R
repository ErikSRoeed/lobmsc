#---- GENIND_UTILS -------------------------------------------------------------
#>                                                                            <#
#>    Erik Sandertun Roeed                                                    <#
#>                                                                            <#
#>    Contents                                                                <#
#>    function: drop_genind_NAs()                                             <#
#>    function: get_genind_snp_sequences()                                    <#
#>    function: reorder_genind_pops()                                         <#
#>    function: overwrite_genind_pops_with_groups()                           <#
#>    function: match_pop_to_group()                                          <#
#>    function: write_genind_popsizes_to_table()                              <#
#>                                                                            <#
#------------------------------------------------------------------------------#

drop_genind_NAs <- function(genind, drop = "individuals")
{
  cells_are_na <- genind2df(genind, usepop = FALSE) |> is.na()
  
  if (drop == "individuals")
  {
    ACROSS_ROWS <- 1
    individuals_have_na <- cells_are_na |> apply(
      ACROSS_ROWS, function(loci_are_na_for_individual)
      {
        individual_has_na <- any(loci_are_na_for_individual)
        return(individual_has_na)
      }
    )
    
    drop_count <- sum(individuals_have_na)
    warning(paste("Dropping", drop_count, "individuals."))
    
    individuals_to_keep <- which(! individuals_have_na)
    filtered_genind <- genind[individuals_to_keep]
    
    return(filtered_genind)
  }
  
  if (drop == "loci")
  {
    ACORSS_COLUMNS <- 2
    loci_have_na <- cells_are_na |> apply(
      ACROSS_COLUMNS, function(individuals_are_na_for_locus)
      {
        locus_has_na <- any(individuals_are_na_for_locus)
        return(locus_has_na)
      }
    )
    
    drop_count <- sum(loci_have_na)
    warning(paste("Dropping", drop_count, "loci."))
    
    loci_to_keep <- which(! loci_have_na)
    filtered_genind <- genind[loc = loci_to_keep]
    
    return(filtered_genind)
  }
  
  warning("None dropped as 'drop' != 'individuals' or 'loci'.")
}



get_genind_snp_sequences <- function(genind, drop_na_by = "individuals")
{
  ploidy <- unique(genind$ploidy)
  all_individuals_share_ploidy <- length(ploidy) == 1
  stopifnot(all_individuals_share_ploidy)
  
  genind <- genind |> drop_genind_NAs(drop = drop_na_by)
  
  haploid_genotypes <- genind |> genind2df(usepop = FALSE, oneColPerAll = TRUE)
  haplotype_numbers <- seq_len(ploidy)
  
  snp_sequences <- matrix(
    dimnames = list(rownames(haploid_genotypes), haplotype_numbers),
    nrow = nrow(haploid_genotypes),
    ncol = ploidy
  )
  
  for (genome_copy in haplotype_numbers)
  {
    genome_copy_suffix_regex <- paste("\\.", genome_copy, sep = "")
    
    which_haploid_genotypes_in_genome_copy <- names(haploid_genotypes) |> 
      grep(pattern = genome_copy_suffix_regex)
    
    haploid_genotypes_of_genome_copy <- haploid_genotypes |>
      select(all_of(which_haploid_genotypes_in_genome_copy))
    
    ACROSS_ROWS <- 1
    haplotype_snp_sequences <- haploid_genotypes_of_genome_copy |>
      apply(MARGIN = ACROSS_ROWS, FUN = str_flatten)
    
    snp_sequences[, genome_copy] <- haplotype_snp_sequences
  }
  
  return(snp_sequences)
}



reorder_genind_pops <- function(genind, new_pop_order)
{
  for (pop in new_pop_order)
  {
    pop_inds <- poppr::popsub(genind, sublist = pop)
    if (! exists("ordered_genind"))
    {
      ordered_genind <- pop_inds
      next
    }
    ordered_genind <- repool(ordered_genind, pop_inds)
  }
  return(ordered_genind)
}



overwrite_genind_pops_with_groups <- function(genind, group_list)
{
  genind$pop <- genind$pop |>
    sapply(match_pop_to_group, group_list = group_list) |> 
    factor(levels = names(group_list))
  return(genind)
}



match_pop_to_group <- function(pop, group_list)
{
  pop_group <- sapply(group_list, function(group) pop %in% group) |> 
    which() |> 
    names()
  
  if (length(pop_group) == 0)
  {
    stop(paste("No group to assign pop", pop, "in group list!"))
  }
  
  return(pop_group)
}



write_genind_popsizes_to_table <- function(genind, output_dir)
{
  genind_name <- deparse(substitute(genind))
  output_filepath <- paste0(output_dir, "/", genind_name, ".csv")
  genind |>
    pop() |>
    table() |>
    write.table(file = output_filepath, sep = ";", row.names = FALSE)
}
