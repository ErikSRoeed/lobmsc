#---- lobsim ------------------------------------------------------------------#
#                                                                              #
#     || LOBSTER SIMULATIONS                                                   #
#     || Project main script                                                   #
#                                                                              #
#     Erik Sandertun Røed                                                      #
#                                                                              #
#     Project includes data provided by, and reproduced here with permission   #
#     of, Dr. Charlie Ellis (University of Exeter)                             #
#                                                                              #
#------------------------------------------------------------------------------#

#---- 0 SETUP: USER INPUT ------------------------------------------------------
#>    Set user-defined values                                                  
#------------------------------------------------------------------------------#

# Set these as per your installation of SLiM and available CPU nodes
YOUR_SLIM_PATH <- "C:/msys64/mingw64/bin/slim.exe" # Default for Windows/MINGW64
PARALLEL_NODES <- 10 # Expect up to around 1GB RAM used per node

# Lobster for reference alleles (arbitrarily Sin55 from Singlefjord, Norway)
REFERENCE_LOBSTER <- "Sin55" 

# Locales/pops grouped as per Ellis et al. (2020)
GENOTYPE_GROUPS <- list(
  # Atlantic locales/pops (broadly defined)
  `H_gammarus` = c(
    "Tro", "Ber", "Flo", "Sin", "Gul", "Kav", "Lys", "Hel",
    "Oos",
    "Cro", "Brd", "Eye", "She", "Ork", "Heb", "Sul", "Arr", "Don", "Iom",
    "Cor", "Hoo", "Mul", "Kil", "Ven", "Lyn", "Pem", "Pad", "EuroCook",
    "Ios", "Loo", "Jer", "Sbs", "Bre", "Idr", "Vig", "Pen", "Far", "Tan"
  ),
  # Mediterranean locales/pops
  `Mediterrenean_gammarus` = c(
    "Csa", "Sar", "Laz", "Ion",
    "Adr", "Ale", "Sky", "Spo", "Chi", "The", "Tor"
  ),
  # Hybrid and americanus
  `Hybrid` = c("HybridX"),
  `H_americanus` = c("Americanus", "AmerCook")
)

# Remove Mediterranean samples for simplicity
GENOTYPE_GROUPS_KEPT <- GENOTYPE_GROUPS[-2]

#---- 1 SETUP: HOUSEKEEPING ----------------------------------------------------
#>    Source all scripts and packages, and check installation of slimmr        
#------------------------------------------------------------------------------#

# Echo for background jobs
cat("Step 1: Housekeeping")

SOURCE_FILES <- dir("r", full.names = TRUE)
DEPENDENCIES <- c(
  # For simulations (devtools::install_github can be used to install slimmr)
  "slimmr", # devtools,
  # Core dependencies for genotype storage and wrangling, I/O, analyses
  "adegenet", "poppr", "vcfR",
  # For data wrangling
  "dplyr", "tidyr", "purrr", "stringr", "forcats", "magrittr",
  # For plotting
  "ggplot2", "patchwork", "viridis", "RColorBrewer"
)

sapply(SOURCE_FILES, source)
stopifnot(check_slimmr_installed())
sapply(DEPENDENCIES, library, character.only = TRUE)

#---- 2 SETUP: PARSE EMPIRICAL EXETER DATA -------------------------------------
#>    Import data, get names of SNPs, generate SLiM SNP positions, extract
#>    SNPs in a randomly chosen genome copy from reference individual as an
#>    ancestral "sequence", and export a VCF file for each genotype
#>    group defined in the user input above.
#>    
#>    Also test DAPC and Snapclust on raw data with and without NAs.
#------------------------------------------------------------------------------#

# Echo for background jobs
cat("Step 2: Parsing empirical data")

# Where to obtain input and place output
EXETER_DATA_RDS <- "exeter/data/exeter_lobster_1591ind_79snps_52pop.rds"
EXETER_VCF_DIR <- "output/00_sourcedata_as_vcf_files/vcf"
dir.create(EXETER_VCF_DIR, recursive = TRUE)

# Geninds including NA values (but excluding Mediterranean)
lobsters_with_NA <- readRDS(EXETER_DATA_RDS) |> 
  poppr::popsub(unlist(GENOTYPE_GROUPS_KEPT)) |> 
  reorder_genind_pops(GENOTYPE_GROUPS_KEPT)
  
# Geninds excluding NA values
lobsters <- lobsters_with_NA |> drop_genind_NAs()

# Extract info on SNPs
original_snp_ids <- locNames(lobsters)
snp_positions_in_slim <- seq_len(nLoc(lobsters))
reference_genotype_as_sequence <- lobsters[REFERENCE_LOBSTER] |> 
  get_genind_snp_sequences() |>
  sample(1)

# Export all genind data from kept genotype groups to VCF for import in SLiM
for (pop_name in unlist(GENOTYPE_GROUPS_KEPT))
{
  export_genind_as_vcf(
    genind = lobsters[pop = pop_name],
    snp_positions = snp_positions_in_slim,
    reference_snp_sequence = reference_genotype_as_sequence,
    vcf_filepath = paste0(EXETER_VCF_DIR, "/", pop_name, ".vcf")
  )
}

# Run DAPC and Snapclust on empirical data, plots show effect of excluding NAs
empirical_data_dapcs <- list(
  with_NA = apply_exeter_dapc_pipeline(lobsters_with_NA),
  without_NA = apply_exeter_dapc_pipeline(lobsters)
)

empirical_data_snapclusts <- list(
  with_NA = apply_exeter_snapclust_pipeline(lobsters_with_NA),
  without_NA = apply_exeter_snapclust_pipeline(lobsters)
)

empirical_data_plots <- (
  guide_area() /
    (
      (
        (
          custom_dapc_plot(empirical_data_dapcs[[1]]) +
            custom_snapclust_plot(empirical_data_snapclusts[[1]])
        ) + plot_layout(widths = c(9, 1))
      ) | (
        (
          custom_dapc_plot(empirical_data_dapcs[[2]]) +
            custom_snapclust_plot(empirical_data_snapclusts[[2]])
        ) + plot_layout(widths = c(9, 1))
      ) 
    ) +
    plot_annotation(tag_level = "a") +
    plot_layout(heights = c(0.5, 4.5), guides = "collect")
)

empirical_data_plots |> 
  ggsave(
    filename = "output/00_sourcedata_as_vcf_files/00_figure.png",
    width = 10,
    height = 5,
    dpi = 600
  )

#---- 3 SETUP: IMPORT MAIN SLiM MODEL ------------------------------------------
#>    This is the model that will branch out with slimmr. Remove SLiMGUI
#>    code blocks that are used for prototyping or one-off runs in SLiMGUI.
#------------------------------------------------------------------------------#

# For background jobs
cat("Step 3: Import main SLiM model")

# Import
slim_model_main <- slimmr::import_slim_model(
  script_path = "slim/main_model.slim",
  name = "Main"
)

# Remove SLiMGUI blocks
slim_model_main |> slimmr::remove_blocks(
  c(
    1, # SLiMGUI initialize() callback
    6, # SLiMGUI modifyChild() callback for inheriting SLiMGUI colours
    9  # SLiMGUI 1 early() callback for setting SLiMGUI colours
  )
)

#---- 4 VALIDATION: TEST BASIC INPUT AND OUTPUT FROM SLiM ----------------------
#>    For each (Atlantic) gammarus, americanus, and hybrid pop in the empirical
#>    data, simulate analogous "simulobster" population of size N. Compare
#>    population assignment between simulobster and empirical populations to
#>    verify that the SLiM model works as expected.
#------------------------------------------------------------------------------#

# For background jobs
cat("Step 4: Test basic input and output from SLiM")

# Set up directories for output
TEST_DIR <- "output/01_test_slim_input_output/"
TEST_VCF_DIR <- "output/01_test_slim_input_output/vcf"
TEST_RDS_DIR <- "output/01_test_slim_input_output/rds"
dir.create(TEST_VCF_DIR, recursive = TRUE)
dir.create(TEST_RDS_DIR, recursive = TRUE)

# Branch and manipulate main SLiM model (slim_model_main)
slim_model_test <- slimmr::branch_model(
  slim_model = slim_model_main,
  branch_name = "Test I/O"
)

# No need for multiple populations, parentage tracking, or related outputs
slim_model_test |> slimmr::remove_blocks(c(3, 4, 6, 7, 8))

# Set last generation to generation 100 for output and simulation end blocks
slim_model_test |> slimmr::change_block_callback(
  block_index = 4, # 125 late()
  new_callback = "100 late()"
)

slim_model_test |> slimmr::change_block_callback(
  block_index = 5, # 125 late()
  new_callback = "100 late()"
)

# Newly simulated simulobsters will be bound to this instead of "lobsters":
lobsters_and_simulobsters <- lobsters

# BEGIN SIMULATIONS LOOP ACROSS INCLUDED EMPIRICAL POPS
for (pop_name in rev(levels(lobsters$pop))) # rev() for plotting convenience
{
  pop_real_lobsters <- poppr::popsub(
    lobsters_and_simulobsters,
    sublist = pop_name
  )
  input_pop_size <- nInd(pop_real_lobsters)
  
  # The SLiM script outputs VCF files
  slim_model_test |> slimmr::run_slim(
    slim_command = YOUR_SLIM_PATH,
    parallel_nodes = 2,
    reps = 2,
    # SLiM constant definitions below
    WORKING_DIR = paste0("'", getwd(), "'"),
    OUTPUT_DIR = paste0("'", TEST_VCF_DIR, "'"),
    OUTPUT_VCF_SUFFIX = paste0("'_", pop_name, "_SIM'"),
    REFERENCE_SNPS = paste0("'", reference_genotype_as_sequence, "'"),
    INPUT_VCF = paste0("'", EXETER_VCF_DIR, "/", pop_name, ".vcf'"),
    INPUT_N = paste0('"', input_pop_size, '"'),
    N = "250",
    NU = "79",
    MU = "0.0",
    RHO = "0.5"
  )
  
  # Collect the VCF file paths
  simulated_vcfs <- dir(
    TEST_VCF_DIR,
    pattern = pop_name,
    full.names = TRUE
  )
  
  # Import VCFs as a pooled genind
  simulobsters <- simulated_vcfs |> 
    collect_vcfs_as_pooled_genind(
      reference_snp_sequence = reference_genotype_as_sequence,
      snp_positions = snp_positions_in_slim,
      position_ids = original_snp_ids
    )
  
  simulobster_pop_name <- stringr::str_c(pop_name, "_SIM")
  
  simulobsters$pop <- factor(
    rep(simulobster_pop_name, nInd(simulobsters)),
    levels = simulobster_pop_name
  )
  
  # Bind to lobsters_and_simulobsters
  lobsters_and_simulobsters <- lobsters_and_simulobsters |> repool(simulobsters)
  
  # Finally save
  simulobsters |> saveRDS(
    file = paste0(
      TEST_RDS_DIR,
      "/",
      simulobster_pop_name, ".rds"
    )
  )
}
# END OF SIMULATIONS LOOP

# Update the genotype group list with simulated individuals (interleaved)
genotype_groups_with_simulobsters <- GENOTYPE_GROUPS_KEPT |> 
  append(
    list(
      `H_gammarus*` = unlist(GENOTYPE_GROUPS_KEPT[1]) |> paste0("_SIM")
    ),
    after = 0
  ) |>
  append(
    list(
      `Hybrid*` = unlist(GENOTYPE_GROUPS_KEPT[2]) |> paste0("_SIM")
    ),
    after = 2
  ) |> 
  append(
    list(
      `H_americanus*` = unlist(GENOTYPE_GROUPS_KEPT[3]) |> paste0("_SIM")
    ),
    after = 4
  )

# Colour palettes for plotting
colour_fills <- RColorBrewer::brewer.pal(n = 6, "Paired") |> rev()

# Match the "Paired" colours up with "Set1" colours (used earlier)
matched_fills <- colour_fills[c(2,1,6,5,4,3)]

# Point outline colours
outline_colours <- c(
  matched_fills[1], "gray10",
  matched_fills[3], "gray10",
  matched_fills[5], "gray10"
)

# Run analyses
test_dapc <- apply_exeter_dapc_pipeline(
  genind_data = lobsters_and_simulobsters,
  group_list = genotype_groups_with_simulobsters
)

test_snapclust <- apply_exeter_snapclust_pipeline(
  genind_data = lobsters_and_simulobsters,
  true_group_list = genotype_groups_with_simulobsters
)

test_plot <- (
  custom_dapc_plot(
    test_dapc,
    fill_colours = matched_fills,
    outline_colours = outline_colours
  ) |
    custom_snapclust_plot(
      test_snapclust,
      true_group_labels = c(
        "H_gammarus*" = "H gammarus*",
        "H_gammarus" = "H gammarus",
        "Hybrid*" = "Hybrid*",
        "Hybrid" = "Hybrid",
        "H_americanus*" = "H americanus*",
        "H_americanus" = "H americanus"
      )
    )) +
  patchwork::plot_layout(widths = c(9, 1)) +
  patchwork::plot_annotation(tag_levels = "a")

test_plot |> ggsave(
  filename = paste0(TEST_DIR, "/01_figure.png"),
  width = 10,
  height = 7.5,
  dpi = 600
)

#---- 5 SIMULATION: RUN SIMULATED CROSSING EXPERIMENT --------------------------
#>    Run a simulated crossing experiment with reciprocal migration between
#>    one population initated from gammarus data and one initiated from 
#>    americanus data. Save the data for comparisons and analyses below.
#------------------------------------------------------------------------------#

# For background jobs
cat("Step 5: Run simulated crossing experiment")

# Choice of empirical data to initialise simulations
CROSS_AMERICANUS_POP <- "Americanus" # Live Americanus (not AmerCook)
CROSS_GAMMARUS_POP   <- "Sin"        # Singlefjord, Norway

# Directories for output
SIM_CROSS_VCF_DIR <- "output/02_simulated_cross/vcf"
SIM_CROSS_RDS_DIR <- "output/02_simulated_cross/rds"
dir.create(SIM_CROSS_VCF_DIR, recursive = TRUE)
dir.create(SIM_CROSS_RDS_DIR, recursive = TRUE)

# Branch main model
slim_model_crossing <- slimmr::branch_model(
  slim_model = slim_model_main,
  branch_name = "Simulate crossing"
)

# Change last generation to 105 (to retain some purebred genotypes)
slim_model_crossing |> slimmr::change_block_callback(10, "105 late()")
slim_model_crossing |> slimmr::change_block_callback(8, "100:105 late()")

# Add reciprocal migration (from gammarus to americanus)
slim_model_crossing |> slimmr::add_lines(
  lines = "\tp2.setMigrationRates(sourceSubpops = p1, rates = M);",
  after_line = 52
)

# Remove regular bulk output of all individuals
slim_model_crossing |> slimmr::remove_blocks(9) # 125 late()

# Replace with targeted output of only pure gammarus and americanus at end
slim_model_crossing |> slimmr::add_blocks(
  blocks_script = c(
    "105 late() {",
    "\toutputHybrids(0.0, '_simcross_0.000');",
    "\toutputHybrids(1.0, '_simcross_1.000');",
    "}",
    ""
  ),
  after_block = 8
)

# Hardcode output of F1 individuals
slim_model_crossing |> slimmr::replace_eidos_pattern(
  pattern = "OUTPUT2_VCF_SUFFIX",
  replacement = "'_simcross_0.500'",
  in_lines = 57
)

# Add output of F2/3
slim_model_crossing |> slimmr::add_lines(
  lines = c(
    "\toutputHybrids(0.125, '_simcross_0.125');",
    "\toutputHybrids(0.25, '_simcross_0.250');",
    "\toutputHybrids(0.75, '_simcross_0.750');",
    "\toutputHybrids(0.875, '_simcross_0.875');"
  ),
  after_line = 57
)

# Collect information from choice of empirical data
cross_americanus_lobsters <- lobsters |> poppr::popsub(CROSS_AMERICANUS_POP)
cross_gammarus_lobsters <- lobsters |> poppr::popsub(CROSS_GAMMARUS_POP)
n_americanus <- nInd(cross_americanus_lobsters)
n_gammarus <- nInd(cross_gammarus_lobsters)

# Run simulations to simulate backcrosses; as above, but output F1-F3 hybrids
slim_model_crossing |> slimmr::run_slim(
  slim_command = YOUR_SLIM_PATH,
  parallel_nodes = PARALLEL_NODES,
  reps = 10,
  # SLiM constant definitions below
  WORKING_DIR = paste0("'", getwd(), "'"),
  OUTPUT_DIR = paste0("'", SIM_CROSS_VCF_DIR, "'"),
  REFERENCE_SNPS = paste0("'", reference_genotype_as_sequence, "'"),
  INPUT_VCF = paste0("'", EXETER_VCF_DIR, "/", CROSS_GAMMARUS_POP,".vcf'"),
  INPUT2_VCF = paste0("'", EXETER_VCF_DIR, "/", CROSS_AMERICANUS_POP,".vcf'"),
  INPUT_N = as.character(n_gammarus),
  INPUT2_N = as.character(n_americanus),
  N = "1000",
  N2 = "1000",
  M = "0.01",
  NU = "79",
  MU = "0.0",
  RHO = "0.5"
)

# Collect VCF files per level of Americanus parentage (0 - 1)
# Escape mismatched files with \\ before . (which is a regex special symbol)
simcross_vcf_files <- list(
  `0.000` = dir(SIM_CROSS_VCF_DIR, pattern = "0\\.000", full.names = TRUE),
  `0.125` = dir(SIM_CROSS_VCF_DIR, pattern = "0\\.125", full.names = TRUE),
  `0.250` = dir(SIM_CROSS_VCF_DIR, pattern = "0\\.250", full.names = TRUE),
  `0.500` = dir(SIM_CROSS_VCF_DIR, pattern = "0\\.500", full.names = TRUE),
  `0.750` = dir(SIM_CROSS_VCF_DIR, pattern = "0\\.750", full.names = TRUE),
  `0.875` = dir(SIM_CROSS_VCF_DIR, pattern = "0\\.875", full.names = TRUE),
  `1.000` = dir(SIM_CROSS_VCF_DIR, pattern = "1\\.000", full.names = TRUE)
)

# Produce list of geninds (one per level of Americanus parentage)
simcross_geninds <- simcross_vcf_files |>
  lapply(
    collect_vcfs_as_pooled_genind,
    reference_snp_sequence = reference_genotype_as_sequence,
    snp_positions = snp_positions_in_slim,
    position_ids = original_snp_ids
  )

# Set genind population names to Americanus parentage level
simcross_geninds <- names(simcross_geninds) |> 
  lapply(
    function(americanus_fraction) {
      backcross_level_genind <- simcross_geninds[[americanus_fraction]]
      backcross_level_genind$pop <- rep(
        factor(americanus_fraction, levels = names(simcross_geninds)),
        nInd(backcross_level_genind)
      )
      return(backcross_level_genind)
    }
  )

# Save
simcross_geninds |> saveRDS(
  file = paste0(SIM_CROSS_RDS_DIR, "/simcross_geninds.rds")
)

#---- 6 ANALYSIS: COMPARE SIMULATED F1 HYBRID ASSIGNMENT TO EMPIRICAL DATA -----
#>    Compare the assignment of simulated F1 hybrids (empirical data-based) to
#>    F1 hybrids "bred in the simulation" from crosses of simulated americanus
#>    and gammarus data (i.e. from previous step).
#------------------------------------------------------------------------------#

# For background jobs
cat("Step 6: Check simulated hybrid assignment")

# Directories to save output
SIMCROSS_F1_PLOT_OUTPUT_DIR <- "output/03_check_simcross_f1_assignment"
dir.create(SIMCROSS_F1_PLOT_OUTPUT_DIR, recursive = TRUE)

# Obtain F1 lobsters from previous step
simcross_f1_simulobsters <- simcross_geninds[[4]]

# Set pop name for SIMcross lobsters
simcross_f1_simulobsters$pop <- factor(
  rep("simcross_F1", nInd(simcross_f1_simulobsters)),
  levels = "simcross_F1"
)

# Colour for SIMcross genotypes
simcross_fill <- RColorBrewer::brewer.pal(n = 7, "Paired")[7]
simcross_outline <- simcross_fill
plot_fills <- append(matched_fills, simcross_fill, after = 2)
plot_outlines <- append(outline_colours, simcross_outline, after = 2)

# Repeat analyses from validation step (4) but now with simcross F1 hybrids
simcross_f1_check_dapc <- apply_exeter_dapc_pipeline(
  genind_data = repool(simcross_f1_simulobsters, lobsters_and_simulobsters),
  group_list = append(
    x = genotype_groups_with_simulobsters,
    list(`Simcross F1*` = c("simcross_F1")),
    after = 2
  )
)

simcross_f1_check_snapclust <- apply_exeter_snapclust_pipeline(
  genind_data = repool(simcross_f1_simulobsters, lobsters_and_simulobsters),
  true_group_list = append(
    x = genotype_groups_with_simulobsters,
    list(`Simcross F1*` = c("simcross_F1")),
    after = 2
  )
)

simcross_f1_check_figure <- (
  custom_dapc_plot(
    simcross_f1_check_dapc,
    fill_colours = plot_fills,
    outline_colours = plot_outlines
  ) |
    custom_snapclust_plot(
      simcross_f1_check_snapclust,
      true_group_labels = c(
        "H_gammarus*" = "H gammarus*",
        "H_gammarus" = "H gammarus",
        "Hybrid*" = "Hybrid*",
        "Hybrid" = "Hybrid",
        "simcross_F1" = "Simcross F1*",
        "H_americanus*" = "H americanus*",
        "H_americanus" = "H americanus"
      )
    )) +
  patchwork::plot_layout(widths = c(9, 1)) +
  patchwork::plot_annotation(tag_levels = "a")

simcross_f1_check_figure |> ggsave(
  filename = paste0(SIMCROSS_F1_PLOT_OUTPUT_DIR, "/03_figure.png"),
  width = 10,
  height = 7.5,
  dpi = 600
)

#---- 7 ANALYSIS: EXAMPLE APPLICATION: F2/3 DETECTION --------------------------
#>    As an example of applying the SLiM simulations, test if the panel is
#>    likely to detect F2/F3 hybrids. A somewhat contrived example, given that
#>    you could probably do the same with adegenet::hybridise(), but this
#>    solution implemented in SLiM is much more flexible and can be extended.
#------------------------------------------------------------------------------#

# For background jobs
cat("Step 7: Example application - can the panel detect backcrosses?")

# Create directory for output
BACKCROSS_DETECTION_DIR <- "output/04_backcross_detection"
dir.create(BACKCROSS_DETECTION_DIR, recursive = TRUE)

# Create combined geninds for analyses
simcross_geninds_f1f2 <- repool(simcross_geninds[-c(2,6)])
simcross_geninds_f1f2f3 <- repool(simcross_geninds)

# Group list for ordering populations in analyses / plots
simcross_genotype_group_list <- list(
  `0`     = "0.000",
  `0.125` = "0.125",
  `0.25`  = "0.250",
  `0.5`   = "0.500",
  `0.75`  = "0.750",
  `0.875` = "0.875",
  `1`     = "1.000"
)

snapclust_f1f2 <- apply_exeter_snapclust_pipeline(
  genind_data = simcross_geninds_f1f2,
  hybrid_coef = c(0.25, 0.5),
  true_group_list = simcross_genotype_group_list
)

snapclust_f1f2f3 <- apply_exeter_snapclust_pipeline(
  genind_data = simcross_geninds_f1f2f3,
  hybrid_coef = c(0.125, 0.25, 0.5) ,
  true_group_list = simcross_genotype_group_list
)

# Plot panels
p1 <- custom_snapclust_assignment_dot_plot(
  snapclust_f1f2,
  actual_group_title = "Actual H. americanus* parentage",
  assigned_group_title = "Assigned H. americanus* parentage"
)
p2 <- custom_snapclust_plot(
  snapclust_f1f2,
  suppress_ticks = TRUE,
  actual_group_label = "Actual H. americanus* parentage"
)
p3 <- custom_snapclust_assignment_dot_plot(
  snapclust_f1f2f3,
  actual_group_title = "Actual H. americanus* parentage",
  assigned_group_title = "Assigned H. americanus* parentage"
)
p4 <- custom_snapclust_plot(
  snapclust_f1f2f3,
  suppress_ticks = TRUE,
  actual_group_label = "Actual H. americanus* parentage"
)

f1f2_backcross_assignment_plot <- (p1 | p2) +
  plot_layout(widths = c(4.25, 0.75))

f1f2f3_backcross_assignment_plot <- (p3 | p4) +
  plot_layout(widths = c(4.25, 0.75))

combined_backcross_assignment_plot <- (p1 | p2 | p3 | p4) +
  plot_layout(widths = c(4.25, 0.75, 4.25, 0.75)) +
  plot_annotation(tag_levels = "a")

f1f2_backcross_assignment_plot |> ggsave(
  filename = paste0(BACKCROSS_DETECTION_DIR, "/04_figure1.png"),
  width = 6,
  height = 6,
  dpi = 600
)

f1f2f3_backcross_assignment_plot |> ggsave(
  filename = paste0(BACKCROSS_DETECTION_DIR, "/04_figure2.png"),
  width = 6,
  height = 6,
  dpi = 600
)

combined_backcross_assignment_plot |> ggsave(
  filename = paste0(BACKCROSS_DETECTION_DIR, "/04_figure3.png"),
  width = 12.5,
  height = 7.5,
  dpi = 600
)

#---- 7 UTILITY: COUNT (SIMU)LOBSTERS IN VARIOUS GROUPS ------------------------
#>    For ease of writing about these results, output the number of lobsters in
#>    each group, across all steps!
#------------------------------------------------------------------------------#

# Create directory for output
COUNTS_DIR <- "output/05_lobster_simulobster_counts"
dir.create(COUNTS_DIR, recursive = TRUE)

grouped_lobsters_with_NA <- lobsters_with_NA |>
  overwrite_genind_pops_with_groups(GENOTYPE_GROUPS_KEPT)
grouped_lobsters <- lobsters |>
  overwrite_genind_pops_with_groups(GENOTYPE_GROUPS_KEPT)
grouped_lobsters_and_simulobsters <- lobsters_and_simulobsters |> 
  overwrite_genind_pops_with_groups(genotype_groups_with_simulobsters)

lobsters_with_NA |> write_genind_popsizes_to_table(COUNTS_DIR)
grouped_lobsters_with_NA |> write_genind_popsizes_to_table(COUNTS_DIR)
lobsters |> write_genind_popsizes_to_table(COUNTS_DIR)
grouped_lobsters |> write_genind_popsizes_to_table(COUNTS_DIR)
lobsters_and_simulobsters |> write_genind_popsizes_to_table(COUNTS_DIR)
grouped_lobsters_and_simulobsters |> write_genind_popsizes_to_table(COUNTS_DIR)
simcross_geninds_f1f2f3 |> write_genind_popsizes_to_table(COUNTS_DIR)

#---- END OF SCRIPT -----------------------------------------------------------#
