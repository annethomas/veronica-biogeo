# Veronica uplift BGB analyses protocols

Done externally: phylogeny (StarBEAST2 on CIPRES), band assignment (done manually based on floras and expert opinion), and TTR runs (Matt Larcombe)

Bold scripts (.R) are sequential steps in a pipeline (or sourced by these scripts); some were run on departmental cluster, others on local machine, so filepaths may vary

## BioGeoBEARS + binning (dir: biogeobears)
- **run_bgb_analysis_cluster.R** (cluster only)
  - run BioGeoBEARS and BSMs and/or bin and plot cladogenetic/anagenetic events and rates for tree set
  - input: posterior tree set file, habitat file, a few other BGB input files
  - options: 
    - designate input files for two-band (merged) or three-band dataset
    - PREP_TREES=T/F
      - selects 100 random trees and formats using **prep_trees.R** (only do once)
    - BGB=T/F; BSM=T/F 
      - run Biogeobears and BSMs for tree set with run_BioGeoBEARS_plusBSM_multitrees() from biogeobears_functions.R 
        - runs run_biogeobears() and run_bgb_BSM() from **biogeobears_functions.R** (see below)
      - will check if BGB (and BSM) output files already exist, so if you want to run BSMs on existing files, both must be true (but it won’t rerun the models). BSM=F will only run the BGB models. 
      - comment toggle: parallel run function (cluster; haven’t tested lately) or linear run function
    - BIN=T/F; PLOTONLY=T/F
      - will generate binned anangentic/cladogenetic rate/event plots from the BSMs pooled across all the trees and BSMs (if BSMs have already been run, choose PLOTONLY)
      - runs run_bin_plot_BSM_multi() from binning/**bin_plot_BSMs.R** 
      - includes sortdate50_ASL and sortadate50_merged; also includes a couple of other scenarios that may or may not be tested
      - plots are basic; more refined plots on laptop binning_figures.R and final versions in uplift_biogeo_ms_plots.R
- run_bgb_single.R 
  - optional/util; run the 6 BGB models for a single (usually the MCC) tree 
  - runs run_biogeobears() from biogeobears_functions.R
- **biogeobears_functions.R** (sourced)
  - run_BioGeoBEARS_plusBSM_multitrees()
    - run 6 BioGeoBEARS models (Nick Matzke’s scriptrun_biogeobears()) and generate BSMs for posterior distribution (or other set) of trees in a directory (also plots output)
    - also checks for existing model files and saves AIC scores in table
  - input: tree directory; run_BioGeoBEARS_plusBSM_multitrees_listed() takes vector of tree filenames
- binning/**bin_plot_BSMs.R** (sourced)
  - binning functions, binning modification functions (e.g. padNAs), plotting functions
- binning/run_bin_plot_BSMs.R (laptop)
  - optional/util; commands for run_bin_plot_BSM_multi() from bin_plot_BSMs.R for several run directories
  - option for 100 trees (if downloaded from cluster) but file paths need to be updated; prefer cluster
- binning/binning_figures.R (laptop)
  - manuscript plots for binning figures, two-band and three-band (final version in uplift_biogeo_ms_plots.R) 
  - individual plots and paneled
  - plot_dir=file.path(bgb_dir,”biogeobears_output_100trees/sortadate50_merged/binned_plots/bluegreen/”); sortadate50_ASL/binned_plots/plots_2.0 for three-band
  - deals with median line before habitat is available by changing 0s to NAs (instead of doing this in the binning calculations to allow the median to reflect sample size in outlier tree ages, see next point)
  - deals with outlier tree ages by truncating the plots where the median is 0
  - where it works, uses plot_binned_multi_var() from bin_plot_BSMs.R to combine clado and dispersal plots (only two-band events currently); otherwise plot_binned_single_var()
  - includes some exploratory plotting and color manipulation from past (labelled ‘obsolete’) 
- my_plot_BioGeoBEARS_results.R (laptop)
  - modified version of BioGeoBEARS function for tree plotting with customized arguments/visuals
  - sourced in TTR/tip_trait_plotting.R, biogeobears/ms_tree_fig.R, and biogeobears/ms_ASL_tree_fig.R

## TTR traits (dir: TTR)
- **PCA_select_TTRtraits.R**: phylogenetic PCA to reduce dimensionality of TTR traits prior to mvMORPH
  - automatic eigenvector threshold followed by manual filtering based on PC weights and biplot correlation between variables
- prep_trait_tables.R
  - generate several modified/expanded trait/species tables from existing ones and save/resave csv (only run once)
  - alt_lat_bands.csv combines altitude/latitude info with assigned habitat bands
  - Veronica.sdm.stats.corrected.csv has unconstrained trait values (outside the transformed window and thus realistic biological values) replaced with the max value for both transformed and back-transformed columns
  - traits_bands_back.csv has the selected back-transformed traits joined to assigned bands
- compare_TTR_traits.R
  - look at the distribution of traits for different habitat groups in boxplots, phylo anova, phylo regression
  - draft manuscript boxplot of significant traits
- false_neg_pos_rates.R
  - calculate overall TTR model fit based on false positive and negative rates (from existing stats table) 
- tip_trait_plotting.R
  - plot BGB tree with selected traits colored on tips
  - sources my_plot_BioGeoBEARS_results.R, nodiv_plotting_functions_edited.R, nodiv_color_functions.R
  - includes some exploratory attempts for reference

## mvMORPH (dir: niche_evolution)
- **prep_mvMORPH.R** (laptop): Prepare tree + trait input to run mvMORPH on the cluster, including simmaps for OUM
  - Format TTR traits from manually saved (?) nexus file alongside tree and subset based on traits selected by eigenvector loadings in PCA_select_TTRtraits.R
    - (using Claddis functions in time_slice/**ancestral_slice_functions.R** to prepare for downstream analyses)
    - Sources make.simmap.BGB.original.R and rpanda_stratified_BGB_to_tables.R 
  - choose single or multi-regime model type
  - single: make simple input list with tree and trait matrix
  - multi: make simmaps from BGB BSMs and create input list for each BSM
  - upload input to cluster directories described below
- run_mvMORPH_single.R (cluster)
  - args: 
    - traitset, eg “TTR9”
    - model, ie “BM”, “OU”, or “EB”
  - reads input file from prep_mvMORPH.R at file.path(TTR_dir,paste0("mvmorph_input_",traitset,".Rdata"))
  - output: list with fit and ancestral estimations in file.path(TTR_dir,"mvmorph_output",paste0("mvmorph_output_",traitset,model,".Rdata"))
- **run_mvMORPH_multiBSM.R** (cluster)
  - args
    - traitset=commandArgs(trailingOnly=TRUE)[1] ##TTR4, TTR9
    - model=commandArgs(trailingOnly=TRUE)[2] ## BM, OU
    - dirname=commandArgs(trailingOnly=TRUE)[3] ## e.g. sortadate50_merged_TTRtips
    - bsm_dirname=commandArgs(trailingOnly=TRUE)[4] ## e.g. BSMs_run1
    - BSM_range=commandArgs(trailingOnly=TRUE)[5] ## BSM number(s), comma separated in quotes
  - BSM_range allows for subsets and for jobs to be run in batches, as it runs in a single loop through named BSMs
  - cluster currently doesn’t support ape on all nodes so best to run via tmux on node13; run_mvMORPH_multiBSM.sh has a toggle for node7 or node13 and has two arguments for start/end of BSMrange 
- compare_mvMORPH_models.R
  - look at AICc and plot values on trees
  - check BSMs from OUM for reliability

## Time slice disparity (dir: time_slice)
- **time_slice_BSMs.R**
  - choose a dataset (merged or ASL) and mvMORPH model fit (OU or OUM, hypothetically BM or EB as well) and run the time slice pipeline for each BSM of MCC tree
    - note: the median distance calculations and summary stats/models in run_slice() etc are currently calculated from the whole set of nodes in the designated focal areas (that is for the whole dist_df). If only interested in node-by-node distances/disparity in dist_df, set distdf_only to TRUE 
    - If interested in plotting both mountain and lowland disparity but not combined stats, set ML_distdf to TRUE
    - For a full run focused on mountains, set ML_distdf to FALSE.
  - if OU, loads a single model and morphospace for all BSMs; if OUM, will load individual mvMORPH files and generate morphospace for each BSM
  - run_slice() for each BSM from ancestral_slice_functions.R:
    - slice_distance() does the main time slice analysis on the input tree
    - shuffle_slice() generates null distrubtion (100 tip-shuffles and slice_distance() for each)
      - currently has the same seed set for each function call
    - plot_dist_obs_null() generates p-values from the null distribution for various stats—lots of toggle options in the ancestral_slice_functions.R arguments, as well as options whether to plot and/or save plots to pdf; may need to go directly to change defaults in function or add arguments to this script
    - includes output examining relationships of node distance/disparity with time for cladogenesis vs dispersal with scatter/boxplots, linear model, and correlation coefficients
    - manuscript only uses median niche distance and the difference of this between clado and dispersal, and overall disparity
- plot_niche_distance_disparity.R: generates the two main analyses/figures from the time slice output for the paper, median niche distance between clado and dispersal, and disparity through time (final versions in uplift_biogeo_ms_plots.R)
  - currently designed for OUM BSMs
  - niche distance and pvalues are extracted from the BSMs (for mountain-only run, currently a separate output Rdata file from the ML one) into a data frame, saved in slice_results<>.csv (long and wide), median BSM values plotted for clado and dispersal 
  - time slice disparity currently plotted for mountains and lowlands (the ML output file) across BSMs
  - also has exploratory options for looking at the null distributions for the BSMs (including M alone)

## Range and niche overlap (dir: niche_overlap)
- generate TTR projections from Matt Larcombe’s TTR fits: 01-subset_extent_latlon.R, 02-prepare_env_grid.R, 08-predict_species_dist.R on cluster
- below scripts source overlap_util_functions.R
- **rangeOverlap_jaccard_phyloclim.R**: calculate buffers around observed occurrences and overlap between buffer polygons
  - phyloclim calculates Schoener’s D from rasterized buffers (necessary for age-range correlation)
  - jaccard calculates spatial overlap ratio from original polygons (just a possible alternative for non phyloclim calculations)
- **nicheOverlap_spaa_phyloclim.R**: calculate niche overlap from TTR projections (rasters)
  - phyloclim calculates Schoener’s D
  - spaa does exactly the same but doesn’t produce the output format needed for ARC (niolap)
  - **process_ARC_range_niche.R**
    - first create input list of overlap matrix and tree pruned as desired for each of niche and range
    - **cluster_run_age_range_correlation.R** on the cluster with input list; takes hours to do the 100 null simulations
      - uses age_range_correlation_warren.R with modifications to the info saved
    - various plots (traditional scatterplots, histograms compared to null medians, range vs niche correlation)
    - sister pairs: another option for scatter plot comparison (no null simulation)

## dtt 
see uplift_biogeo_ms_plots.R

## rpanda
- rpanda_models.R: run various rpanda models with elevation and time

## util
tree plotting functions

