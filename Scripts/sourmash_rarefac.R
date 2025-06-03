# Load libraries
library(dplyr)
library(sourmashconsumr)
library(ggplot2)
library(RColorBrewer)

setwd("/Users/nastassiapatin/OneDrive - SCCWRP/RREAS Time Series/Metagenomes/sourmash/v3_nocontrols_outliers/sigs/RL2019")

## One file
# specify path 
signature_df <- read_signature("1903_buk5-1_interleaved.sig_sub.sig")

# change 'filename' column to 'name'
sig_df <- rename(signature_df, name = filename)

# rarefaction of one file
rarefaction_df <- from_signatures_to_rarefaction_df(signatures_df = sig_df)

### Multiple files
# read in paths
sigs <- list.files("RL2019", pattern="sub")

# rarefaction of multiple files
setwd("/Users/nastassiapatin/OneDrive - SCCWRP/RREAS Time Series/Metagenomes/sourmash/v3_nocontrols_outliers/sigs/RL2019")
signature_df <- read_signature(sigs)

sig_df <- rename(signature_df, name = filename)

rarefaction_df <- from_signatures_to_rarefaction_df(signatures_df = sig_df)

# plot
rarefaction_plt <- plot_signatures_rarefaction(rarefaction_df = rarefaction_df,
                                               fraction_of_points_to_plot = 1)

p <- rarefaction_plt +
  ggplot2::theme_classic() +
  ggplot2::geom_point(size=0.5)

ggsave("LASKER2019_sourmash_rarefac.png",
       p, width=9, height=5, dpi=300, units="in")

# multiple rarefaction dfs
rarefaction_df_CANON2018 <- rarefaction_df
rarefaction_df_Lasker2018 <- rarefaction_df
rarefaction_df_Lasker2019 <- rarefaction_df



