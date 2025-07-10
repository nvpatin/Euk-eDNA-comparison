## Katie Pitz
## 01/03/25
## Rarefaction curves

# Set directory to save plots
directory <- '/Users/kpitz/github/nvpatin/Euk-eDNA-comparison/katie_sandbox/'


# Load Libraries -----------------------------------------------------------------

# install ranacapa
# library(devtools)
# devtools::install_github("gauravsk/ranacapa")

library(phyloseq)
library(vegan)
library(ranacapa)  #for rarefaction curve

# Create Phyloseq Object and rarefaction curves ----------------------------

## 18S -------------------------------------
asv_file <- '/Users/kpitz/github/nvpatin/Euk-eDNA-comparison/Dada2_asv_data/Euk-eDNA_18S_asv_match.csv'
taxa_file <- '/Users/kpitz/github/nvpatin/Euk-eDNA-comparison/Dada2_asv_data/Euk-eDNA_18S_taxa_match.csv'
meta_file <- '/Users/kpitz/github/nvpatin/Euk-eDNA-comparison/Dada2_asv_data/Euk-eDNA_18S_meta_match.csv'

phylo_object <- merge_phyloseq(otu_table(as.matrix(read.csv(file = asv_file, row.names = 1,check.names = FALSE)),taxa_are_rows = TRUE),
                           tax_table(as.matrix(read.csv(file = taxa_file, row.names = 1))),
                           sample_data(read.csv(file = meta_file, row.names = 1)))

#ranacapa rarefaction curves
p <- ggrare(phylo_object, step = 1000, label = "SAMPLING_cruise", color = "PlateID",
            plot = TRUE, parallel = FALSE, se = FALSE)
p

p <- ggrare(phylo_object, step = 1000, color = "PlateID",
            plot = TRUE, parallel = FALSE, se = FALSE)
p
filename = paste(directory, '18S_Ranacapa_rarefaction_curve_PlateID.png', sep='')
print(filename)
ggsave(filename,height = 6, width =10, units = 'in')

p <- ggrare(phylo_object, step = 1000, color = "SAMPLING_cruise",
            plot = TRUE, parallel = FALSE, se = FALSE)
p
filename = paste(directory, '18S_Ranacapa_rarefaction_curve_cruise.png', sep='')
print(filename)
ggsave(filename,height = 6, width =10, units = 'in')

p <- ggrare(phylo_object, step = 1000,
            plot = TRUE, parallel = FALSE, se = FALSE)
p
filename = paste(directory, '18S_Ranacapa_rarefaction_curve.png', sep='')
print(filename)
ggsave(filename,height = 6, width =10, units = 'in')

### Plots by Cruise ----------------------------
lasker_2019 = subset_samples(phylo_object, SAMPLING_cruise %in% c('1903', '1903.0'))
CN18_phylo = subset_samples(phylo_object, SAMPLING_cruise %in% c('CN18S','CN18F'))
lasker_2018 = subset_samples(phylo_object, SAMPLING_cruise %in% c('LASKER2018'))

p <- ggrare(lasker_2019, step = 1000,
            plot = TRUE, parallel = FALSE, se = FALSE)
p
filename = paste(directory, '18S_Ranacapa_rarefaction_curve_Lasker2019.png', sep='')
print(filename)
ggsave(filename,height = 6, width =10, units = 'in')

p <- ggrare(lasker_2018, step = 1000,
            plot = TRUE, parallel = FALSE, se = FALSE)
p
filename = paste(directory, '18S_Ranacapa_rarefaction_curve_Lasker2018.png', sep='')
print(filename)
ggsave(filename,height = 6, width =10, units = 'in')

p <- ggrare(CN18_phylo, step = 1000,
            plot = TRUE, parallel = FALSE, se = FALSE)
p
filename = paste(directory, '18S_Ranacapa_rarefaction_curve_CANON2018.png', sep='')
print(filename)
ggsave(filename,height = 6, width =10, units = 'in')

## COI -------------------------------------
asv_file <- '/Users/kpitz/github/nvpatin/Euk-eDNA-comparison/Dada2_asv_data/Euk-eDNA_COI_asv_match.csv'
taxa_file <- '/Users/kpitz/github/nvpatin/Euk-eDNA-comparison/Dada2_asv_data/Euk-eDNA_COI_taxa_match.csv'
meta_file <- '/Users/kpitz/github/nvpatin/Euk-eDNA-comparison/Dada2_asv_data/Euk-eDNA_COI_meta_match.csv'

phylo_object <- merge_phyloseq(otu_table(as.matrix(read.csv(file = asv_file, row.names = 1,check.names = FALSE)),taxa_are_rows = TRUE),
                               tax_table(as.matrix(read.csv(file = taxa_file, row.names = 1))),
                               sample_data(read.csv(file = meta_file, row.names = 1)))

#ranacapa rarefaction curves
p <- ggrare(phylo_object, step = 1000, label = "SAMPLING_cruise", color = "PlateID",
            plot = TRUE, parallel = FALSE, se = FALSE)
p

p <- ggrare(phylo_object, step = 1000, color = "PlateID",
            plot = TRUE, parallel = FALSE, se = FALSE)
p
filename = paste(directory, 'COI_Ranacapa_rarefaction_curve_PlateID.png', sep='')
print(filename)
ggsave(filename,height = 6, width =10, units = 'in')

p <- ggrare(phylo_object, step = 1000, color = "SAMPLING_cruise",
            plot = TRUE, parallel = FALSE, se = FALSE)
p
filename = paste(directory, 'COI_Ranacapa_rarefaction_curve_cruise.png', sep='')
print(filename)
ggsave(filename,height = 6, width =10, units = 'in')

p <- ggrare(phylo_object, step = 1000,
            plot = TRUE, parallel = FALSE, se = FALSE)
p
filename = paste(directory, 'COI_Ranacapa_rarefaction_curve.png', sep='')
print(filename)
ggsave(filename,height = 6, width =10, units = 'in')

p <- ggrare(phylo_object, step = 1000, label = 'ASV.sample', 
            plot = TRUE, parallel = FALSE, se = FALSE)
p
filename = paste(directory, 'COI_Ranacapa_rarefaction_curve_sampleID.png', sep='')
print(filename)
ggsave(filename,height = 6, width =10, units = 'in')

# remove sample with large number of reads, CN19_Buck6_Dav_Night_0m_1

newPhylo = subset_samples(phylo_object, ASV.sample != 'CN19_Buck6_Dav_Night_0m_1')
#newPhylo = subset_samples(phylo_object, ASV.sample %in% c('CN19_Buck6_Dav_Night_0m_1','CN19_Buck6_Dav_Night_0m_2'))
p <- ggrare(newPhylo, step = 1000, label = 'ASV.sample', 
            plot = TRUE, parallel = FALSE, se = FALSE)
p

# plots without high read # sample:
p <- ggrare(newPhylo, step = 1000, color = "PlateID",
            plot = TRUE, parallel = FALSE, se = FALSE)
p
filename = paste(directory, 'COI_Ranacapa_rarefaction_curve_PlateID_samplim.png', sep='')
print(filename)
ggsave(filename,height = 6, width =10, units = 'in')

p <- ggrare(newPhylo, step = 1000, color = "SAMPLING_cruise",
            plot = TRUE, parallel = FALSE, se = FALSE)
p
filename = paste(directory, 'COI_Ranacapa_rarefaction_curve_cruise_samplim.png', sep='')
print(filename)
ggsave(filename,height = 6, width =10, units = 'in')

p <- ggrare(newPhylo, step = 1000,
            plot = TRUE, parallel = FALSE, se = FALSE)
p
filename = paste(directory, 'COI_Ranacapa_rarefaction_curve_samplim.png', sep='')
print(filename)
ggsave(filename,height = 6, width =10, units = 'in')

### Plots by Cruise ----------------------------
lasker_2019 = subset_samples(phylo_object, SAMPLING_cruise %in% c('1903', '1903.0'))
CN18_phylo = subset_samples(phylo_object, SAMPLING_cruise %in% c('CN18S','CN18F'))
lasker_2018 = subset_samples(phylo_object, SAMPLING_cruise %in% c('LASKER2018'))

p <- ggrare(lasker_2019, step = 1000,
            plot = TRUE, parallel = FALSE, se = FALSE)
p
filename = paste(directory, 'COI_Ranacapa_rarefaction_curve_Lasker2019.png', sep='')
print(filename)
ggsave(filename,height = 6, width =10, units = 'in')

p <- ggrare(lasker_2018, step = 1000,
            plot = TRUE, parallel = FALSE, se = FALSE)
p
filename = paste(directory, 'COI_Ranacapa_rarefaction_curve_Lasker2018.png', sep='')
print(filename)
ggsave(filename,height = 6, width =10, units = 'in')

p <- ggrare(CN18_phylo, step = 1000,
            plot = TRUE, parallel = FALSE, se = FALSE)
p
filename = paste(directory, 'COI_Ranacapa_rarefaction_curve_CANON2018.png', sep='')
print(filename)
ggsave(filename,height = 6, width =10, units = 'in')
