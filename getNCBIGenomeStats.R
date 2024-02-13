library(here)
library(taxizedb)
library(biomartr)

readRenviron("~/.Renviron")

# Get a list of all accession numbers for species-level Molluscs
id_txdb <- taxizedb::name2taxid('Phoronida', db = "ncbi")

x <- taxizedb::downstream(id_txdb, db = "ncbi", downto = "species")

ids <- x[[id_txdb]][["childtaxa_id"]]

# download genome assembly stats file for Homo sapiens
for (s in ids) {
  getAssemblyStats(db  = "genbank", 
                   organism = s, 
                   reference = FALSE,
                   path = file.path("GenBank","Phoronida_AssemblyStats"))
}

 # to check position of id last downloaded, e.g.
which(ids=='1068534')
[1] 91171

# use output to re-start download list from the middle, e.g.
ids_v2 <- ids[91172:182437]

# Scatter plots of genome vs gene nucleotide content in NCBI 
library(ggplot2)
library(ggrepel)
library(wesanderson)

nts <- read.csv(file=here("genomestats-v2.csv"))

# Plot all phyla
p <- ggplot(nts, aes(x=Total.NCBI.genome.nts.average.metagenome.nts, 
                     y=Total.NCBI.18S.nts.Average.sample.18S.nts, 
                     label=Phylum)) + 
  geom_point(color="blue") +
  geom_text_repel(data=subset(nts, 
                              Total.NCBI.18S.nts.Average.sample.18S.nts > 0.1 |
                                Total.NCBI.genome.nts.average.metagenome.nts > 160)) +
  coord_cartesian(expand=FALSE, ylim=c(0, 0.4), xlim=c(0, 170)) + #,  +
  theme_classic()
 
ggsave(here("genome_vs_18Sgene_nt_comparison-norm.png"), p,
        width=7, height=4, dpi=300, units="in")

# Plot phyla with lower nt content
p <- ggplot(nts, aes(x=Total.NCBI.genome.nts.average.metagenome.nts,
                     y=Total.NCBI.COI.nts.Average.sample.COI.nts, 
                     label=Phylum)) + 
  geom_point(color="blue") +
  geom_text_repel() +
  scale_x_continuous(limits=c(0, 3.5)) + 
  scale_y_continuous(limits=c(0, 0.05)) +
  coord_cartesian(expand=FALSE) +
  theme_classic()

ggsave(here("genome_vs_COIgene_nt_comparison-norm-smallsubset-1.png"), p,
       width=7, height=4, dpi=300, units="in")

# Plot phyla with even lower nt content
p <- ggplot(nts, aes(x=Total.NCBI.genome.bps, y=Total.NCBI.18S.bps, label=Phylum)) + 
  geom_point(color="blue") +
  geom_text_repel() + 
  scale_x_continuous(limits=c(0, 1E10)) + 
  scale_y_continuous(limits=c(0,50000)) +
  theme_classic()

ggsave(here("genome_vs_18Sgene_nt_comparison_smallsubset-2.png"), p,
       width=7, height=4, dpi=300, units="in")

# Plot ratio of average sample bps:available NCBI bps
p <- ggplot(nts, aes(x=Average...sample.COI.bps.NCBI.COI.bps, 
                     y=Average...sample.18S.bps.NCBI.18S.bps, 
                     label=Phylum)) + 
  geom_point(color="blue") +
  geom_text_repel() + 
  scale_x_continuous(limits=c(0, 20000)) + 
  scale_y_continuous(limits=c(0, 20000)) +
  theme_classic()

ggsave(here("COI_vs_18S_samplebps_over_NCBIbps-smallsubset.png"), p,
       width=7, height=4, dpi=300, units="in")

# Violin plots of NCBI genome bps vs NCBI gene bps 

# Need to convert to long format
nts_long <- nts %>% 
  gather(key="Ratio", 
         value="Val", 
         Total.NCBI.genome.COI.bps, 
         Total.NCBI.genome.18S.bps)

p1 <- ggplot(nts_long, aes(x=Ratio, 
                     y=Val, 
                     fill=Ratio)) + # fill=name allow to automatically dedicate a color for each group
  geom_violin() +
  scale_y_continuous(name="Ratio", 
                     trans=scales::pseudo_log_trans(base = 10),
                     breaks=c(100000, 2000000, 10000000),
                     labels = scales::comma ) +
  geom_point(position = position_jitter(seed = 1, width = 0.1)) +
  coord_cartesian(expand=FALSE) +
  theme_classic()

p2 <- ggplot(nts_long, aes(x=Ratio, 
                     y=Val, 
                     fill=Ratio)) + # fill=name allow to automatically dedicate a color for each group
  geom_violin() +
  geom_point(position = position_jitter(seed = 1, width = 0.1)) +
  scale_y_continuous(limits=c(0,2000000)) +
  coord_cartesian(expand=FALSE) +
  theme_classic()

p3 <- ggplot(nts_long, aes(x=Ratio, 
                     y=Val, 
                     fill=Ratio,
                     label=Phylum)) + # fill=name allow to automatically dedicate a color for each group
  geom_violin() +
  geom_point(position = position_jitter(seed = 1, width = 0.1)) +
  geom_text_repel() +
  scale_y_continuous(limits=c(750000,20000000)) +
  coord_cartesian(expand=FALSE) +
  theme_classic()

ggsave(here("NCBI_genome-to-gene-ratio_top.png"), p3,
       width=7, height=4, dpi=300, units="in")


