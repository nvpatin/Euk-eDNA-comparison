library(here)
library(taxizedb)
library(biomartr)

readRenviron("~/.Renviron")

# Get a list of all accession numbers for species-level Molluscs
id_txdb <- taxizedb::name2taxid('Streptophyta', db = "ncbi")

x <- taxizedb::downstream(id_txdb, db = "ncbi", downto = "species")

ids <- x[[id_txdb]][["childtaxa_id"]]

# download genome assembly stats file for Homo sapiens
for (s in ids_v2) {
  getAssemblyStats(db  = "genbank", 
                   organism = s, 
                   reference = FALSE,
                   path = file.path("GenBank","Streptophyta_AssemblyStats"))
}

 # to check position of id last downloaded, e.g.
which(ids=='1068534')
[1] 91171

# use output to re-start download list from the middle, e.g.
ids_v2 <- ids[91172:182437]

