
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(stringr)
library(glue)
library(usethis)

# Import HCOP orthologs ---------------------------------------------------

# Import the human and ortholog data from all HCOP species
hcop_txt_url <- "ftp://ftp.ebi.ac.uk/pub/databases/genenames/hcop/human_all_hcop_sixteen_column.txt.gz"
hcop <- read_tsv(hcop_txt_url, col_types = cols())
format(object.size(hcop), units = "Mb")
nrow(hcop)

# Check if the HCOP table is too small
if (nrow(hcop) < 865000) stop("HCOP file may be truncated")

# Clean up the table of orthologs
orthologs <-
  hcop %>%
  select(
    human_symbol,
    human_entrez = human_entrez_gene,
    human_ensembl = human_ensembl_gene,
    taxon_id = ortholog_species,
    symbol = ortholog_species_symbol,
    entrez = ortholog_species_entrez_gene,
    ensembl = ortholog_species_ensembl_gene,
    support
  ) %>%
  mutate(
    human_symbol = na_if(human_symbol, "-"),
    human_entrez = na_if(human_entrez, "-"),
    human_ensembl = na_if(human_ensembl, "-"),
    symbol = na_if(symbol, "-"),
    entrez = na_if(entrez, "-"),
    ensembl = na_if(ensembl, "-")
  ) %>%
  mutate(
    human_entrez = as.integer(human_entrez),
    entrez = as.integer(entrez),
    taxon_id = as.integer(taxon_id)
  ) %>%
  arrange(human_symbol, human_ensembl, taxon_id, symbol, ensembl) %>%
  distinct()
nrow(orthologs)

# Remove repeating databases (some may be listed multiple times)
orthologs <- mutate(orthologs, support = strsplit(support, ",", fixed = TRUE))
orthologs <- mutate(orthologs, support = map(support, sort))
orthologs <- mutate(orthologs, support = map(support, unique))
orthologs <- mutate(orthologs, support = map_chr(support, paste, collapse = "|"))

# Check the number of supporting databases per ortholog pair
orthologs <- mutate(orthologs, support_n = str_count(support, "\\|") + 1)
orthologs %>% count(support_n)

# Filter ortholog pair with only one supporting database
orthologs <- filter(orthologs, support_n > 1)

# Check orthologs stats
n_distinct(orthologs$human_entrez)
n_distinct(orthologs$human_ensembl)
orthologs %>%
  group_by(taxon_id) %>%
  summarize(n_distinct(human_symbol), n_distinct(symbol), max(support_n))

# Import MGI orthologs ----------------------------------------------------

# Import the human and mouse data from MGI
mgi_txt_url <- "http://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt"
mgi <- read_tsv(mgi_txt_url, col_names = FALSE, col_types = cols())
format(object.size(mgi), units = "Mb")
nrow(mgi)

# Clean up the table of orthologs
mgi_orthologs <-
  mgi %>%
  select(human_symbol = X1, mouse_symbol = X5) %>%
  distinct()

# Check the number of genes
n_distinct(mgi_orthologs$human_symbol)
n_distinct(mgi_orthologs$mouse_symbol)

# Import taxonomy data ----------------------------------------------------

# Download the taxonomy zip file
taxdmp_url <- "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip"
download.file(url = taxdmp_url, destfile = "taxdmp.zip")

# Extract the taxonomy names file
unzip("taxdmp.zip", files = "names.dmp")

taxdmp <- read_delim("names.dmp", delim = "|", trim_ws = TRUE, col_names = FALSE, col_types = cols())

# Delete the taxonomy zip file and its contents since they are no longer needed
file.remove("taxdmp.zip")
file.remove("names.dmp")

# Extract the names for the relevant species
species <-
  taxdmp %>%
  select(taxon_id = X1, name = X2, name_class = X4) %>%
  mutate(taxon_id = as.integer(taxon_id)) %>%
  filter(taxon_id %in% orthologs$taxon_id) %>%
  filter(str_detect(name_class, "name")) %>%
  arrange(name)

# Check that all HCOP species are present
n_distinct(species$taxon_id)
if (length(setdiff(orthologs$taxon_id, species$taxon_id))) stop("species mismatch")

# Check the number of genes per species
species %>%
  filter(name_class == "scientific name") %>%
  full_join(orthologs, by = "taxon_id") %>%
  group_by(name, taxon_id) %>%
  summarize(n_distinct(ensembl))

# Prepare package ---------------------------------------------------------

# Convert to standard data frames (dplyr pipeline sets them as tibbles)
orthologs_df <- as.data.frame(orthologs)
mgi_orthologs_df <- as.data.frame(mgi_orthologs)
species_df <- as.data.frame(species)

# Check the size of final tables
format(object.size(orthologs_df), units = "Mb")
format(object.size(mgi_orthologs_df), units = "Mb")
format(object.size(species_df), units = "Kb")

# Create package data
use_data(
  orthologs_df,
  mgi_orthologs_df,
  species_df,
  internal = TRUE,
  overwrite = TRUE,
  compress = "xz"
)
