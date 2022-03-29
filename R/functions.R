#' Retrieve gene orthologs/homologs
#'
#' Retrieve gene orthologs/homologs for a set of genes.
#' Converts between human and non-human analogs.
#'
#' @param genes A vector of gene symbols or Entrez/Ensembl IDs.
#' @param species Species name, such as \code{Mus musculus} or \code{mouse} (see \code{\link[=species]{species()}} for options).
#' @param human A logical scalar indicating if the input genes are human. If \code{TRUE}, the input genes are human. If \code{FALSE}, the input genes correspond to the non-human species and the output will be the human equivalents.
#' @param min_support Minimum number of supporting source databases.
#' @param top For each gene, output only the match with the highest support level if there are multiple hits.
#'
#' @return A data frame of gene pairs (human and given species).
#'
#' @references
#' Wright MW, Eyre TA, Lush MJ, Povey S, Bruford EA. HCOP: the HGNC comparison of orthology predictions search tool. \emph{Mamm Genome}. 2005 Nov;16(11):827-8. \doi{10.1007/s00335-005-0103-2}
#'
#' Eyre TA, Wright MW, Lush MJ, Bruford EA. HCOP: a searchable database of human orthology predictions. \emph{Brief Bioinform}. 2007 Jan;8(1):2-5. \doi{10.1093/bib/bbl030}
#'
#' Seal RL, Gordon SM, Lush MJ, Wright MW, Bruford EA. genenames.org: the HGNC resources in 2011. \emph{Nucleic Acids Res}. 2011 Jan;39:D514-9. \doi{10.1093/nar/gkq892}
#'
#' @importFrom methods is
#' @importFrom dplyr distinct group_by inner_join top_n ungroup
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' orthologs(genes = "TP53", species = "mouse", human = TRUE)
#' orthologs(genes = c("Ptprc", "Cd34"), species = "mouse", human = FALSE)
orthologs <- function(genes, species, human = TRUE, min_support = 3, top = TRUE) {
  # check if inputs are valid
  if (!is.vector(genes)) {
    stop("`genes` is not a character vector (can be a single gene)")
  }
  if (!is(species, "character")) {
    stop("`species` is not a character string")
  }
  if (!is.logical(human)) {
    stop("`human` is not TRUE/FALSE")
  }
  if (!is.numeric(min_support)) {
    stop("`min_support` is not a number")
  }
  if (!is.logical(top)) {
    stop("`top` is not TRUE/FALSE")
  }

  # subset the orthologs table to the relevant species
  species <- species(species)$taxon_id
  odf <- orthologs_df[orthologs_df$taxon_id == species, ]

  # determine the gene type used (for determining the top match)
  if (human) {
    if (any(odf$human_entrez %in% genes)) {
      gene_col <- "human_entrez"
      odf <- odf[odf$human_entrez %in% genes, ]
    } else if (any(odf$human_symbol %in% genes)) {
      if (any(odf$human_ensembl %in% genes)) {
        stop("both gene symbols and Ensembl IDs are used")
      }
      gene_col <- "human_symbol"
      odf <- odf[odf$human_symbol %in% genes, ]
    } else if (any(odf$human_ensembl %in% genes)) {
      gene_col <- "human_ensembl"
      odf <- odf[odf$human_ensembl %in% genes, ]
    } else {
      stop("no orthologs found or the genes are not valid human genes")
    }
  } else {
    if (any(odf$entrez %in% genes)) {
      gene_col <- "entrez"
      odf <- odf[odf$entrez %in% genes, ]
    } else if (any(odf$symbol %in% genes)) {
      if (any(odf$ensembl %in% genes)) {
        stop("both gene symbols and Ensembl IDs are used")
      }
      gene_col <- "symbol"
      odf <- odf[odf$symbol %in% genes, ]
    } else if (any(odf$ensembl %in% genes)) {
      gene_col <- "ensembl"
      odf <- odf[odf$ensembl %in% genes, ]
    } else {
      stop("no orthologs found or the genes are not valid")
    }
  }

  # filter by the number of supporting databases
  odf <- odf[odf$support_n >= min_support, ]

  # select only the best match for each gene
  if (top) {
    odf <- dplyr::group_by(odf, .data[[gene_col]], .data[["taxon_id"]])
    odf <- dplyr::top_n(odf, n = 1, wt = .data[["support_n"]])
    odf <- dplyr::ungroup(odf)
  }

  # clean up for output
  odf <- dplyr::distinct(odf)
  odf <- as.data.frame(odf, stringsAsFactors = FALSE)

  return(odf)
}

#' Retrieve the available species
#'
#' List the species with available human orthologs.
#'
#' @param species Species name, such as \code{Mus musculus} or \code{mouse}. If specified, will return results for the given species only.
#'
#' @return A data frame of the available species.
#'
#' @importFrom methods is
#' @importFrom dplyr distinct full_join group_by mutate
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' species()
#' species("Mus musculus")
#' species("mouse")
#' species("rat")
species <- function(species = NULL) {
  # generate a more readable version of the species table
  sci_names <- species_df[species_df$name_class == "scientific name", 1:2]
  colnames(sci_names) <- c("taxon_id", "scientific_name")
  common_names <- species_df[species_df$name_class != "scientific name", 1:2]
  common_names <- dplyr::group_by(common_names, .data[["taxon_id"]])
  common_names <- dplyr::mutate(common_names, common_name = toString(.data[["name"]]))
  common_names <- dplyr::distinct(common_names[, c("taxon_id", "common_name")])
  sdf <- dplyr::full_join(sci_names, common_names, by = "taxon_id")
  sdf <- as.data.frame(sdf)

  # return the full table if the species name is not specified
  if (is.null(species)) {
    return(sdf)
  }

  # check if the species name is valid if it is specified
  if (!is(species, "character")) {
    stop("`species` is not a character string")
  }
  if (!species %in% species_df$name) {
    stop("unknown species, must be one of: ", toString(species_df$name))
  }

  # subset to the given species
  species_id <- species_df[species_df$name == species, "taxon_id"]
  sdf <- sdf[sdf$taxon_id == species_id, ]

  return(sdf)
}
