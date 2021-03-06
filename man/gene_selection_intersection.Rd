% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biological_analysis_utils.R
\name{gene_selection_intersection}
\alias{gene_selection_intersection}
\title{Selects significant genes}
\usage{
gene_selection_intersection(
  full_gene_exprs,
  pseudo_time,
  correlation_vector,
  correlation_vector2 = NULL,
  threshold = 0.5
)
}
\arguments{
\item{full_gene_exprs}{A gene (row) by cell (col) gene expression matrix with all genes.}

\item{pseudo_time}{The pseudo development time}

\item{correlation_vector}{The calculated correlation vector of each gene (from calc_correlation)}

\item{correlation_vector2}{Genes will be selected if the threshold satisfies BOTH correlation vectors.}

\item{threshold}{The threshold cut-off to retain significant genes. Defaults to 0.5}
}
\description{
This function selects the significant genes from 2 correlation vectors (intersecting genes)
}
\section{Biological Analysis}{

}

\examples{
dcor_result <- calc_correlation(full_exprs_matrix, redpath_pseudotime, type = "dcor") 
mic_result <- calc_correlation(full_exprs_matrix, redpath_pseudotime, type = "mic")
selected_result <- gene_selection(full_exprs_matrix, redpath_pseudotime, dcor_result, mic_result, threshold = 0.5)
}
\keyword{gene_selection_intersection}
