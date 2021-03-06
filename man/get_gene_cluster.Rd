% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biological_analysis_plots.R
\name{get_gene_cluster}
\alias{get_gene_cluster}
\title{Gene Clustering}
\usage{
get_gene_cluster(sorted_matrix, cls_num = 8)
}
\arguments{
\item{sorted_matrix}{This can be either the HMM matrix or the sorted gene selection matrix obtained from earlier functions.}

\item{cls_num}{The number of gene cluster desired}
}
\description{
Simple gene cluster on the significant genes sorted matrix (based on dCor and MIC). Can be either HMM matrix or the gene expression matrix
}
\section{Biological Analysis}{

}

\examples{
selected_result <- gene_selection(full_exprs_matrix, redpath_pseudotime, dcor_result, mic_result, threshold = 0.5)
gene_state_result <- calc_all_hmm_diff_parallel(selected_result, cls_num = 2, threadnum = 4) 
inferred_gene_cluster <- get_gene_cluster(sorted_matrix = gene_state_result, cls_num = 2)
heatmap_plot <- plot_heatmap_analysis(sorted_matrix = inferred_gene_cluster, redpath_pseudotime, cell_labels = NULL, gene_cluster = inferred_gene_cluster, input_type = "hmm")
}
\keyword{get_gene_cluster}
