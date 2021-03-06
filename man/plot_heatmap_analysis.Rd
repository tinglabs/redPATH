% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biological_analysis_plots.R
\name{plot_heatmap_analysis}
\alias{plot_heatmap_analysis}
\title{Heatmap Analysis}
\usage{
plot_heatmap_analysis(
  sorted_matrix,
  pseudo_time,
  cell_labels = NULL,
  gene_cluster = NULL,
  heading = "HEATMAP",
  input_type = "hmm"
)
}
\arguments{
\item{sorted_matrix}{This can be either the HMM matrix or the sorted gene selection matrix obtained from earlier functions.}

\item{pseudo_time}{Pseudo development time}

\item{cell_labels}{Defaults to NULL. The cell type labels.}

\item{gene_cluster}{Defaults to NULL. It can be automatically caluclated with defualt gene cluster of 8. Or user can define the inferred gene cluster results using hierarchical clustering (get_gene_cluster) or by other means. Note that the input must match for the input matrix (ie gene / hmm)}

\item{heading}{Plot title}

\item{input_type}{"hmm" or "gene". Please specify type of sorted matrix input to the function.}
}
\description{
Comprehensive heatmap plot function with gene cluster and cell type labels.
}
\section{Biological Analysis}{

}

\examples{
selected_result <- gene_selection(full_exprs_matrix, redpath_pseudotime, dcor_result, mic_result, threshold = 0.5)
gene_state_result <- calc_all_hmm_diff_parallel(selected_result, cls_num = 2, threadnum = 4) 
inferred_gene_cluster <- get_gene_cluster(sorted_matrix = gene_state_result, cls_num = 8)
heatmap_plot <- plot_heatmap_analysis(sorted_matrix = inferred_gene_cluster, redpath_pseudotime, cell_labels = NULL, gene_cluster = inferred_gene_cluster, input_type = "hmm")
}
\keyword{plot_heatmap_analysis}
