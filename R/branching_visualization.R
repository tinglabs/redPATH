# redPATH for Branching Processes:

# Idea:
# Record all the possible transitions from cell to cell / group of cells to group of cells (from each HMT solution)
# Normalize and assign weight
# Define an N by N transition matrix for visualization
theme_redpath <- function(base_size=16, base_family="sans serif") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(size = rel(1.2)),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text.y = element_text(),#element_text(angle=75, vjust = 0.5), 
            axis.text.x = element_text(),
            axis.line = element_line(colour="black", size = 1),
            axis.ticks.x = element_line(),
            axis.ticks.y = element_line(),
            panel.grid.major = element_blank(),#element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            #legend.position = c(.9,.15), #"bottom", #(.9, .25)
            #legend.position = "bottom",
            #legend.direction = "horizontal",
            legend.key.size= unit(0.25, "cm"), # unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            #legend.title = element_text(face="italic"),
            legend.title = element_blank(),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background = element_blank(),
            #strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0", size = 0.5),
            
            strip.text = element_text(face="bold", size = 12)
    ))
  
}


#########################################################
# Plotting the trajectory:

#' Prepares the data for visualizing linear or branching process
#' 
#' The function takes the matrix of multiple Hamiltonian path solutions, which is automatically saved to the global environment as 'norm_res_new' after running the redpath function.
#' It calculates the transition matrix and prepares the PCA coordinates for visualizing the linear / branching process.
#' 
#' 
#'@param all_hmt_soln The Hamiltonian path solutions of your data, stored as 'norm_res_new' after running the redpath() function.
#'@param pseudotime The inferred pseudo development time
#'@param labels The cell type / time-point labels.
#'@export
#'@keywords get_branch_viz_data
#'@section Biological Analysis:
#'@examples
#'redpath_pseudotime <- redpath(data, threadnum = 4) 
#'# Note that norm_res_new is stored as a global variable and it will be overwritten each time when you perform redpath. 
#'all_hamiltonian_solutions <- norm_res_new
#'branching_visualization_data <- get_branch_viz_data(all_hamiltonian_solutions, pseudotime = redpath_pseudotime, labels = NULL)
#'# Plot to visualize the trajectory
#'traj_plots <- plot_trajectory(branching_visualization_data, redpath_pseudotime, path = F, rev = F)
#'library(gridExtra)
#'grid.arrange(traj_plots[[1]], traj_plots[[2]], layout = matrix(c(1:2), 1, 2)) 
get_branch_viz_data <- function(all_hmt_soln, pseudotime, labels){
  processed_branch_data <- process_branch(all_hmt_soln, multiple = T, n_cells = length(pseudotime))
  normalized_data <- processed_branch_data * pseudotime
  normalized_data <- t(t(normalized_data) * pseudotime)
  
  prepared_data <- get_pca_data(normalized_data, pseudotime, labels = labels)
  return(prepared_data)
}

library(ggplot2)
library(ggalt)

#' Plotting function for visualization of the trajectory
#' 
#' This function takes the prepared data from the 'get_branch_viz_data' function. 
#' It produces two PCA plots, showing the pseudotime information and cell type information separately.
#' 
#' 
#'@param prepared_data The data returned from the get_branch_viz_data() function.
#'@param pseudotime The inferred pseudo development time
#'@param path Defaults to F. Decides whether to apply the geom_path function on the plots.
#'@param rev Defaults to F. Sets the forward or reverse direction of the pseudotime.
#'@export
#'@keywords plot_trajectory
#'@section Biological Analysis:
#'@examples
#'redpath_pseudotime <- redpath(data, threadnum = 4) 
#'# Note that norm_res_new is stored as a global variable and it will be overwritten each time when you perform redpath. 
#'all_hamiltonian_solutions <- norm_res_new
#'branching_visualization_data <- get_branch_viz_data(all_hamiltonian_solutions, pseudotime = redpath_pseudotime, labels = NULL)
#'# Plot to visualize the trajectory
#'traj_plots <- plot_trajectory(branching_visualization_data, redpath_pseudotime, path = F, rev = F)
#'
#'library(gridExtra)
#'grid.arrange(traj_plots[[1]], traj_plots[[2]], layout = matrix(c(1:2), 1, 2)) 
plot_trajectory <- function(prepared_data, pseudotime, path = F, rev = F){
  mid <- mean(pseudotime)
  # 
  # circle_section <- kmeans(pseudotime, k)
  # 
  # t1 <- prepared_data[which(circle_section$cluster==1), ]
  # t2 <- prepared_data[which(circle_section$cluster==2), ]
  # t3 <- prepared_data[which(circle_section$cluster==3), ]
  # t4 <- prepared_data[which(circle_section$cluster==4), ]
  # t5 <- prepared_data[which(circle_section$cluster==5), ]
  # 
  # g1 <- ggplot(prepared_data, aes(X, Y, color = as.factor(labels)))+geom_point(size = 3) + 
  #   scale_color_manual(name = "", values = as.character(redpath_colorset[1:length(unique(prepared_data$labels))])) + redpath_theme()+
  #   geom_encircle(aes(x = X, y = Y),
  #                 data = t1,
  #                 color = "red", size = 1, expand = 0.01)+
  #   geom_encircle(aes(x = X, y = Y),
  #                 data = t2,
  #                 color = "red", size = 1, expand = 0.01)+
  #   geom_encircle(aes(x = X, y = Y),
  #                 data = t3,
  #                 color = "red", size = 1, expand = 0.01)+
  #   geom_encircle(aes(x = X, y = Y),
  #                 data = t4,
  #                 color = "red", size = 1, expand = 0.01)+
  #   geom_encircle(aes(x = X, y = Y),
  #                 data = t5,
  #                 color = "red", size = 1, expand = 0.01)
  if(rev == T){
    prepared_data$PT = 1 - prepared_data$PT
  }
  g1 <- ggplot(prepared_data, aes(X, Y, color = PT))+geom_point(size = 5) + 
    theme_redpath() + scale_color_gradient2(midpoint=mid, low="blue", mid="white", high = "red")+
    #scale_shape_manual(values = c(0, 6, 3, 2, 25))+
    xlab("Component 1") + ylab("Component 2") + labs(colour = "Pseudotime", title = "Pseudotime trend on PCA plot")
  if(path == F){
    if(length(unique(prepared_data$labels)) < 10){
      g2 <- ggplot(prepared_data, aes(X, Y, color = as.factor(labels)))+geom_point(size = 5) + 
        scale_color_manual(name = "Cell Type", values = as.character(redpath_colorset)) + 
        theme_redpath()+
        #scale_shape_manual(values = c(0, 6, 3, 2, 25))+
        xlab("Component 1") + ylab("Component 2") + labs(colour = "Cell_Type", title = "Cell Type PCA plot")
    }else{
      g2 <- ggplot(prepared_data, aes(X, Y, color = as.factor(labels)))+geom_point(size = 5) + 
        #scale_color_manual(name = "Cell Type", values = as.character(redpath_colorset)) + 
        theme_redpath()+
        #scale_shape_manual(values = c(0, 6, 3, 2, 25))+
        xlab("Component 1") + ylab("Component 2") + labs(colour = "Cell_Type", title = "Cell Type PCA plot")
    }
    
  }else{
    if(length(unique(prepared_data$labels)) < 10){
      g2 <- ggplot(prepared_data, aes(X, Y, color = as.factor(labels)))+geom_point(size = 5) + 
        scale_color_manual(name = "Cell Type", values = as.character(redpath_colorset)) + 
        theme_redpath()+geom_path(aes(color = labels))+
        #scale_shape_manual(values = c(0, 6, 3, 2, 25))+
        xlab("Component 1") + ylab("Component 2") + labs(colour = "Cell_Type", title = "Cell Type PCA plot")
    }else{
      g2 <- ggplot(prepared_data, aes(X, Y, color = as.factor(labels)))+geom_point(size = 5) + 
        #scale_color_manual(name = "Cell Type", values = as.character(redpath_colorset)) + 
        theme_redpath()+geom_path(aes(color = labels))+
        #scale_shape_manual(values = c(0, 6, 3, 2, 25))+
        xlab("Component 1") + ylab("Component 2") + labs(colour = "Cell_Type", title = "Cell Type PCA plot")
    }
    
  }
  g_list <- list()
  g_list[[1]] <- g1
  g_list[[2]] <- g2
  
  return(g_list)
}
#########################################################
library(combinat)


get_pca_data <- function(processed_branch_distance_matrix, pseudotime, labels = NULL){
  
  #hclust_results <- cutree(hclust(as.dist(processed_branch_distance_matrix)), k)
  #center_points <- aggregate(processed_branch_distance_matrix, list(hclust_results), mean)
  #center_points <- center_points[,-1]
  
  #ew_matrix <- rbind(processed_branch_distance_matrix, center_points)
  pca_data <- prcomp(processed_branch_distance_matrix, center = T, scale. = T)
  #k1 <- kmeans(pca_data$rotation[, c(1:100)], k, iter.max = 500, nstart = 10)
  
  if(is.null(labels)){
    tmp_df <- data.frame(X = pca_data$rotation[,1], Y = pca_data$rotation[,2], PT = pseudotime)
  }else{
    tmp_df <- data.frame(X = pca_data$rotation[,1], Y = pca_data$rotation[,2], labels = labels, PT = pseudotime)
  }
  return(tmp_df)
}


process_branch <- function(all_hmt_soln, multiple = F, n_cells){
  if(multiple == T){
    total_soln <- nrow(all_hmt_soln)
    #n_cells <- ncol(all_hmt_soln)
    
    trans_matrix <- matrix(0, n_cells, n_cells)
    for(i in c(1:total_soln)){
      ith_soln <- all_hmt_soln[i, ]
      unique_ <- unique(ith_soln)
      unique_order <- order(unique_)
      
      uniq_len <- length(unique_)
      for(j in c(1:(uniq_len-1))){
        start_pts <- which(ith_soln == unique_[unique_order[j]])
        end_pts <- which(ith_soln == unique_[unique_order[j+1]])
        trans_matrix[start_pts, end_pts] <- trans_matrix[start_pts, end_pts] + (1 / uniq_len) #(unique_[unique_order[j+1]] * unique_[unique_order[j]])#
      }
    
    }
    return(trans_matrix)
  }else{
    ith_soln <- all_hmt_soln
    unique_ <- unique(ith_soln)
    unique_order <- order(unique_)
    
    trans_matrix <- matrix(0, n_cells, n_cells)
    
    uniq_len <- length(unique_)
    for(j in c(1:(uniq_len-1))){
      start_pts <- which(ith_soln == unique_[unique_order[j]])
      end_pts <- which(ith_soln == unique_[unique_order[j+1]])
      trans_matrix[start_pts, end_pts] <- trans_matrix[start_pts, end_pts] + 1 / uniq_len
    }
    return(trans_matrix)
    
  }

}
