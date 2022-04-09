library(monocle3)
######monocle3###
Thy_c<-readRDS("Thy.rds")
Thy_cs<-subset(Thy_c, idents = c(0:14))

expression_matrix <- Thy_cs@assays$RNA@counts
cell_metadata <- Thy_cs@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(Thy_cs@assays$RNA), row.names = rownames(Thy_cs@assays$RNA))
# Seurat-derived CDS
my.cds3 <- new_cell_data_set(expression_matrix,
                             cell_metadata = cell_metadata,
                             gene_metadata = gene_annotation)

# Transfer Seurat embeddings
reducedDim(my.cds3, type = "PCA") <- Thy_cs@reductions$pca@cell.embeddings 
my.cds3@preprocess_aux$prop_var_expl <- Thy_cs@reductions$pca@stdev

my.cds3@int_colData@listData$reducedDims$UMAP <- Thy_cs@reductions$umap@cell.embeddings
plot_cells(my.cds3)

my.cds3 <- cluster_cells(my.cds3, reduction_method = "UMAP")

# transfer cluster
Thy_cs[[sprintf("Cluster")]] <- Idents(object = Thy_cs)
list_cluster <- Thy_cs@meta.data[[sprintf("Cluster")]]
names(list_cluster) <- Thy_cs@assays[["RNA"]]@data@Dimnames[[2]]
my.cds3@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster

plot_cells(my.cds3,group_label_size = 3.5, cell_size = 0.6)
plot_cells(my.cds3, color_cells_by = "partition")

my.cds3@clusters$UMAP$partitions[my.cds3@clusters$UMAP$partitions == "2"] <- "1"
my.cds3@clusters$UMAP$partitions[my.cds3@clusters$UMAP$partitions == "3"] <- "1"
my.cds3 <- learn_graph(my.cds3, use_partition = TRUE)
my.cds3 <- choose_cells(my.cds3)
colData(my.cds3)$cell_type <- as.character(clusters(my.cds3))
colData(my.cds3)$cell_type=dplyr::recode(colData(my.cds3)$cell_type,
                                         "0"="0",
                                         "1"="1",
                                         "2"="2",
                                         "3"="3",
                                         "4"="4",
                                         "5"="5",
                                         "6"="6",
                                         "7"="7",
                                         "8"="8",
                                         "9"="9",
                                         "10"="10",
                                         "11"="11",
                                         "12"="12",
                                         "13"="13",
                                         "14"="14")

# a helper function to identify the root principal points:
get_earliest_principal_node <- function(my.cds3,time_bin="0"){
  cell_ids <- colData(my.cds3)[,"cell_type"] == time_bin
  closest_vertex <-
    my.cds3@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(my.cds3), ])
  root_pr_nodes <-
    igraph::V(principal_graph(my.cds3)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

n=get_earliest_principal_node(my.cds3) 

my.cds3 <- order_cells(my.cds3, root_pr_nodes=c(n))
plot_cells(my.cds3,
           color_cells_by = "pseudotime",
           label_groups_by_cluster = TRUE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           cell_size = 0.5,
           group_label_size = 4,
           trajectory_graph_segment_size = 0.8)

plot_genes_in_pseudotime(my.cds3[genes,colData(my.cds3)$cell_type %in% c(6:15)], color_cells_by = "cell_type")


##########Finding genes that change as a function of pseudotime########
cds_pr_test_res2 <- graph_test(my.cds3, neighbor_graph="principal_graph")
pr_deg_ids <- row.names(subset(my.cds3, q_value < 0.05))


gene_module <- find_gene_modules(my.cds3[pr_deg_ids,], resolution = c(10^seq(-6,-1)))


cell_group <- tibble::tibble(cell=row.names(colData(my.cds3)), 
                             cell_group=colData(my.cds3)$cell_type)
agg_mat <- aggregate_gene_expression(my.cds3, gene_module, cell_group)
row.names(agg_mat) <- stringr::str_c("Module", row.names(agg_mat))
pheatmap::pheatmap(agg_mat,
                   scale="column", clustering_method="ward.D2")


