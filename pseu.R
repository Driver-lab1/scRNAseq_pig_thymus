
#########SCORPIUS
library(SCORPIUS)

topg = VariableFeatures(object = NKT_n)
group <- Idents(NKT_n)
expression_a <- t(as.matrix(NKT_n[topg,]@assays$RNA@data))
head(expression_a)

space_a <- reduce_dimensionality(expression_a, dist = "spearman", ndim = 3)
draw_trajectory_plot(space_a, progression_group = group, contour = TRUE,point_size = -1)
traj_a <- infer_trajectory(space_a)
draw_trajectory_plot(
  space_a, 
  progression_group = group,
  path = traj_a$path,
  contour = T,
  point_size = 0.7,
  progression_group_palette = palette
)


gimp_a <- gene_importances(expression_a, traj_a$time, num_permutations = 0, num_threads = 8)
gene_sel_a <- gimp_a[1:50,]
expr_sel_a <- expression_a[,gene_sel_a$gene]
draw_trajectory_heatmap(expr_sel_a, traj_a$time, group, show_labels_row = T)

modules_a <- extract_modules(scale_quantile(expr_sel_a), traj_a$time, verbose = FALSE)
draw_trajectory_heatmap(expr_sel_a, traj_a$time, group, modules, fontface="bold",show_labels_row = T)

########slingshot

library(slingshot)
library(Polychrome)
library(ggbeeswarm)
library(ggthemes)
library(SingleCellExperiment)

nkt<-as.SingleCellExperiment(nkt)
table(nkt$ident)
sce <- slingshot(nkt, clusterLabels = nkt$ident, reducedDim = "PCA",
                 allow.breaks = FALSE)
# get the lineages:
lnes <- getLineages(reducedDim(sce,"PCA"), sce$ident)
lnes@lineages

slingshot_df <- data.frame(colData(sce))

nkt_cols <- c('1'='#E41A1C','2'='#377EB8','3'='#4DAF4A','4'='#984EA3','5'='#FF7F00',
              '6'='#FFFF33','7'='#A65628','8'='#F781BF','9'='#999999')

plot(reducedDims(sce)$PCA, col = nkt_cols[as.character(sce$ident)],
     pch=16, 
     asp = 1,cex=1)
legend("bottomright",legend = names(nkt_cols[levels(sce$ident)]),  
       fill = nkt_cols[levels(sce$ident)], border = F,cex=0.6)
lines(SlingshotDataSet(lnes), lwd=2, type = 'lineages', col = c("black"), cex=2)

