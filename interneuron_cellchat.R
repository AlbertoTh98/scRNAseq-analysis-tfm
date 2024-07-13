rm(list=ls()) #  --> clean your environment

library(CellChat)
library(reticulate)
library(anndata)
library(SingleCellExperiment)
library(zellkonverter)
library(Matrix)
library(NMF)
library(ggalluvial)

py_config()
py_module_available('matplotlib')
scipy_sparse = import('scipy.sparse')
numba = import('numba')
sc = import('scanpy')
scipy <- import("scipy")
scipy_sparse <- scipy$sparse
np = import('numpy')

# Análisis de poblaciones de interneuronas en diferentes condiciones


# 1. Preparar los sce
adata_control <- read_h5ad("C:/Users/Usuario/OneDrive/Desktop/interneurons/representacion_resultados/subset_control.h5ad")
adata_injection <- read_h5ad("C:/Users/Usuario/OneDrive/Desktop/interneurons/representacion_resultados/subset_injection.h5ad")

# control
counts_control <- as(adata_control$X, "CsparseMatrix")
counts_control = t(counts_control)
rownames(counts_control) <- as.character(py_to_r(adata_control$var_names))
colnames(counts_control) <- as.character(py_to_r(adata_control$obs_names))
col_data_control <- as.data.frame(py_to_r(adata_control$obs))
row_data_control <- as.data.frame(py_to_r(adata_control$var))

# En principio no es necesario pca ni umap coords
sce_control <- SingleCellExperiment(assays = list(counts = counts_control), 
                            rowData = row_data_control,
                            colData = col_data_control)
rownames(sce_control) = rowData(sce_control)$var_names
rm(counts_control, adata_control, col_data_control, row_data_control)

# injection
counts_injection <- as(adata_injection$X, "CsparseMatrix")
counts_injection = t(counts_injection)
rownames(counts_injection) <- as.character(py_to_r(adata_injection$var_names))
colnames(counts_injection) <- as.character(py_to_r(adata_injection$obs_names))
col_data_injection <- as.data.frame(py_to_r(adata_injection$obs))
row_data_injection <- as.data.frame(py_to_r(adata_injection$var))

# En principio no es necesario pca ni umap coords
sce_injection <- SingleCellExperiment(assays = list(counts = counts_injection), 
                            rowData = row_data_injection,
                            colData = col_data_injection)
rownames(sce_injection) = rowData(sce_injection)$var_names
rm(counts_injection, adata_injection, col_data_injection, row_data_injection)




# 2. Input de cellchat ############################################################################################

# objeto singlecellexperiment injection 
data_injection.input <- SingleCellExperiment::counts(sce_injection) # normalized data matrix
meta_injection <- as.data.frame(SingleCellExperiment::colData(sce_injection)) # extract a dataframe of the cell labels
meta_injection$labels <- meta_injection["leiden_r0.01"]
# Convertir a caracteres 
meta_injection$obs_names <- as.character(meta_injection$obs_names)
meta_injection$obs_names <- make.unique(meta_injection$obs_names)
rownames(meta_injection) <- meta_injection$obs_names
colnames(data_injection.input) = rownames(meta_injection)

# Input de cellchat: objeto singlecellexperiment control
data_control.input <- SingleCellExperiment::counts(sce_control) # normalized data matrix
meta_control <- as.data.frame(SingleCellExperiment::colData(sce_control)) # extract a dataframe of the cell labels
meta_control$labels <- meta_control["leiden_r0.01"]
# Convertir a caracteres 
meta_control$obs_names <- as.character(meta_control$obs_names)
meta_control$obs_names <- make.unique(meta_control$obs_names)
rownames(meta_control) <- meta_control$obs_names
colnames(data_control.input) = rownames(meta_control)



# 3. create cellchat object for injection
cellchat_injection <- createCellChat(object = data_injection.input, meta = meta_injection, group.by = "leiden_r0.01_sub_PV2_SST")
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB. We do not suggest to use it in this way because CellChatDB v2 includes "Non-protein Signaling" (i.e., metabolic and synaptic signaling). 
# set the used database in the object
cellchat_injection@DB <- CellChatDB.use

# create cellchat object for control
cellchat_control <- createCellChat(object = data_control.input, meta = meta_control, group.by = "leiden_r0.01_sub_PV2_SST")
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB. We do not suggest to use it in this way because CellChatDB v2 includes "Non-protein Signaling" (i.e., metabolic and synaptic signaling). 
# set the used database in the object
cellchat_control@DB <- CellChatDB.use



# subset the expression data of signaling genes for saving computation cost
cellchat_injection <- subsetData(cellchat_injection) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat_injection <- identifyOverExpressedGenes(cellchat_injection)
cellchat_injection <- identifyOverExpressedInteractions(cellchat_injection)
# communication network
cellchat_injection <- computeCommunProb(cellchat_injection, type = "triMean")
cellchat_injection <- filterCommunication(cellchat_injection, min.cells = 10)
df_injection.net <- subsetCommunication(cellchat_injection)
# Infer cell-cell communication at a signaling pathway level.
cellchat_injection <- computeCommunProbPathway(cellchat_injection)
# Calculate aggregated cell-cell communication network.
cellchat_injection <- aggregateNet(cellchat_injection)


# subset the expression data CONTROL of signaling genes for saving computation cost
cellchat_control <- subsetData(cellchat_control) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat_control <- identifyOverExpressedGenes(cellchat_control)
cellchat_control <- identifyOverExpressedInteractions(cellchat_control)
# communication network
cellchat_control <- computeCommunProb(cellchat_control, type = "triMean")
cellchat_control <- filterCommunication(cellchat_control, min.cells = 10)
df_control.net <- subsetCommunication(cellchat_control)
# Infer cell-cell communication at a signaling pathway level.
cellchat_control <- computeCommunProbPathway(cellchat_control)
# Calculate aggregated cell-cell communication network.
cellchat_control <- aggregateNet(cellchat_control)





################################################ ANALISIS
object.list <- list(INJECTION = cellchat_injection, CONTROL = cellchat_control)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

#When comparing cell-cell communication among multiple biological conditions, it can
#answer the following biological questions:
#• Whether the cell-cell communication is enhanced or not;
#• The interaction between which cell types is significantly changed;
#• How the major sources and targets change from one condition to another.



# 17. Compare the total number of interactions and interaction strength. To answer the
# question on whether the cell-cell communication is enhanced or not.
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2


# 18. Compare the number of interactions and interaction strength among different cell
# populations. To identify the interaction between cell populations showing significant
# changes.

# (A) Circle plot showing differential number of interactions or interaction strength
# among different cell populations across two datasets
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

# (B) Heatmap showing differential number of interactions or interaction strength
# among different cell populations across two datasets
gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2

# (C) Circle plot showing the number of interactions or interaction strength among
# different cell populations across multiple datasets
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T,
                   label.edge= F, edge.weight.max = weight.max[2], edge.width.max =
                     12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}



# (D) Circle plot showing the differential number of interactions or interaction
# strength among coarse cell types
# Here, CellChat categorize the cell populations into three cell types, and then re-merge 
# the list of CellChat objects.
group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4))
group.cellType <- factor(group.cellType, levels = c("FIB", "DC", "TC"))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

# Show the number of interactions or interaction strength between any two cell types in each dataset.
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "
net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
# Similarly, CellChat can also show the differential number of interactions or 
# interaction strength between any two cell types using circle plot.
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "
count.merged", label.edge = T)



# 19. Compare major sources and targets in 2D space.

# Calcular la centralidad de la red para cada objeto CellChat en object.list
for (i in seq_along(object.list)) {
  object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]])
}

# (A) Identify cell populations with significant changes in sending or receiving signals
num.link <- sapply(object.list, function(x) {rowSums(x@net$count)
  + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]],
                                               title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)

# (B) Identify the signaling changes of specific cell populations
gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Cck", signaling.exclude = "MIF")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Cck/Vip", signaling.exclude = c("MIF"))
patchwork::wrap_plots(plots = list(gg1,gg2))






# 20. Identify signaling networks with larger (or smaller) difference as well as signaling 
# groups based on their functional/structure similarity.
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)
# Compute and visualize the pathway distance in the learned joint manifold
rankSimilarity(cellchat, type = "functional")



# 21. Identify altered signaling with distinct interaction strength.

# (A) Compare the overall information flow of each signaling pathway.
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

# (B) Compare outgoing (or incoming) signaling patterns associated with each cell population
library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))




# 22. Identify dysfunctional signaling by comparing the communication probabilities.
# Compare the communication probabilities from certain cell groups to other cell groups
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2), angle.x = 45)
# Identify the up-regulated (increased) and down-regulated (decreased) signaling ligand-receptor pairs in one dataset compared to the other dataset.
gg1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling", angle.x = 45, remove.isolate = T)
gg2 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling", angle.x = 45, remove.isolate = T)
gg1 + gg2


# 23. Identify dysfunctional signaling by using differential expression analysis.
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "INJECTION"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, 
                                       thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "INJECTION",ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "INJECTION",ligand.logFC = -0.1, receptor.logFC = -0.1)
# do further deconvolution to obtain the individual signaling genes
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
# Users can also find all the significant outgoing/incoming/both signaling according to the customized DEG features and cell groups of interest
df <- findEnrichedSignaling(object.list[[2]], features = c("CCL19", "CXCL12"), idents = c("Inflam. FIB", "COL11A1+ FIB"), pattern ="outgoing")


# 24. Visualize the identified up-regulated and down-regulated signaling ligand-receptor
# pairs using bubble plot (option A), chord diagram (option B) or wordcloud (option C).

# (A) Bubble plot
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2), angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2), angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg1 + gg2

# (B) Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = 4, targets.use = c(5:11), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = 4, targets.use = c(5:11), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))

# (C) Wordcloud plot


library("wordcloud")

# visualize the enriched ligands in the first condition
computeEnrichmentScore(net.down, species = 'mouse')
# visualize the enriched ligands in the second condition
computeEnrichmentScore(net.up, species = 'mouse')

########## Visually compare inferred cell-cell communication networks


# 25. Visualize inferred cell-cell communication networks.

# (A) Circle plot.
pathways.show <- c("CXCL")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "c
ircle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name =
                        paste(pathways.show, names(object.list)[i]))
}


# (B) Heatmap plot
pathways.show <- c("CXCL")
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))



# 26. Visualize gene expression distribution. 
# Plot the gene expression distribution of signaling genes related to L-R pairs or 
# signaling pathway using a Seurat wrapper function `plotGeneExpression`.
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("NL", "LS")) # set factor level
plotGeneExpression(cellchat, signaling = "CXCL", split.by = "datasets", colors.ggplot = T)



save(object.list, file = "cellchat_object.list_humanSkin_NL_LS.RData")
save(cellchat, file = "cellchat_merged_humanSkin_NL_LS.RData")















