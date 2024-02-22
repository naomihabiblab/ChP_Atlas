library(dplyr)
library(reshape2)
library(tidyverse)
library(pheatmap)
library(cowplot)
library(dplyr)
library(plotly)
library(RColorBrewer)
library(tibble)

##############################################################################################
##                                     Utill Functions                                      ##
##############################################################################################
SymmetrizeGraph <- function(graph) {
    d <- diag(graph)
    graph <- graph + t(graph)
    diag(graph) <- d
    return(graph)
}


CreateGraph <- function(edges, sym=F, reps=10000, verbose=T) {
  graph <- table(edges$lcell, edges$rcell, dnn=c("Ligand Cell", "Receptor Cell"))
  if(sym) graph <- SymmetrizeGraph(graph)
  
  pvals <- apply(simplify2array(lapply(1:reps, function(.)  {
    e <- edges
    # e$rcell <- sample(e$rcell)
    # e$lcell <- sample(e$lcell)
    e$rcell <- sample(unique(e$rcell), replace = T, size = nrow(e))
    e$lcell <- sample(unique(e$lcell), replace = T, size = nrow(e))
    p.graph <- table(e$lcell, e$rcell)
    if(sym) p.graph <- SymmetrizeGraph(p.graph)
    return(p.graph >=  graph)
  })), 1:2, mean)
  qvals <- matrix(p.adjust(pvals, method = "BY"), nrow = nrow(pvals), dimnames = list(`Ligand Cell`=rownames(graph), `Receptor Cell`=colnames(graph)))
  
  if(verbose) {
    message("Created Graph:");   print(graph)
    message("\nCalculated p-values:");   print(pvals)
    message("\nCalculated adjusted p-values:");   print(qvals)
  }
  
  return(list(edges=edges, graph=graph, pvals=pvals, qvals=qvals,
              disp=matrix(cut(qvals, c(-.1, 0.001, 0.01, 0.5, Inf), c("***", "**", "*", "")),nrow = nrow(graph))))
}


CreateDifferentialGraph <- function(edges.a, edges.b, sym=F, reps=10000, verbose=T) {
  edges <- setdiff(edges.a, edges.b)
  graph <- table(edges$lcell, edges$rcell, dnn=c("Ligand Cell", "Receptor Cell"))
  if(sym) graph <- SymmetrizeGraph(graph)
  
  pvals <- apply(simplify2array(lapply(1:reps, function(.)  {
    e <- edges.b
    e$rcell <- sample(e$rcell)
    e$lcell <- sample(e$lcell)
    e <- setdiff(edges.a, e)
    p.graph <- table(e$lcell, e$rcell)
    if(sym) p.graph <- SymmetrizeGraph(p.graph)
    return(p.graph >=  graph)
  })), 1:2, mean)
  
  qvals <- matrix(p.adjust(pvals, method = "BY"), nrow = nrow(pvals), dimnames = list(`Ligand Cell`=rownames(graph), `Receptor Cell`=colnames(graph)))
  
  if(verbose) {
    message("Created Graph:");   print(graph)
    message("\nCalculated p-values:");   print(pvals)
    message("\nCalculated adjusted p-values:");   print(qvals)
  }
  
  return(list(edges=edges, graph=graph, pvals=pvals, qvals=qvals,
              disp=matrix(cut(qvals, c(-.1, 0.001, 0.01, 0.5, Inf), c("***", "**", "*", "")),nrow = nrow(graph))))
}


.SankeyPlot <- function(nodes, links, title=NULL) {
  library(networkD3)
  library(htmlwidgets)
  library(htmltools)
  colors <- nord::nord_palettes$aurora
  sn <- sankeyNetwork(links, nodes, 
                      Source="source", Target="target", Value="perc",
                      NodeID="name", 
                      colourScale = paste0('d3.scaleOrdinal() .domain([0,1,2,3,4]) .range(["', paste0(colors, collapse = '","'), '"])'),
                      LinkGroup = "group", NodeGroup = "group",
                      iterations = 0,
                      fontSize = 12,
                      margin = list("left"=100))
  if(!is.null(title))
    sn <- htmlwidgets::prependContent(sn, htmltools::tags$h3(title))
  onRender(
    sn,  'function(el,x){
      // select all our node text
      var node_text = d3.select(el)
        .selectAll(".node text")
        //and make them match
        //https://github.com/christophergandrud/networkD3/blob/master/inst/htmlwidgets/sankeyNetwork.js#L180-L181
        .attr("x", 6 + x.options.nodeWidth)
        .attr("text-anchor", "start");
    }')
  return(sn)
}

TwoWay.SankeyPlot <- function(obj, to.return=F, title=NULL) {
  links <- melt(obj$graph)
  colnames(links) <- c("source","target","value")
  links <- data.frame(links %>% arrange(target, source) %>% group_by(target) %>% mutate(perc=round(100*value/sum(value), 2)))
  links$target <- paste(links$source, links$target, sep="_")
  rownames(links) <- links$target
  
  nodes <- data.frame(name=c(colnames(obj$graph), links$target))
  rownames(nodes) <- nodes$name
  nodes$id <- 1:nrow(nodes) -1
  nodes$group <- as.character(nodes[gsub("_.*$", "", nodes$name), "id"])
  nodes[-(1:5),]$name <- paste0(links[nodes[-(1:5), "name"], "perc"], "%")
  
  links$group  <- as.character(nodes[links$source, "group"])
  links$source <- nodes[links$source, "id"]
  links$target <- nodes[links$target, "id"]
  
  if(to.return)
    return(list(nodes=nodes, links=links)) 
  .SankeyPlot(nodes, links, title)
}

ThreeWay.SankeyPlot <- function(obj1, obj2, to.return=F, title=NULL) {
  sankey1 <- TwoWay.SankeyPlot(obj1, T)
  sankey2 <- TwoWay.SankeyPlot(obj2, T)
  
  sankey1$nodes <- sankey1$nodes[-(1:5),]
  sankey1$nodes$id <- max(sankey1$nodes$id) + (1:nrow(sankey1$nodes))
  sankey1$links$target <- max(sankey1$links$target) + (1:nrow(sankey1$links))
  colnames(sankey1$links)[1:2] <- colnames(sankey1$links)[2:1]
  
  nodes <- rbind(sankey1$nodes, sankey2$nodes) %>% arrange(id)
  links <- rbind(sankey1$links, sankey2$links) %>% arrange(source)
  
  if(to.return)
    return(list(nodes=nodes, links=links))  
  .SankeyPlot(nodes, links, title)
}
##############################################################################################





##############################################################################################
##                                   Interaction Graph                                      ##
##############################################################################################
data.sets <- c("Aged_TCP", "Aged_DCP", "Aged_HCP", "Adult_TCP","Adult_DCP", "Adult_HCP",  "Embryo_TCP", "Embryo_DCP",
               "Embryo_HCP", "SCEmbryo_TCP", "SCEmbryo_DCP", "SCEmbryo_HCP")

results <- lapply(data.sets, function(data.set) {
  message(paste0("\n#############################################################\n","\t\t\t",data.set,
                 "\n#############################################################\n"))
  message(paste(Sys.time(),"\tBeginning", data.set))

  ix <- read.csv(file.path(data.set, "out/significant_means.txt"), header = T, sep="\t", stringsAsFactors = F)
  ix$receptor_a <- ix$receptor_a == "True"; 
  ix$receptor_b <- ix$receptor_b == "True";

  ix[ix$gene_a == "", "gene_a"] <- gsub("complex:","", ix[ix$gene_a == "","partner_a"])
  ix[ix$gene_b == "", "gene_b"] <- gsub("complex:","", ix[ix$gene_b == "","partner_b"])
  ix <- ix[,-c(1,2,3,4,7,11,12)]
  
  ix <- ix[xor(ix$receptor_a, ix$receptor_b) | ix$annotation_strategy == "curated",-5]
  ix <- melt(ix, id.vars = c("gene_a","gene_b","receptor_a","receptor_b"), na.rm = T)
  ix$variable <- as.character(ix$variable)
  ix <- ix %>% separate(variable, c("ident_a", "ident_b"), "\\.")
  
  edges <- as.data.frame(mapply(c, 
              ix[!ix$receptor_a & ix$receptor_b, c("ident_a","ident_b","gene_a","gene_b")],
              ix[ix$receptor_a & !ix$receptor_b, c("ident_b","ident_a","gene_b","gene_a")],
              ix[ix$receptor_a == ix$receptor_b, c("ident_a","ident_b","gene_a","gene_b")]))
  
  colnames(edges) <- c("lcell", "rcell", "lig", "rec")
  edges$lweight <- 1; edges$rweight <- 1

  edges$lcell <- plyr::mapvalues(edges$lcell, c("Endothelial","Epithelial", "Fibroblast"), c("Endo","Epi","Mes"))
  edges$rcell <- plyr::mapvalues(edges$rcell, c("Endothelial","Epithelial", "Fibroblast"), c("Endo","Epi","Mes"))
  
  edges$lcell <- factor(edges$lcell, levels = c("Epi", "Mes", "Endo", "NeuroGlia", "Immune"))
  edges$rcell <- factor(edges$rcell, levels = c("Epi", "Mes", "Endo", "NeuroGlia", "Immune"))
  
  message(paste(Sys.time(),"\tCreating graph and shuffled graphs"))
  obj <- CreateGraph(edges)
  obj$setting <- gsub("_", " ", data.set)
  
  
  setwd("graphs/")
  TwoWay.SankeyPlot(obj, title = obj$setting) %>% saveNetwork(paste0("single.", obj$setting, ".sankey.html"))
  webshot(paste0("single.", obj$setting, ".sankey.html"), paste0("single.", obj$setting, ".sankey.png"))
  setwd("../")
  
  return(obj)
})
names(results) <- data.sets
##############################################################################################





##############################################################################################
##                               Interaction Graph Heatmaps                                 ##
##############################################################################################
scale.min <- 0 #Reduce(min, lapply(results, function(r) min(r$graph)))
scale.max <- Reduce(max, lapply(results, function(r) max(r$graph)))
breaksList <- seq(scale.min, scale.max, length.out=10)
cols <- (colorRampPalette(brewer.pal(length(breaksList), "Reds"))(length(breaksList)))

ggsave(width=9, height = 11, dpi=400, filename = "graphs/heatmaps.pdf", plot=
  plot_grid(plotlist = lapply(seq_along(results), function(i) {
    res <- results[[i]]
    r <-floor((i-1)/3) + 1; c <- (i-1) %% 3 + 1
    
    # scale.min <- min(res$graph)
    # scale.max <- max(res$graph)
    # breaksList <- seq(scale.min, scale.max, length.out=9)
    # cols <- (colorRampPalette(brewer.pal(length(breaksList), "Reds"))(length(breaksList)))
    p <- pheatmap(res$graph, display_numbers = res$disp, main = res$setting, 
             cluster_cols = F, cluster_rows = F, 
             legend = (r==2 && c==3), 
             show_rownames = c==3, show_colnames = r==4,
             cellwidth = 25, cellheight = 25, fontsize_number = 14,
             breaks = breaksList, color=cols)[[4]]
    p}),
    ncol = 3, nrow = 4,
    rel_heights = c(1,1,1,1.2),
    rel_widths = c(1,1,1,1.2)))

rm(scale.min, scale.max, breaksList, cols)
##############################################################################################





##############################################################################################
##                               Interaction Graph Bar-Plots                               ##
##############################################################################################
df <- do.call(rbind, lapply(seq_along(results), function(i) 
    data.frame(case=names(results)[i], melt(results[[i]]$graph / sum(results[[i]]$graph) )))) %>%  
  separate(case, c("Age", "Region"), "_")
df$Age <- factor(df$Age, levels = c("SCEmbryo","Embryo","Adult","Aged"))
df$Region <- factor(df$Region, c("TCP","DCP","HCP"))
df <- df[df$Age != "SCEmbryo",]


ggsave(width=10, height=10, dpi=300, filename="graphs/interactions.bar.pdf", plot=plot_grid(plotlist=
  lapply(levels(df$Region), function(r) 
    ggplot(df[df$Region == r,] %>% group_by(Age, Region, Ligand.Cell) %>% mutate(value=sum(value)), aes(x=Ligand.Cell, y=value, fill=Age)) + 
      geom_bar(stat = "identity", position = position_dodge()) +
      #facet_wrap("Age", ncol = 1, scales = "free_y") +
      theme_classic() + 
      theme(legend.position = "bottom", legend.justification = c(0, 0)) + 
      labs(x="", y="Number of Interactions", fill="") + 
      ggtitle(paste0(r," - Interactions Across Age Categories")) + 
      scale_fill_brewer(palette="Dark2")), ncol = 1))
rm(df)
##############################################################################################




##############################################################################################
##                                Differential Interaction Graphs                           ##
##############################################################################################
pairs <- data.frame(t(combn(data.sets, 2)))  %>% 
  separate(X1, c("X1_a", "X1_b"), "_", remove = F)  %>% 
  separate(X2, c("X2_a", "X2_b"), "_", remove = F) %>% 
  filter(X1_a == X2_a | X1_b == X2_b) %>% 
  mutate(title=ifelse(X1_a == X2_a, paste0(X1_a, " - ", X1_b, " vs. ", X2_b), paste0(X1_b, " - ", X1_a, " vs. ", X2_a))) %>%
  select(c(X1,X2,title))
names <- pairs$title


library(webshot)
setwd("graphs/")
pairs <- apply(pairs, 1, function(r) {
  message(paste0("\n#############################################################\n","\t\t\t",r[1]," vs. ", r[2],
                 "\n#############################################################\n"))
  df1 <- results[[r[1]]]$edges 
  df2 <- results[[r[2]]]$edges
  
  obj1 <- CreateDifferentialGraph(df1, df2)
  obj2 <- CreateDifferentialGraph(df2, df1)
  
  scale.min <- Reduce(min, lapply(list(obj1, obj2), function(r) min(r$graph)))
  scale.max <- Reduce(max, lapply(list(obj1, obj2), function(r) max(r$graph)))
  breaksList <- seq(scale.min, scale.max, length.out=10)
  cols <- (colorRampPalette(brewer.pal(length(breaksList), "Reds"))(length(breaksList)))
  
  ggsave(width=9, height=4, filename = paste0(r[3], ".heatmap.pdf"), plot=plot_grid(
    pheatmap(obj1$graph, display_numbers = obj1$disp, cluster_rows = F, cluster_cols = F, main = paste(r[1], "Except", r[2]),
             cellwidth = 30, cellheight = 30, legend = F, breaks = breaksList, color = cols, fontsize_number = 14)[[4]], 
    pheatmap(obj2$graph, display_numbers = obj2$disp, cluster_rows = F, cluster_cols = F, main = paste(r[2], "Except", r[1]),
             cellwidth = 30, cellheight = 30, legend = T, breaks = breaksList, color = cols, fontsize_number = 14)[[4]]))
  
  ThreeWay.SankeyPlot(obj1, obj2, title = r[3]) %>% saveNetwork(paste0(r[3], ".sankey.html"))
  webshot::webshot(paste0(r[3], ".sankey.html"), paste0(r[3], ".sankey.png"))
  
  res <- list(obj1, obj2)
  names(res) <- c(r[1], r[2])
  return(res)
})
setwd("../")
names(pairs) <- names
rm(names)
##############################################################################################






##############################################################################################
##                                    DE Genes With DE Edges                                ##
##############################################################################################
# Keep only DE edges where lig-rec are both DE (for Epi and Mes)
conv <- read.csv("Gene_Conversions.txt", sep = "\t")
epi_de <- read.csv("Epi_DESeq2_Results_Embryo_vs_Adulthood_sig_genes_w_gene_annotation.txt", sep = "\t")
mes_de <- read.csv("Mes_DESeq2_Results_Embryo_vs_Adulthood_sig_genes.txt", sep = "\t")
epi_de$X <- plyr::mapvalues(epi_de$X, conv$symbolMM, conv$symbolHS)
mes_de$X <- plyr::mapvalues(mes_de$X, conv$symbolMM, conv$symbolHS)

df <- unique(do.call(rbind, lapply(pairs, function(p) rbind(p[[1]]$edges, p[[2]]$edges))))
df <- df[df$lcell %in% c("Epi","Mes") | df$rcell %in% c("Epi","Mes"),]

df$orig.lig <- df$lig
df$orig.rec <- df$rec

for (term in c(" complex", "_complex1", "_complex2", "_BMPR2"," receptor")) {
  df$lig <- gsub(term, "", df$lig)
  df$rec <- gsub(term, "", df$rec)
}
df <- df[((df$lcell == "Epi" | df$rcell == "Epi") & (df$lig %in% epi_de$X | df$rec %in% epi_de$X)) |
          (df$lcell == "Mes" | df$rcell == "Mes") & (df$lig %in% mes_de$X | df$rec %in% mes_de$X),]

de.genes <- do.call(rbind, 
                    list(unique(df[,c(3,4,7,8)]), 
                         c("IL1B", "IL1R1", "IL1B", "IL1R1"),
                         c("VEGFA","FLT1", "VEGFA","FLT1"),
                         c("CSF1","CSF1R", "CSF1","CSF1R")))

de.genes$ligMM <- plyr::mapvalues(de.genes$lig, conv$symbolHS, conv$symbolMM)
de.genes$recMM <- plyr::mapvalues(de.genes$rec, conv$symbolHS, conv$symbolMM)


# Summarized expression of selected DE genes for DotPlot
library(Seurat)
obj <- readRDS("Seurat_All_Cells_wo_Doublets_filtered_082020.rds")
cell.info <- read.delim('Cell_info.txt', row.names = 1)
cell.info <- subset(cell.info,!(cell_type_subclusters %in% c('group','Immune','Mes','Plasma')))

genes <- unique(c(de.genes$ligMM, de.genes$recMM))
obj@meta.data$cell_type <- cell.info[rownames(obj@meta.data),"cell_type"]
expression <- FetchData(obj, c(genes, "Ventricle","Age","cell_type"))
genes <- genes[genes %in% colnames(expression)]

alpha.summary <- expression %>% 
  dplyr::group_by(Ventricle,Age,cell_type) %>% 
  dplyr::summarise(across(genes, ~ 100*mean(.x > 0))) %>% 
  melt(id.var=c('Ventricle','Age','cell_type'), value.name = "alpha", variable.name="gene") 

mu.summary <- expression %>% 
  dplyr::group_by(Ventricle,Age,cell_type) %>% 
  dplyr::summarise(across(genes, ~ median(.x[.x > 0]))) %>%  # mean(.x)
  melt(id.var=c('Ventricle','Age','cell_type'), value.name = "mu", variable.name="gene") 
mu.summary$mu[mu.summary$mu > 2.5] <- 2.5
mu.summary$mu[is.na(mu.summary$mu)] <- 0

gene.summary <- dplyr::inner_join(alpha.summary, mu.summary,by=c('Ventricle','Age','cell_type','gene'))
gene.summary$cell_type <- factor(gene.summary$cell_type,levels=c('Immune','NeuroGlia','Endo', 'Mes', 'Epi'))
gene.summary$Age <- factor(gene.summary$Age,levels=c('Embryo','Adult','Aged'))
gene.summary$Ventricle <- factor(gene.summary$Ventricle,levels=c('TCP','DCP','HCP'))
gene.summary$gene <- as.character(gene.summary$gene)



material.heat <- function(n) {
  mh = c("#283593", #indigo 800
         "#3F51B5", #indigo
         "#2196F3", # blue
         "#00BCD4", #cyan
         # "#4CAF50", #green
         "#8BC34A", #light green
         "#CDDC39", #lime
         #"#FFEB3B", #yellow
         "#FFC107", #amber
         "#FF9800", #orange
         "#FF5722"#, #deep orange)
         #"#f44336"
         )
  colorRampPalette(mh)(n)
}

pdf(height = 5, width = 5, file = "graphs/examples.exp.all.dotplot.pdf")
for (i in 1:nrow(de.genes)) {
  message(paste(de.genes[i,5:6], collapse='_'))
  s <- subset(gene.summary, gene %in% unname(unlist(de.genes[i,5:6])) )
  message(paste("Got expression for:", paste(unique(s$gene), collapse = ",")))

  if (length(unique(s$gene)) == 2) {
    message("Plotting")
    print(ggplot(s,aes(gene,group=Ventricle,cell_type,color=mu,size=alpha)) +
      geom_point(position=position_dodge(width=0.8)) +
      facet_wrap("Age", ncol = 1) +
      theme_classic() +
      theme(axis.text.x = element_blank()) +
      labs(x="", y="") +
      scale_colour_gradientn(colours=material.heat(50)) + labs(title=paste(de.genes[i,5:6], collapse=' - ')))
  }
  message("")
}
while (!is.null(dev.list())) dev.off()



#Selected Ligand-Receptor pairs for DotPlots
pdf(height = 5, width = 14, file = "graphs/selected.examples.exp.all.sup.pdf")

# Age Differential
selected.genes <- c("Cxcl12", "Ackr3", "Efna5", "Ephb2", "Igf2", "Igf2r")
s <- subset(gene.summary, gene %in% selected.genes)
s$gene <- factor(s$gene, levels = selected.genes)
print(ggplot(s,aes(gene,group=Ventricle,cell_type,color=mu,size=alpha)) +
        geom_point(position=position_dodge(width=0.5)) +
        facet_wrap("Age", ncol = 1) +
        theme_classic() +
        labs(x="", y="") +
        scale_colour_gradientn(colours=material.heat(50)))

# Ventricle Differential
selected.genes <- c("Nov", "Notch1", "Bmp5", "Bmpr1b", "Dlk1", "Notch2", "Notch3", "Notch4","Rspo2","Lgr4")
s <- subset(gene.summary, gene %in% selected.genes)
s$gene <- factor(s$gene, levels = selected.genes)
print(ggplot(s,aes(gene,group=Ventricle,cell_type,color=mu,size=alpha)) +
        geom_point(position=position_dodge(width=0.5)) +
        facet_wrap("Age", ncol = 1) +
        theme_classic() +
        labs(x="", y="") +
        scale_colour_gradientn(colours=material.heat(50)))


selected.genes <- c("Igf2", "Igf1r", "Pdgfa", "Pdgfra", "Pdgfb", "Pdgfrb", "Wnt5a", "Ror1", "Vegfa","Flt1","Efna5","Epha4","Epha5","Epha7","Csf1","Csf1r")
s <- subset(gene.summary, gene %in% selected.genes)
s$gene <- factor(s$gene, levels = selected.genes)
print(ggplot(s,aes(gene,group=Ventricle,cell_type,color=mu,size=alpha)) + 
        geom_point(position=position_dodge(width=0.8)) +
        facet_wrap("Age", ncol = 1) +
        theme_classic() +
        labs(x="", y="") +
        scale_colour_gradientn(colours=material.heat(50)))

while (!is.null(dev.list())) dev.off()



# Over ligand-receptors of interest - Plot edges indicator matrix
all_edges <- do.call(rbind, lapply(seq_along(results), function(i)
  data.frame(case=names(results)[i], results[[i]]$edges) %>% 
    separate(case, c("Age", "Region"), "_")))
#all_edges <- all_edges[all_edges$Age != "SCEmbryo",]
all_edges$Age <- factor(all_edges$Age, levels=c("Embryo","Adult","Aged"))
colnames(all_edges)[2] <- "Ventricle"
all_edges$Ventricle <- factor(all_edges$Ventricle,levels=c('TCP','DCP','HCP'))

write_excel_csv(all_edges[,1:6], "ligand.receptor.edges.csv")


all_tuples <- unique(all_edges[,c("Age","Region","lcell","rcell")])
# Outer-join of all interactions in question with all combinations of interacting cells (by region+age)
df <- base::merge(data.frame(ex.lig=de.genes$orig.lig, ex.rec=de.genes$orig.rec), all_tuples, all=T)
df$interaction <- paste(df$ex.lig, df$ex.rec, sep = " - ")


df <- base::merge(df, all_edges, all=T)
df <- df[mapply(grepl, pattern=df$ex.lig, x=df$lig),]
df <- df[mapply(grepl, pattern=df$ex.rec, x=df$rec),]
df <- base::merge(all_tuples, df, all.x=T)

df[,c("lweight","rweight")][is.na(df[,c("lweight","rweight")])] <- 0
df <- dcast(df, Age+Region+lcell+rcell~interaction, value.var = "lweight", fun.aggregate = sum)
df <- df[,-ncol(df)]
df[,-c(1:4)] <- df[,-c(1:4)] > 0


cols <- brewer.pal(3, "Dark2")
names(cols) <- levels(df$Age)
annotation.colors <- list(Age = cols)

ggsave(width=20, height=10, dpi=400, filename = "graphs/preserved.interactions.pdf", plot=plot_grid(plotlist = 
  lapply(levels(df$Region), function(r) {
    d <- df[df$Region == r,]
    ident <- paste(paste(d$lcell, d$rcell, sep = " - "))
    annotations <- data.frame(Age=d$Age, row.names = rownames(d))
    d <- apply(d[,-c(1:4)], 2, as.numeric)
    rownames(d) <- rownames(annotations)
    pheatmap(d,
             cluster_rows = F, cluster_cols = F, color = c("lightgrey","darkred"), border_color = NA, fontsize_row = 7, cellwidth = 20,
             annotation_row = annotations,  labels_row = ident, angle_col=45, annotation_colors = annotation.colors,
             legend = F, annotation_names_row = F, annotation_legend = r == "TCP", treeheight_row = 0, main = r)[[4]]
  }), nrow=1, rel_widths = c(1,1,1.2)))



s <- df %>% unite(setting, Age, Region, sep=" - ") %>% 
  select(-one_of("lcell","rcell")) %>% 
  group_by(setting) %>% 
  summarise_all(sum) %>% 
  column_to_rownames("setting") %>% as.data.frame()

ggsave(width=10, height=4, dpi=300, filename="graphs/ligand.receptor.settings.count.pdf", plot=
         pheatmap(s, display_numbers = s, cluster_rows = F, cluster_cols = F, color = RColorBrewer::brewer.pal(max(s), "Reds"))[[4]])

pdf("graphs/examples.pdf")
for (c in colnames(df)[-(1:4)]) {
  d <- data.frame(df[1:4], lig.rec=as.numeric(df[,c])) %>%
    dcast(Age+lcell+rcell~Region, value.var = "lig.rec") %>% 
    unite(cells, lcell, rcell, sep=" - ")
  row.annotations <- data.frame(Age=as.character(d[,"Age"]), row.names=rownames(d))
  rownames(d) <- rownames(row.annotations)
  
  pheatmap(d[,3:5], cluster_rows = F, cluster_cols = F, color = c("white","red3"), legend = F, main=c, width = 7,
           labels_row = d[,"cells"], annotation_row = row.annotations, annotation_colors = annotation.colors)
  
}
while (!is.null(dev.list())) dev.off()


rm(df, cols, annotation.colors, all_edges, all_tuples)
##############################################################################################

