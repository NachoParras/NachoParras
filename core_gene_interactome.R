#! /usr/bin/env Rscript

library(RColorBrewer)
library(ggplot2)
library(optparse)
library(STRINGdb)
library(igraph)
library(biomaRt)
library(qgraph)

################################################################
# OPTPARSE
################################################################
option_list <- list(
  make_option(c("-d", "--data_file"), type="character",
              help="Tabulated file as a geneset of interest"),
  make_option(c("-o", "--output"), type="character", default="results",
              help="Output figure file"),
  make_option(c("-c", "--column"), type="character", default="Symbol",
              help="Name of the column with the gene identifiers"),
  make_option(c("-n", "--neighbors"), type="logical", default=FALSE,
	     help="if TRUE, neighbors will be added to the network"),
  make_option(c("-p", "--param"), type="character", default="NA",
             help="Feature to filter genes"),
  make_option(c("-f", "--output_format"), type="character", default="pdf",
              help="pdf or jpeg file output format"),
  make_option(c("-t", "--threshold"), type="numeric", default=400,
              help="threshold to set in the ppi interaction network"),
  make_option(c("-g", "--communities"), type="logical", default=FALSE,
              help="If TRUE, communities will be drown in the graph"),
  make_option(c("-s", "--specie"), type="numeric", default=10090,
              help="specie to download the target interaction network"),
  make_option(c("-u", "--organism"), type="character", default="mmusculus",
              help="organism to match network identifiers")

)
opt <- parse_args(OptionParser(option_list=option_list))

################################################################
## MAIN
################################################################

#1.Load differential gene expression data or a gene list of interest
if(opt$param=="NA"){
data <- read.table(file = opt$data, sep = "\t", header = TRUE)
}else{
data <- read.table(file = opt$data, sep = "\t", header = TRUE)
a <- which(colnames(data)==opt$column)
data <- data[data[a,]==opt$param,]
}

#2. Load Mouse Interactome
string_db <- STRINGdb$new( version="11", species=opt$specie, score_threshold=opt$threshold, input_directory="")

#3. Map loaded data in STRING interaction network
mapped_data <- string_db$map( data, opt$column, removeUnmappedRows = TRUE ) 
hits <- mapped_data$STRING_id
if(opt$neighbors==TRUE){
genes_neighbors <- string_db$get_neighbors( c(hits) )
string_ids <- c(hits, genes_neighbors)
g <- string_db$get_subnetwork( string_ids )
}else{
g <- string_db$get_subnetwork( hits )
}

#4. Transform STRING names (ensembl_peptide_id) into mgi gene symbol
eliminate <- paste(as.character(opt$specie), ".", sep="")
V(g)$name <- sub(eliminate,"",V(g)$name)
mart<- useMart("ensembl")
mart_organism <- paste0(opt$organism, "_gene_ensembl", collapse="")
mart <- useDataset(mart_organism, mart)
mart_results <- getBM(attributes = c("ensembl_peptide_id", "mgi_symbol", "description"), filters = "ensembl_peptide_id", values = V(g)$name, mart = mart)
ix <- match(V(g)$name, mart_results$ensembl_peptide_id)
ix <- ix[!is.na(ix)]
newnames <- V(g)$name
newnames[match(mart_results[ix,'ensembl_peptide_id'], newnames)] <- mart_results[ix, 'mgi_symbol']
V(g)$name  <- newnames
write.table(mart_results, file = paste0(opt$output, "interacting_genes_description.txt", collapse=""), dec=".", sep="\t", quote = FALSE)

#5. Determine output format
if (opt$output_format == "pdf"){
  pdf(paste(opt$output, '.pdf', sep=""))
}else if(opt$output_format == "jpeg"){
  jpeg(paste(opt$output, '.jpeg', sep=""))
}

#6. Plot the network
#6.1 Eliminate nodes with degree=0
deg <- degree(g, mode="all")
sb <- deg[deg>3]
sg <- subgraph(g, names(sb))

#6.2. Eliminate possible duplicated adges
sg <- simplify(sg, remove.multiple=TRUE)

#6.3. Graph appearance
ecol <- rep("gray80", ecount(sg))
ew <- rep (0.3, ecount(sg))
el <- get.edgelist(sg)
v <- vcount(sg)
e <- get.edgelist(sg, names=FALSE)
l <- qgraph.layout.fruchtermanreingold(e, vcount=vcount(sg), area=50*(vcount(sg)^2))
#l <- layout_nicely(sg)
if(opt$communities==FALSE){
plot(sg, layout=l,
     vertex.label.color="gray10", vertex.label.cex=0.2, vertex.label.font=1,
     vertex.size=5, vertex.shape="circle", vertex.frame.color="gray1", 
     edge.color=ecol, edge.width=ew, edge.arrow.mode=0, asp=0)
}else{
gc <- edge.betweenness.community(sg)
gc_names <- V(sg)$name
gc_mem <- gc$membership
#V(gc)$name <- gc_names
df <- data.frame(vertex_name=gc_names, membership=gc_mem) 
write.table(df, file = paste0(opt$output, "genes_communities.txt", collapse=""), dec=".", sep="\t", quote = FALSE)

V(sg)$community <- gc$membership
colors <- brewer.pal(n = 8, name = "RdBu")
colors <- colorRampPalette(colors)(5)
plot(gc, sg, layout=l,
	vertex.label.color="gray10", vertex.label.cex=0.2, vertex.label.font=1, 
	vertex.size=5, vertex.shape="circle", verex.frame.color="gray1", vertex.color=colors[V(sg)$community], 
	edge.color=ecol, edge.width=ew, edge.arrow.mode=0, asp=0)
}
dev.off()
