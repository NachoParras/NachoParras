#! /usr/bin/env Rscript
# x,y graph

#README. It takes a tabulated txt file with the names of the datasets to overlap, next the folder with the symbolic names of them, and the output folder. It will compute the overlap between those genelist and create a venn diagram of the result

library(ggplot2)
library(optparse)
library(RColorBrewer)
library(ggvenn)
library(VennDiagram)

################################################################
# OPTPARSE
################################################################
option_list <- list(
  make_option(c("-i", "--input"), type="character",
    help="Tabulated file with the name of each dataset to analyze"),
  make_option(c("-d", "--input_dir"), type="character",
    help="Directory with the symbolic names of the datasets to analyze"),
  make_option(c("-c", "--column"), type="character", default="genes_tag", 
    help="Name of the column with the parameter to select genes"),
  make_option(c("-p", "--param"), type="character", default="NA",
    help="Name of the feature to filter selected genes"),
  make_option(c("-o", "--output"), type="character", default="results", 
    help="Output figure file"),
  make_option(c("-f", "--output_format"), type="character", default="pdf",
    help="pdf or jpeg file output format")

)
opt <- parse_args(OptionParser(option_list=option_list))

################################################################
## MAIN
################################################################

#1. Indicate the output format
if (opt$output_format == "pdf"){
  pdf(paste(opt$output, '.pdf', sep=""))
}else if(opt$output_format == "jpeg"){
  jpeg(paste(opt$output, '.jpeg', sep=""))
}
#2. Read every dataset and make a list witht the Prevalent DEGs from them 
input <- read.table(opt$input, sep="\t", header = FALSE)
input_names <- input[,1]
n <- length(input_names)
	if(n<1){
	return("It seems there is no input files to read")
	}else{
	generate_input <- function(s){
	i <- 1
	l <- vector(mode = "list", length = n)
	names(l) <- input_names
		for (i in i:n){
		data <- read.table(paste0(opt$input_dir, input_names[i], collapse=""), sep="\t", header = TRUE, row.names=1)
		#df <- data[data$col==param,]
		col <- which(colnames(data)==opt$column)
		if(opt$param=="NA"){
		df <- data
		}else{
		df <- data[data[,col]==opt$param,]
		}
		l[[i]] <- rownames(df)
		#print(l[[i]])
		}
	return(l)
	}}

DEGs_list <- generate_input(n)

#3. Calculate and save te overlap
overlap <- calculate.overlap(x = DEGs_list)
name_overlap <-  paste(opt$output,".txt", sep ="")
lapply(overlap, function(x) write.table( data.frame(x), name_overlap, append = T, sep = "\t"))

#4. Create VennDiagram
create_venny <- function(ls){
	ggvenn(ls, 
      	#fill_color = colors,
      	stroke_size = 0.5, 
      	set_name_size = 4,
      	columns = names(ls),
      	show_percentage = TRUE,
   	 )
}
create_venny(DEGs_list)

dev.off()
       
