args<-commandArgs(TRUE)
#rm(list=ls())
library(ggplot2)

sample_name <- args[0]
dir <- args[1]
coverage_file <- args[2]

#setwd("/Users/aminzia/Desktop/coverage")
setwd(dir)
#cover = read.table("acmg_coverage_genes.txt", header=T, row.names=NULL, as.is=T)
cover = read.table(coverage_file, header=T, row.names=NULL, as.is=T)


myClip <- function(x, a, b) {
  ifelse(x <= a,  a, ifelse(x >= b, b, x))
}

genes <- cover$gene
dindex <- which(duplicated(genes))
genes <- genes[-dindex]

nrow_max = 0
for (gene in genes) {
  coverg <- cover[which(cover$gene == gene), ]
  rownames(coverg) <- NULL

  if (nrow(coverg) > nrow_max) {
    nrow_max = nrow(coverg)
  }   
}

sw = 0
for (gene in genes) {
  coverg <- cover[which(cover$gene == gene), ]
  rownames(coverg) <- NULL
  chrom <- coverg[1,1]
  coverage <- coverg$Depth_for_case0017
  coverage.mean = as.numeric(mean(coverage))
  coverage.std = as.numeric(3*sqrt(var(coverage)))
  
  genome <- data.frame(as.data.frame(matrix(41.97, nrow = length(coverage), ncol = 1)))
  
  new_row <- c(gene, as.numeric(coverage.mean), as.numeric(coverage.std))
  coverage.matrix.col <- data.frame(gene=rep(gene, length(coverage)), coverage=coverage)

  if (sw==1){
    cover_stats <- rbind(cover_stats, new_row)
    coverage.matrix <- rbind(coverage.matrix, coverage.matrix.col)
  } else {
    cover_stats <- data.frame(gene=gene, mean_cover=as.numeric(coverage.mean), std_cover=as.numeric(coverage.std), stringsAsFactors=FALSE)    
    coverage.matrix <- coverage.matrix.col
    sw = 1
  }
}  
  gr<-ggplot() +
    #geom_bar(stat="identity", fill="white", colour="black", position = position_dodge()) +
    geom_point(data=cover_stats, aes(x=gene, y=as.numeric(mean_cover)), size=4, shape=21, fill="white") +
    geom_errorbar(data=cover_stats, aes(x=gene, ymin=myClip(x=as.numeric(mean_cover)-as.numeric(std_cover),a=0, b=200) , ymax=as.numeric(mean_cover)+as.numeric(std_cover)), width=0.5, size=1, color="blue") +
    #geom_errorbar(data=cover_stats, aes(x=gene, ymin=as.numeric(mean_cover), ymax=as.numeric(mean_cover)), width=0.2, size=1, color="blue") +
    theme(axis.text.x=element_text(angle=-90, colour="black", size=9)) +
    #theme(text = element_text(size=10)) +
    ylab("Coverage")+
    xlab("") +
    #theme(plot.margin = unit(c(15, 5, 5, 5), "cm"))
    #scale_y_continuous(limits=c(0, 100), oob=rescale_none) +
    #coord_cartesian(ylim=c(0,100)) +
    ylim(0, 200)+
    #labs(title = "NA12887-CLIA")
    labs(title = sample_name)
  
  ggsave(filename = "coverage_bar_3sigma.png", plot = gr, scale = 1, width = par("din")[1], height = par("din")[2], units = c("in","cm", "mm"), dpi = 300)

