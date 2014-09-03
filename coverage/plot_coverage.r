#rm(list=ls())
library(ggplot2)
#setwd("/Users/aminzia/Desktop/coverage")

sample_name = "case0017"
genome_mean_coverage = 41.97
input_dir = "~/Clinical_Service/cases/case0017/medgap-2.0/QC-0.1/new3/"
#input_file = paste("coverage_hist_Q0_withgenes_". )
input_file = "coverage_hist_Q0_withgenes_dcm.txt"
#input_file = "coverage_hist_Q0_withgenes_acmg.txt"
#input_file = "acmg_genes_coverage_hist.txt"
#input_file = "acmg_coverage_hist_genes.txt"

input_path = paste(input_dir, input_file, sep="")
cover = read.table(input_path, header=T, row.names=NULL, as.is=T)

myClip <- function(x, a, b) {
  ifelse(x <= a,  a, ifelse(x >= b, b, x))
}

#gene1 = subset(cover, cover$gene=="SDHB")
#gene2 = subset(cover, cover$gene=="MUTYH")

genes <- cover$gene
dindex <- which(duplicated(genes))
genes <- genes[-dindex]

#genes <- c("SDHB", "MUTYH")
#gene <- "SDHB"
nrow_max = 0
for (gene in genes) {
  coverg <- cover[which(cover$gene == gene), ]
  rownames(coverg) <- NULL

  if (nrow(coverg) > nrow_max) {
    nrow_max = nrow(coverg)
  }   
}

sw = 0
gene="TTN"
for (gene in genes) {
  coverg <- cover[which(cover$gene == gene), ]
  rownames(coverg) <- NULL
  chrom <- coverg[1,1]
  coverage <- coverg$Depth
  coverage.median = as.numeric(median(coverage))
  coverage.mean = as.numeric(mean(coverage))
  coverage.std = as.numeric(3*sqrt(var(coverage)))
    
  new_row <- c(gene, as.numeric(coverage.median), as.numeric(coverage.mean), as.numeric(coverage.std))
  #coverage <- c(coverage, rep(NA, nrow_max - length(coverage)))
  coverage.matrix.col <- data.frame(gene=rep(gene, length(coverage)), coverage=coverage)

  if (sw==1){
    cover_stats <- rbind(cover_stats, new_row)
    coverage.matrix <- rbind(coverage.matrix, coverage.matrix.col)
  } else {
    cover_stats <- data.frame(gene=gene, median_cover=as.numeric(coverage.median), mean_cover=as.numeric(coverage.mean), std_cover=as.numeric(coverage.std), stringsAsFactors=FALSE)    
    coverage.matrix <- coverage.matrix.col
    sw = 1
  }
}  
  #gf <- coverage.matrix
  #ggplot(gf, aes(x=gene, y=coverage, fill=gene)) + geom_boxplot() +  guides(fill=FALSE) + coord_flip() +
  #ylim(0, 100) +
  #geom_boxplot(stat = "identity")
  #geom_bar("boxplot", "jitter")
  #fill=gear
  #main="Mileage by Gear Number" +
  #xlab="" +
  #ylab="Miles per Gallon"
  #gr
  #gf <- cover_stats
  #ggplot(gf, aes(x=gene, y=mean_cover, fill=gene)) + geom_bar() +  guides(fill=FALSE) + coord_flip()
  #ggplot(gf, aes(x=gene, y=mean_cover, fill=gene)) + geom_bar() +  guides(fill=FALSE)

  #gr<-ggplot(cover_stats, aes(x = gene, y = mean_cover)) +
  gr<-ggplot() +
    #geom_bar(stat="identity", fill="white", colour="black", position = position_dodge()) +
    geom_point(data=cover_stats, aes(x=gene, y=as.numeric(median_cover)), size=2, shape=21, fill="white") +
    geom_errorbar(data=cover_stats, aes(x=gene, ymin=myClip(x=as.numeric(mean_cover)-as.numeric(std_cover),a=0, b=200) , ymax=as.numeric(mean_cover)+as.numeric(std_cover)), width=0.5, size=1, color="blue") +
    #geom_errorbar(data=cover_stats, aes(x=gene, ymin=as.numeric(mean_cover), ymax=as.numeric(mean_cover)), width=0.2, size=1, color="blue") +
    theme(axis.text.x=element_text(angle=-90, colour="black", size=9)) +
    #theme(text = element_text(size=10)) +
    ylab("Coverage (median, 3-sigma)")+
    xlab("") +
    #theme(plot.margin = unit(c(15, 5, 5, 5), "cm"))
    #scale_y_continuous(limits=c(0, 100), oob=rescale_none) +
    #coord_cartesian(ylim=c(0,100)) +
    ylim(0, 100)+
    labs(title = sample_name)
  #gr
  #ggsave(filename = paste("graphs/coverage_bar_3sigma_",sample_name,".png", sep=""), plot = gr, scale = 1, width = par("din")[1], height = par("din")[2], units = c("in","cm", "mm"), dpi = 200)
  ggsave(filename = paste("graphs/coverage_bar_3sigma_",sample_name,".png", sep=""), plot = gr, scale = 1, width = 10, height = 4, units = "in", dpi = 200)
#file <- "graphs/coverage_bar_3sigma.png"
  #png(filename=file)
  #plot(gr)
  #dev.off()

#}
#file <- "graphs/acmg_genes_ave_coverage.png"
  #png(filename=file)
  #plot(gr)
  #dev.off()
gene="TTN"
#for (gene in genes) {
  coverg <- cover[which(cover$gene == gene), ]
  rownames(coverg) <- NULL
  chrom <- coverg[1,1]
  coverage <- coverg$Depth
  genome <- data.frame(as.data.frame(matrix(genome_mean_coverage, nrow = length(coverage), ncol = 1)))
  locas <- coverg$Locus
  Exon <- factor(coverg$transcript)
  gr2<-ggplot(fill=gene) +
  #geom_point(data=coverg, aes(y=coverage, x=factor(locas), color=Exon), size=1, shape=2, fill="white") +
  geom_point(data=coverg, aes(y=coverage, x=factor(locas), color=Exon), size=1, fill="white") +
  #geom_line() +
  #geom_point(data=coverg, aes(y=coverage, x=as.factor(locas))) #+
  #geom_line(data=coverg, aes(y=coverage, x=locas, colour=gene)) +
  #geom_line(data=genome, aes(y=genome$V1, x=as.factor(locas), colour="GENOME")) +
  #scale_colour_manual(values=c("red","blue")) +
  #theme(legend.title=element_blank()) +
  #theme(legend.position = "none") +
  #scale_x_discrete(breaks=factor(locas)) +
  #scale_x_continuous(breaks=seq(1:3303)) +
  labs(title = gene) +
  xlab(chrom) +
  #theme(axis.text.x=element_text(angle=-90, colour="black", size=1)) +
  theme(axis.text.x=element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  ylab("Coverage") +
  coord_fixed(ratio=20) +
  ylim(0, 100)
  #gr2
  file <- paste("graphs/", gene, ".png", sep="")
  #png(filename=file, width = 10, height = 4, units = 'in', res=200)
  ggsave(filename = file, plot = gr2, scale = 1, width = par("din")[1], height = par("din")[2], units = c("in","cm", "mm"), dpi = 200)
  #png(filename=file, res=200)
  #plot(gr2)
  #dev.off()
#}

#coverage1 = cover$Depth_for_NA12878
#proportion_coverge = gene1$Average_Depth_sample
#locas <- cover$locas

#
##for(i in 2:length(sample_list)) {
  #print(sample_list[[i]][[2,1]])
  #data <- combineDataSets(data, sample_list[[i]], sample_list[[i]][[2,1]])
#  data <- combineDataSets(data, sample_list[[i]], sample_list[[i]][[2,1]])
#}
#  
## Remove the "gte_" part from the coverage function
#data$coverage<-sapply(data$coverage, function(x) sub("gte_", "",x))
# 
## Convert to reasonable data types...
#data$coverage<-as.numeric(data$coverage)
#data$proportion_coverge<-as.numeric(as.character(data$proportion_coverge))
# 
#ggplot(data=data, aes(x=coverage/max(coverage), y=proportion_coverge, color=sample)) +
#  geom_line() +
#  labs(title = "Cumulative coverage compassion") +
#  scale_y_continuous(breaks=seq(0,1,0.1))#
#
#createInitialDataset <- function(x, sampleName) {
#  extraLineRemoved <- t(removeExtraLine(x))
#  data<-data.frame(extraLineRemoved, rep(sampleName))
#  colnames(data)<-c("coverage","proportion_coverge","sample")
#  data
#}

#combineDataSets <- function(data, x, sampleName) {
#  extraLineRemoved <- t(removeExtraLine(x))
#  x<-data.frame(extraLineRemoved, rep(sampleName))
#  colnames(x)<-c("coverage","proportion_coverge","sample")
#  data<-rbind(data, x)
#  data
#}

#removeExtraLine <- function(x) {
#  x[,2:length(x)]
#}#
