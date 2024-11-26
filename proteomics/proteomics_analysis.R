library(Cairo)
library(dplyr)
library(ggplot2)

set.seed(417)

proteingroups.fi <- "proteinGroups.txt" #provide path to mass spec output file
setwd("") #set to directory in which you want output table and figures saved

##Protein enrichment analysis
#Format gene names
proteinGroups <- read.delim(proteingroups.fi, stringsAsFactors=FALSE)
proteinGroups <- proteinGroups[proteinGroups$Reverse != "+",]
proteinGroups <- proteinGroups[proteinGroups$Only.identified.by.site != "+",]
proteinGroups <- proteinGroups[proteinGroups$Potential.contaminant != "+",]
proteinGroups <- proteinGroups[which(! is.na(proteinGroups$id)),] 

proteinGroups$Gene.names <- sapply(strsplit(sapply(strsplit(proteinGroups$Fasta.headers, "GN="), tail, 1)," "),head,1)
allproteinIDs <- strsplit(proteinGroups$Majority.protein.IDs,"\\;")
allproteinIDs <- lapply(allproteinIDs, tail, 1)
allproteinIDs <- as.character(allproteinIDs)
proteinGroups$Gene.names <- ifelse(proteinGroups$Gene.names == "", proteinGroups$Protein.IDs, proteinGroups$Gene.names)#fill empty gene name entries with UniprotIDs

#Reduce, transform, and impute
proteinGroupsReduced <- subset(proteinGroups, select=c(Gene.names,Majority.protein.IDs,grep("LFQ.intensity.",names(proteinGroups))))

proteinGroupsRed <- subset(proteinGroups, select=grep("LFQ.intensity.",names(proteinGroups)))       
proteinGroupsRed[,grepl("LFQ.intensity.",names(proteinGroupsRed))] <- log2(proteinGroupsRed[,grepl("LFQ.intensity.",names(proteinGroupsRed))]) #log2 transformation
proteinGroupsRed <- do.call(data.frame,lapply(proteinGroupsRed, function(x) replace(x, is.infinite(x),NA)))

#filter rows based on valid values
proteinGroupsRedlogical <- !is.na(proteinGroupsRed) #make logical dataframe for valid (non NA) entries

results <- apply(proteinGroupsRedlogical,1,function(x){
  groups <- unique(gsub("_[0-9]$","",names(x)))
  tempsum <- numeric()
  for (group in 1:length(groups)){
    tempsum[group] <- sum(x[grep(paste0("^",groups[group]),names(x))])
  }
  if (sum(grepl("[2,3]",tempsum)>0)) {y=T} else {y=F}
  return(y)
})

filtered <- proteinGroupsRed[results,]

#imputation
filtered <- apply(filtered[,],2, function(x){
  standd <- sd(x,na.rm = T) #standard deviation of the measured values
  standd1 <- standd*0.25 #downscaling of the sd for imputing values
  meand <- mean(x,na.rm = T) #mean of the measured values
  meand1 <- meand - (1.8*standd) #downscaling of the mean for imputing values
  sapply(x, function(y){
    if(is.na(y)){
      y <- rnorm(1, mean= meand1, sd=standd1)
      return(y)
    } else  if(!is.na(y)){return(y)}
  })
})
filtered <- data.frame(filtered)

gene.name <- proteinGroupsReduced$Gene.names[results]
filtered1 <- cbind(gene.name,filtered)
filtered_final <- cbind(proteinGroups$Protein.IDs[results],filtered1)
names(filtered_final)[1] <- "Protein.IDs"
names(filtered_final)[2] <- "Gene.names"

#get ttest results for DYRK3 vs. GFP and PP1-NIPP1 vs. GFP

ttest <- function(df,g1,g2){
  tet <- sapply(1:nrow(df), function(i){ 
    ttest <- t.test(df[i,grepl(g1,names(df))], filtered[i,grepl(g2,names(df))], var.equal = T) 
    log2fc <- ttest$estimate[2]-ttest$estimate[1]
    pval <- ttest$p.value
    #pval <- (-log10(pval))
    out <- c(log2fc,pval)
    return(out)
  })
  tet <- t(data.frame(tet))
  colnames(tet) <- c("log2FC","pval")
  rownames(tet) <- NULL
  return(tet)
}

dyrk3.df <- ttest(filtered,"GFP_Control","DYRK3_")
pp1.df <- ttest(filtered,"GFP_Control","PP1.NIPP1")
plotdf <- as.data.frame(cbind(dyrk3.df,pp1.df))
colnames(plotdf) <- c("x","p.x","y","p.y")

#save supplementary table
supple <- cbind(filtered_final[,c("Protein.IDs","Gene.names")],plotdf)
names(supple) <- c("Protein.IDs","Gene.names","Log2FC.DYRK3_vs_GFP","p.DYRK3_vs_GFP",
                   "Log2FC.PP1-NIPP1_vs_GFP","p.PP1-NIPP1_vs_GFP")
write.table(supple,"SuppTable1_Proteomics.csv",row.names = FALSE,sep=',',quote=FALSE)

#plot DYRK3 vs. PP1-NIPP1 enrichment
plotdf$sig <- plotdf$p.x < 0.05 & plotdf$p.y < 0.05 & plotdf$x >= 2 & plotdf$y >= 2
plotdf$sig[which(plotdf$sig == "TRUE")] <- "enriched"
plotdf$sig[which(plotdf$sig == "FALSE")] <- "non-sig" #p > 0.05 or |l2fc| < 2"

#plot l2FC for (DYRK3 vs. GFP) vs. (PP1-NIPP1 vs. GFP), 
#with speckle genes labelled and coloring by significance (l2FC >= 2, p < 0.05) in both vs. not both 
g <-  ggplot(plotdf, ggplot2::aes(x, y),echo=T) +
  geom_point(aes(colour=sig),alpha=0.6) + scale_color_manual(values = c("#E69F00","#999999")) + theme_bw() + theme_classic(base_size = 16) + 
  xlab("Log2FC(DYRK3/GFP)") + ylab("Log2FC(PP1-NIPP1/GFP)") +
  geom_vline(xintercept=2, linetype="dashed", color = "#E69F00") + geom_hline(yintercept=2, linetype="dashed", color = "#E69F00") +
  ggrepel::geom_text_repel(data=plotdf[grepl("SRS|SRRM|SON",filtered1$gene.name),],
                           min.segment.length = unit(0, 'lines'),
                  aes(label = filtered1$gene.name[grepl("SRS|SRRM|SON",filtered1$gene.name)]),max.overlaps=Inf)
print(g)  

# Define export parameters
ty <- "pdf" # File type
wi <- 400 # Width in pixels
he <- 320 # Height in pixels
re <- 75 # Resolution in dpi
mte <- 10 # Main title text size
# Create and export graph
Cairo(file=paste0("l2fc_scatter.pdf"), type=ty, width=wi, height=he, pointsize=mte, dpi=re)
print(g)
dev.off()


#plot volcano plots
library(Cairo)
library(ggplot2)
group1 <- "DYRK3" #"GFP-PP1m-NIPP1"
group2 <- "DYRK3+GSK626616" 
tet <- sapply(1:nrow(filtered), function(i){ 
  ttest <- t.test(filtered[i,grepl("DYRK3.GSK",names(filtered))], filtered[i,grepl("DYRK3_",names(filtered))], var.equal = T) # .PP1m.NIPP1
  log2fc <- ttest$estimate[2]-ttest$estimate[1]
  pval <- ttest$p.value
  pval <- (-log10(pval))
  out <- c(log2fc,pval)
  return(out)
})
tet <- t(data.frame(tet))

colnames(tet) <- c("log2FC","pval") #turned around, careful, check following code
rownames(tet) <- NULL
plotdf <- data.frame(tet)

positive <- plotdf[,2]>=(-log10(0.05)) & plotdf[,1] >=2

colnames(plotdf) <- c("x","y")


g <-  ggplot2::ggplot(plotdf, ggplot2::aes(x, y),echo=T) +
  ggplot2::geom_point(colour=densCols(plotdf$x,plotdf$y,colramp=colorRampPalette(c("grey","black"))),size=2)+
  ggplot2::theme_bw()+
  ggplot2::xlab(paste0("Log2FC (",group1,"/",group2,")")) +
  ggplot2::ylab("-log[10] p-value") +
  ggplot2::ggtitle(paste("AP",group1,"vs.",group2,sep = " "))+
  ggplot2::theme(plot.title = ggplot2::element_text(size=20, face="bold", vjust=2))+
  ggrepel::geom_text_repel(data=plotdf[grepl("SRS|SRRM|SF3B|SON",filtered1$gene.name),],
                           ggplot2::aes(label = filtered1$gene.name[grepl("SRS|SRRM|SF3B|SON",filtered1$gene.name)]),max.overlaps=Inf) + 
  ggplot2::theme_classic(base_size = 16)+
  ggplot2::geom_hline(yintercept=(-log10(0.05)), linetype="dashed", color = "red")+
  ggplot2::geom_vline(xintercept=2, linetype="dashed", color = "red")+
  ggplot2::geom_vline(xintercept=-2, linetype="dashed", color = "red")+
  ggplot2::scale_y_continuous(limits=c(0,max(plotdf$y)))+
  geom_point(data=plotdf[grepl("SRS|SRRM|SF3B|SON",filtered1$gene.name),], 
             color='red',
             size=2)

# Define export parameters
ty <- "pdf" # File type
wi <- 500 # Width in pixels
he <- 500 # Height in pixels
re <- 75 # Resolution in dpi
mte <- 10 # Main title text size
# Create and export graph
Cairo(file=paste0("volcano_",group1,"v",group2), type=ty, width=wi, height=he, pointsize=mte, dpi=re)
print(g)
dev.off()
