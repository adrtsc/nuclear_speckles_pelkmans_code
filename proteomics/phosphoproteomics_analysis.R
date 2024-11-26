library(Cairo)
library(lmtest)
library(dplyr)
library(ggplot2)
library(eulerr) #for Venn diagrams
library(stringr) #to count R/S residues

set.seed(417)

setwd("") #set working directory in which output table and figures will be saved

#provide input files:
peptides.fi <- "peptides.txt" 
proteingroups.fi <- "proteinGroups.txt"
phosphopeptides.fi <- "Phospho\ (STY)Sites.txt"
speckleproteins.fi <- "QuickGO_nuclear_speck_GO_0016607.tsv"
nucleusproteins.fi <- "QuickGO_nucleus_GO_0005634.tsv"
phosphosite.db.fi <- "phosphoSite_phosphorylation_sites.txt"

filter.proteomics <- function(df){
  df <- df[df$Reverse != '+',]
  df <- df[df$Potential.contaminant != '+',]
  return(df)
}

log2.intensities <- function(df){
  df <- df[,colnames(df)[grep("^Intensity[:.:]",colnames(df))]]
  #keep only intensities summed across phosphorylation states for each peptide
  df <- df[,!grepl("___",colnames(df))] 
  df <- log2(as.matrix(df))
  df[is.infinite(df)] <- NA
  return(df)
}

filter.by.detection <- function(df){
  df.logical <- !is.na(df) #make logical dataframe for valid (non NA) entries
  #ensure measured in at least one condition
  results <- apply(df.logical,1,function(x){
    groups <- unique(gsub("_[0-9].[cP]","",names(x)))
    tempsum <- numeric()
    for (group in 1:length(groups)){
      tempsum[group] <- sum(x[grep(paste0("^",groups[group]),names(x))])
    }
    if (sum(grepl("[2,3]",tempsum)>0)) {y=T} else {y=F}
    return(y)
  })
  return(results)
}

impute.data <- function(df){
  df <- apply(df[,],2, function(x){
    x <- as.numeric(x)
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
  return(df)
}

run.lrt <- function(df,GNetc){
  condition <- gsub("_\\d.[cP]", "", colnames(df))
  phosphorylation <- sub('.*(?=.$)', '', colnames(df), perl=T)
  lrt <- sapply(1:nrow(df), function(i){ 
    peptide <- cbind(as.data.frame(as.matrix(df[i,])), condition,phosphorylation)
    colnames(peptide)[1] <- 'intensity'
    peptide$intensity <- as.numeric(peptide$intensity)
    null.model <- glm(intensity ~ condition + phosphorylation,data=peptide) #,family="Gamma")
    interaction.model <- glm(intensity ~ condition + phosphorylation + condition:phosphorylation,data=peptide) #,family="Gamma")
    lrt.results <- lrtest(interaction.model,null.model)
    pval <- lrt.results$`Pr(>Chisq)`[2]
    log2fc <- (mean(peptide$intensity[1:3])-mean(peptide$intensity[4:6]))-(mean(peptide$intensity[7:9])-mean(peptide$intensity[10:12]))
    out <- c(log2fc,pval)
    return(out)
  })
  lrt <- t(data.frame(lrt))
  lrt <- as.data.frame(cbind(GNetc,lrt))
  colnames(lrt) <- c("Gene.names","Sequence","Amino.acid","Position","log2fc","pval")
  rownames(lrt) <- NULL
  lrt$pval <- as.numeric(lrt$pval)
  lrt$log2fc <- as.numeric(lrt$log2fc)

  return(lrt)
}

plot.volcano <- function(df,group1,group2,geneDoe){
  plotdf <- data.frame(df[,c("log2fc","pval")])
  plotdf$pval <- -log10(plotdf$pval)
  colnames(plotdf) <- c("x","y")
  
  g <-  ggplot2::ggplot(plotdf, ggplot2::aes(x, y),echo=T) +
    ggplot2::geom_point(colour=densCols(plotdf$x,plotdf$y,colramp=colorRampPalette(c("grey","black"))),size=2)+
    ggplot2::theme_bw()+
    ggplot2::xlab(paste0("Difference in Phos vs. Unphos\nLog2FC (",group1,"/",group2,")")) +
    ggplot2::ylab("-log[10] p-value") +
    ggplot2::theme(plot.title = ggplot2::element_text(size=20, face="bold", vjust=2))+
    #ggrepel::geom_text_repel(data=plotdf[grepl(geneDoe,df$Gene.names),],ggplot2::aes(label = df$Gene.names[grepl(geneDoe,df$Gene.names)]),max.overlaps=Inf) + 
    ggplot2::theme_classic(base_size = 16)+
    ggplot2::geom_hline(yintercept=(-log10(0.05)), linetype="dashed", color = "red")+
    ggplot2::geom_vline(xintercept=2, linetype="dashed", color = "red")+
    ggplot2::geom_vline(xintercept=-2, linetype="dashed", color = "red")+
    ggplot2::scale_y_continuous(limits=c(0,max(plotdf$y))) +
    geom_point(data=plotdf[grepl(geneDoe,df$Gene.names),], 
               color='red',
               size=2)
  
  # Create and export graph
  Cairo(file=paste0("volcano_lrt_",group1,"v",group2,"_",geneDoe,".pdf"), type='pdf', width=500, height=500, pointsize=10, dpi=75)
  print(g)
  dev.off()
}

#Processing and imputation separately for pulldown of GFP-DYRK3 and GFP-PP1(m)-NIPP1, 
#but together for phosphorylated and unphosphorylated peptides

#Unphosphorylated peptides
peptides <- read.delim(peptides.fi, stringsAsFactors=FALSE)
peptides <- filter.proteomics(peptides)
unphos.intensities <- log2.intensities(peptides)
unphos.intensities <- data.frame(cbind(unphos.intensities,peptides$Sequence))
colnames(unphos.intensities) <- c(paste0(colnames(unphos.intensities)[1:dim(unphos.intensities)[2]-1],".c"),"Sequence")

#Phosphorylated peptides
pSTY.sites <- read.delim(phosphopeptides.fi, stringsAsFactors=FALSE)
pSTY.sites <- filter.proteomics(pSTY.sites)
#parse fasta column for gene names, add column with peptide sequence
pSTY.sites$Gene.names <- sapply(strsplit(as.character(sapply(strsplit(pSTY.sites$Fasta.headers, "GN="), tail, 1))," "),head,1)
pSTY.sites$Sequence <- gsub("\\s*\\([^\\)]+\\)","",as.character(pSTY.sites$Phospho..STY..Probabilities))
phos.intensities.log2 <- log2.intensities(pSTY.sites)
phos.intensities <- data.frame(cbind(phos.intensities.log2,pSTY.sites[,c("Sequence","Gene.names","Amino.acid","Position")]))
colnames(phos.intensities) <- c(paste0(colnames(phos.intensities.log2),".P"),"Sequence","Gene.names","Amino.acid","Position")

intensities <- inner_join(phos.intensities,unphos.intensities,by='Sequence',copy=TRUE)
dyrk3.l2fc <- c()
pp1.l2fc <- c()
dyrk3.p <- c()
pp1.p <- c()

#repeat imputation and statistical calculations 50x to then take medians
for(i in 1:50){
  #for comparison to expression of unphosphorylated peptide
  dyrk3.intensities <- intensities[,grepl("DYRK3",colnames(intensities))]
  dyrk3.detection.mask <- filter.by.detection(dyrk3.intensities)
  dyrk3.intensities <- dyrk3.intensities[dyrk3.detection.mask,]
  dyrk3.intensities.imputed <- impute.data(dyrk3.intensities)
  
  pp1.intensities <- intensities[,grepl("PP1",colnames(intensities))]
  pp1.detection.mask <- filter.by.detection(pp1.intensities)
  pp1.intensities <- pp1.intensities[pp1.detection.mask,]
  pp1.intensities.imputed <- impute.data(pp1.intensities)
  
  #for comparison to expression of protein
  #intensities <- left_join(phos.filtered,unphos.filtered,by='Gene.names')
  extra_cols <- as.data.frame(intensities)[,c("Gene.names","Sequence","Amino.acid","Position")]
  
  #fit GLMs with and without interaction terms for each peptide and run likelihood ratio test
  dyrk3.lrt.results <- run.lrt(dyrk3.intensities.imputed,extra_cols[dyrk3.detection.mask,])
  pp1.lrt.results <- run.lrt(pp1.intensities.imputed,extra_cols[pp1.detection.mask,])
  
  dyrk3.p <- cbind(dyrk3.p,dyrk3.lrt.results$pval)
  pp1.p <- cbind(pp1.p,pp1.lrt.results$pval)
  dyrk3.l2fc <- cbind(dyrk3.l2fc,dyrk3.lrt.results$log2fc)
  pp1.l2fc <- cbind(pp1.l2fc,pp1.lrt.results$log2fc)
}

dyrk3.lrt <- cbind(dyrk3.lrt.results[,c("Gene.names","Sequence","Amino.acid","Position")],
                   apply(dyrk3.l2fc, 1, median, na.rm=T),apply(dyrk3.p, 1, median, na.rm = T))
pp1.lrt <- cbind(pp1.lrt.results[,c("Gene.names","Sequence","Amino.acid","Position")],
                 apply(pp1.l2fc, 1, median, na.rm=T),apply(pp1.p, 1, median, na.rm = T))
dyrk3.lrt$seq_pos <- paste0(dyrk3.lrt$Sequence,"_",dyrk3.lrt$Position)
pp1.lrt$seq_pos <- paste0(pp1.lrt$Sequence,"_",pp1.lrt$Position)
colnames(dyrk3.lrt) <- c("Gene.names","Sequence","Amino.acid","Position","log2fc","pval","seq_pos")
colnames(pp1.lrt) <- c("Gene.names","Sequence","Amino.acid","Position","log2fc","pval","seq_pos")

#Summarize results
print.summary <- function(df,lrt){
  condition <- gsub("_\\d.[cP]", "", colnames(df))
  condition1 <- gsub("^.*?\\.","",condition[1]) 
  print(paste('A total of',length(lrt$Sequence),'peptides were analyzed, of which',length(unique(lrt$Sequence)),'were unique'))
  print(paste(dim(lrt[which(lrt$pval < 0.05 & abs(lrt$log2fc) >= 2),])[1],'sites/(non-unique) peptides were significant (p < 0.05) with difference in log2 fold change >= 2'))
  print(paste('Of these,',dim(lrt[which(lrt$pval < 0.05 & lrt$log2fc >= 2),])[1],'showed greater phosphorylation in',condition1))
  print(paste(length(unique(lrt[which(lrt$pval < 0.05 & abs(lrt$log2fc) >= 2),"Sequence"])),'unique peptides were significant (p < 0.05) with difference in log2 fold change >= 2'))
  print(paste('Of these,',length(unique(lrt[which(lrt$pval < 0.05 & lrt$log2fc >= 2),"Sequence"])),'showed greater phosphorylation in',condition1))
}
print.summary(dyrk3.intensities.imputed,dyrk3.lrt)
print.summary(pp1.intensities.imputed,pp1.lrt)

#There may be multiple entries per sequence where different sequence windows were considered for phosphorylation analysis
#But these correspond to different sites
#Results may differ depending on the imputation therefore repeat imputation and statistical calculations 50x
print(paste(length(intersect(dyrk3.lrt$seq_pos,pp1.lrt$seq_pos)),"peptides detected with DYRK3 & PP1 pulldown"))
#consider only peptides detected in both pulldowns 
p.sig.mask <- which(pp1.lrt$pval < 0.05 & pp1.lrt$log2fc >= 2 & pp1.lrt$seq_pos %in% intersect(dyrk3.lrt$seq_pos,pp1.lrt$seq_pos))
p.sig <- unique(pp1.lrt$seq_pos[p.sig.mask])
d.sig.mask <- which(dyrk3.lrt$pval < 0.05 & dyrk3.lrt$log2fc >= 2 & dyrk3.lrt$seq_pos %in% intersect(dyrk3.lrt$seq_pos,pp1.lrt$seq_pos))
d.sig <- unique(dyrk3.lrt$seq_pos[d.sig.mask])
pd.overlap <- length(intersect(p.sig,d.sig))
print(paste(pd.overlap,"significant peptides detected with both pulldowns"))

g <- plot(euler(c("PP1" = length(p.sig)-pd.overlap, "DYRK3" = length(d.sig)-pd.overlap, "PP1&DYRK3" = pd.overlap)), 
     fills = list(fill = c("#2626AC", "#AE7100"), alpha = 0.5),quantities = TRUE)
Cairo(file=paste0("venn_PP1_DYRK3_phosdiff_peptidelevel.pdf"), type='pdf', width=250, height=250, pointsize=10, dpi=75)
print(g)
dev.off()

#plot Venn for proteins-level, rather than peptide-level overlap
p.sig <- unique(pp1.lrt$Gene.names[p.sig.mask])
d.sig <- unique(dyrk3.lrt$Gene.names[d.sig.mask])
pd.overlap <- length(intersect(p.sig,d.sig))

g <- plot(euler(c("PP1" = length(p.sig)-pd.overlap, "DYRK3" = length(d.sig)-pd.overlap, "PP1&DYRK3" = pd.overlap)), 
          fills = list(fill = c("#2626AC", "#AE7100"), alpha = 0.5),quantities = TRUE)
Cairo(file=paste0("venn_PP1_DYRK3_phosdiff_proteinlevel.pdf"), type='pdf', width=250, height=250, pointsize=10, dpi=75)
print(g)
dev.off()

#plot Venn for Nuclear speckle proteins specifically
speckle.df = read.delim(speckleproteins.fi,sep='\t')
speckle.proteins <- speckle.df$SYMBOL
p.sig.mask <- which(pp1.lrt$pval < 0.05 & pp1.lrt$log2fc >= 2 & 
                      pp1.lrt$seq_pos %in% intersect(dyrk3.lrt$seq_pos,pp1.lrt$seq_pos) &
                      pp1.lrt$Gene.names %in% speckle.proteins)
p.sig <- unique(pp1.lrt$seq_pos[p.sig.mask])
d.sig.mask <- which(dyrk3.lrt$pval < 0.05 & dyrk3.lrt$log2fc >= 2 & 
                      dyrk3.lrt$seq_pos %in% intersect(dyrk3.lrt$seq_pos,pp1.lrt$seq_pos) &
                      dyrk3.lrt$Gene.names %in% speckle.proteins)
d.sig <- unique(dyrk3.lrt$seq_pos[d.sig.mask])
pd.overlap <- length(intersect(p.sig,d.sig))

g <- plot(euler(c("PP1" = length(p.sig)-pd.overlap, "DYRK3" = length(d.sig)-pd.overlap, "PP1&DYRK3" = pd.overlap)), 
          fills = list(fill = c("#2626AC", "#AE7100"), alpha = 0.5),quantities = TRUE)
Cairo(file=paste0("venn_PP1_DYRK3_phosdiff_peptidelevel_speckles.pdf"), type='pdf', width=250, height=250, pointsize=10, dpi=75)
print(g)
dev.off()

#plot Venn for proteins-level, rather than peptide-level overlap (speckle-specific)
p.sig <- unique(pp1.lrt$Gene.names[p.sig.mask])
d.sig <- unique(dyrk3.lrt$Gene.names[d.sig.mask])
pd.overlap <- length(intersect(p.sig,d.sig))

g <- plot(euler(c("PP1" = length(p.sig)-pd.overlap, "DYRK3" = length(d.sig)-pd.overlap, "PP1&DYRK3" = pd.overlap)), 
          fills = list(fill = c("#2626AC", "#AE7100"), alpha = 0.5),quantities = TRUE)
Cairo(file=paste0("venn_PP1_DYRK3_phosdiff_proteinlevel_speckles.pdf"), type='pdf', width=250, height=250, pointsize=10, dpi=75)
print(g)
dev.off()

#plot Venn for Nuclear proteins specifically
nucleus.df = read.delim(nucleusproteins.fi,sep='\t')
nucleus.proteins <- nucleus.df$SYMBOL
p.sig.mask <- which(pp1.lrt$pval < 0.05 & pp1.lrt$log2fc >= 2 & 
                      pp1.lrt$seq_pos %in% intersect(dyrk3.lrt$seq_pos,pp1.lrt$seq_pos) &
                      pp1.lrt$Gene.names %in% nucleus.proteins)
p.sig <- unique(pp1.lrt$seq_pos[p.sig.mask])
d.sig.mask <- which(dyrk3.lrt$pval < 0.05 & dyrk3.lrt$log2fc >= 2 & 
                      dyrk3.lrt$seq_pos %in% intersect(dyrk3.lrt$seq_pos,pp1.lrt$seq_pos) &
                      dyrk3.lrt$Gene.names %in% nucleus.proteins)
d.sig <- unique(dyrk3.lrt$seq_pos[d.sig.mask])
pd.overlap <- length(intersect(p.sig,d.sig))

g <- plot(euler(c("PP1" = length(p.sig)-pd.overlap, "DYRK3" = length(d.sig)-pd.overlap, "PP1&DYRK3" = pd.overlap)), 
          fills = list(fill = c("#2626AC", "#AE7100"), alpha = 0.5),quantities = TRUE)
Cairo(file=paste0("venn_PP1_DYRK3_phosdiff_peptidelevel_nucleus.pdf"), type='pdf', width=250, height=250, pointsize=10, dpi=75)
print(g)
dev.off()

#plot Venn for proteins-level, rather than peptide-level overlap (nucleus-specific)
p.sig <- unique(pp1.lrt$Gene.names[p.sig.mask])
d.sig <- unique(dyrk3.lrt$Gene.names[d.sig.mask])
pd.overlap <- length(intersect(p.sig,d.sig))

g <- plot(euler(c("PP1" = length(p.sig)-pd.overlap, "DYRK3" = length(d.sig)-pd.overlap, "PP1&DYRK3" = pd.overlap)), 
          fills = list(fill = c("#2626AC", "#AE7100"), alpha = 0.5),quantities = TRUE)
Cairo(file=paste0("venn_PP1_DYRK3_phosdiff_proteinlevel_nucleus.pdf"), type='pdf', width=250, height=250, pointsize=10, dpi=75)
print(g)
dev.off()

#Save results for supplementary table
dyrk3.lrt$KinPhos <- "DYRK3"
pp1.lrt$KinPhos <- "PP1-NIPP1"
lrt <- rbind(dyrk3.lrt,pp1.lrt)
lrt <- lrt[,c("Gene.names","Sequence","Amino.acid","Position","log2fc","pval","KinPhos")]
write.table(lrt,"SuppTable1_Phosphoproteomics.csv",row.names = FALSE,sep=',',quote=FALSE)

#plot volcanoes 
plot.volcano(dyrk3.lrt,"DYRK3","DYRK3+GSK","DYRK3")
plot.volcano(pp1.lrt,"PP1m-NIPP1","PP1-NIPP1","DYRK3")
plot.volcano(dyrk3.lrt,"DYRK3","DYRK3+GSK","SRRM2")
plot.volcano(pp1.lrt,"PP1m-NIPP1","PP1-NIPP1","SRRM2")
plot.volcano(dyrk3.lrt,"DYRK3","DYRK3+GSK","THRAP3")
plot.volcano(pp1.lrt,"PP1m-NIPP1","PP1-NIPP1","THRAP3")
plot.volcano(dyrk3.lrt,"DYRK3","DYRK3+GSK","SRRM1")
plot.volcano(pp1.lrt,"PP1m-NIPP1","PP1-NIPP1","SRRM1")
plot.volcano(dyrk3.lrt,"DYRK3","DYRK3+GSK","SON")
plot.volcano(pp1.lrt,"PP1m-NIPP1","PP1-NIPP1","SON")

srrm2.intersect <- unique(pp1.lrt[which((pp1.lrt$Gene.names == "SRRM2") & (pp1.lrt$seq_pos %in% intersect(p.sig,d.sig))),])
print(paste('# SRRM2 sites =',dim(srrm2.intersect)[1]))
intersect.genes <- pp1.lrt$Gene.names[which(pp1.lrt$seq_pos %in% intersect(p.sig,d.sig))]
intersect.genes.any.pep <- intersect(pp1.lrt$Gene.names[p.sig.mask],dyrk3.lrt$Gene.names[d.sig.mask])
srrm2.intersect.positions <- srrm2.intersect[,c("Gene.names","Amino.acid","Position")]

#lollipop plot for phosphosites with GFP-DYRK3 vs. GFP-DYRK3+GSK626616 and GFP-PP1m-NIPP1 vs. GFP-PP1-NIPP1
plot.lollipop <- function(gene,genelen,combine=FALSE){
  sites.p <- cbind(pp1.lrt[which(pp1.lrt$Gene.names == gene),],"PP1")
  colnames(sites.p) <- c(colnames(pp1.lrt),"pulldown")
  print(paste(dim(sites.p[which((sites.p$log2fc >= 2 & sites.p$pval < 0.05)),])[1],"sites in",gene,"sig for PP1"))
  if(combine){
    sites.d <- cbind(dyrk3.lrt[which(dyrk3.lrt$Gene.names == gene),],"DYRK3")
    colnames(sites.d) <- c(colnames(dyrk3.lrt),"pulldown")
    sites <- rbind(sites.d,sites.p)
    print(paste(dim(sites.d[which((sites.d$log2fc >= 2 & sites.d$pval < 0.05)),])[1],"sites in",gene,"sig for DYRK3"))
  }else{ sites <- sites.p}
  sites$p <- -log10(sites$pval)
  print(paste('max -log10(p)',max(sites$p)))
  sites$p <- sites$p/14.4 #max(sites$p) all normalized to SRRM2 max
  sites$col <- "#2626AC"
  sites[which(sites$pulldown == "DYRK3"),"col"] <- "#AE7100"
  g <- ggplot(sites, aes(x=Position, y=log2fc)) +
    geom_segment( aes(x=Position, xend=Position, y=0, yend=log2fc)) +
    geom_point( size=1, color=sites$col, fill=alpha(sites$col, sites$p), 
                alpha=sites$p, shape=21, stroke=1) +
    geom_hline(yintercept = 0,color="black") + coord_cartesian(xlim = c(0, genelen),ylim = c(-11,14),expand=F) + 
    xlab(gene) +
    theme_bw()
  Cairo(file=paste0(gene,"_phosphosites_PP1.pdf"), type='pdf', width=400, height=150, pointsize=10, dpi=75, transparent=TRUE)
  print(g)
  dev.off()
}
len.son <- 2426
len.srrm2 <- 2752
len.srsf10 <- 262
len.srsf11 <- 484
len.srsf6 <- 344
len.srsf8 <- 282
plot.lollipop("SON",len.son)
plot.lollipop("SRRM2",len.srrm2)
plot.lollipop("SRSF10",len.srsf10)
plot.lollipop("SRSF11",len.srsf11)
plot.lollipop("SRSF6",len.srsf6)
plot.lollipop("SRSF8",len.srsf8)
genelen <- 2752

#plot RS enrichment of SRRM2
srrm2.seq <- read.delim("SRRM2_uniprot.fasta",sep="\n")
srrm2.seq <- toString(paste0(unlist(srrm2.seq),collapse=''))
window.size = 100
r.enrichment <- c()
s.enrichment <- c()
for (i in 1:(nchar(srrm2.seq)-window.size)){
  x <- substr(srrm2.seq,i,i+50)
  r.enrichment = c(r.enrichment,str_count(x, "R")/50)
  s.enrichment = c(s.enrichment,str_count(x, "S")/50)
}

rs.df <- data.frame("RS"=c(r.enrichment,s.enrichment),"AA"=c(rep("R",length(r.enrichment)),rep("S",length(r.enrichment))),
                    "SRRM2"=rep(1:length(r.enrichment),2))
rs.df["SRRM2"] <- as.numeric(rs.df$SRRM2)
rs.df["R"] <- as.numeric(rs.df$R)
g <- ggplot(rs.df, aes(x=SRRM2, y=RS, color=AA)) +
  geom_line() + scale_color_manual(values=c("#D81E5B","#331832")) + 
  coord_cartesian(xlim = c(0, len.srrm2),expand=FALSE)+theme_bw()+
  theme(rect = element_rect(fill = "transparent"))
Cairo(file=paste0("SRRM2_RS.pdf"), type='pdf', width=400, height=100, pointsize=10, dpi=75, bg="transparent")
print(g)
dev.off()

#plot protein detection vs. number of phosphosites detected 
#compared to number reported in PhosphoSite Plus

phosphosite.plus <- read.delim(phosphosite.db.fi, stringsAsFactors=FALSE)
phosphosite.plus <- phosphosite.plus[phosphosite.plus$ORGANISM == "human",]
phosphosite.plus.counts <- phosphosite.plus %>% count(GENE)
colnames(phosphosite.plus.counts) <- c("Gene.names","reported_phosphosites")

genes.of.interest <- c("SRRM2","SON")
protein.levels <- read.delim(proteingroups.fi)
pSTY.unique.sites <- distinct(pSTY.sites[,c("Gene.names","Amino.acid","Position")])
pSTY.counts <- pSTY.sites %>% count(Gene.names)

#better to impute protein expression levels as well? 
protein.levels <- filter.proteomics(protein.levels)
protein.levels$expression <- log10(rowMeans(protein.levels[,colnames(protein.levels)[grep("^LFQ.intensity.GFP.PP1.NIPP1",colnames(protein.levels))]])+1)
protein.levels$Gene.names <- sapply(strsplit(as.character(sapply(strsplit(as.character(protein.levels$Fasta.headers), "GN="), tail, 1))," "),head,1)
protein.levels <- protein.levels[,c("Gene.names","expression")]

proteins.df <- merge(protein.levels,pSTY.counts,by="Gene.names",all.x=TRUE)
proteins.df <- merge(proteins.df, phosphosite.plus.counts,by="Gene.names",all.x=TRUE)
proteins.df[is.na(proteins.df)] <- 0
proteins.df$n <- proteins.df$n + 1
proteins.df$reported_phosphosites <- proteins.df$reported_phosphosites + 1

g <-  ggplot2::ggplot(proteins.df, ggplot2::aes(x=expression, y=log10(n)),echo=T) +
  ggplot2::geom_point(aes(colour=reported_phosphosites+1,alpha=0.8)) +
  scale_colour_gradient2(midpoint=2,low="yellow", high="#0C2389", mid="#1FFF93",trans="log10",space = "Lab",
                         guide=guide_colourbar(title="# Reported\nPhosphosites")) +
  ggplot2::theme_bw()+
  ggplot2::xlab("log10(Mean Protein Intensity+1)") +
  ggplot2::ylab("log10(# Phosphosites Detected+1)") +
  ggplot2::theme(plot.title = ggplot2::element_text(size=20, face="bold", vjust=2)) +
  ggrepel::geom_text_repel(data=proteins.df[grepl("SRRM2",proteins.df$Gene.names),c("expression","n")],
                           ggplot2::aes(label = "SRRM2"),min.segment.length = unit(0, 'lines'),max.overlaps=Inf) + 
  ggrepel::geom_text_repel(data=proteins.df[grepl("SON",proteins.df$Gene.names),c("expression","n")],
                           ggplot2::aes(label = "SON"),min.segment.length = unit(0, 'lines'),max.overlaps=Inf) + 
  ggplot2::theme_classic(base_size = 16) 

# Create and export graph
Cairo(file=paste0("protein_detection_vs_phosphosites_scatter.pdf"), type='pdf', width=500, height=500, pointsize=10, dpi=75)
print(g)
dev.off()

g <-  ggplot2::ggplot(proteins.df, ggplot2::aes(x=log2(n), y=log2(reported_phosphosites)),echo=T) +
  ggplot2::geom_point(aes(colour=reported_phosphosites+1,alpha=0.8)) +
  scale_colour_gradient2(midpoint=1.2,low="lightpink", high="#0C2389", mid="purple", trans="log10",space = "Lab",
                         guide=guide_colourbar(title="# Reported\nPhosphosites")) +
  ggplot2::theme_bw()+
  ggplot2::xlab("log2(# Phosphosites Detected+1)") +
  ggplot2::ylab("log2(# Phosphosites Reported+1)") +
  ggplot2::theme(plot.title = ggplot2::element_text(size=20, face="bold", vjust=2)) +
  ggrepel::geom_text_repel(data=proteins.df[grepl("SRRM2",proteins.df$Gene.names),c("n","reported_phosphosites")],
                           ggplot2::aes(label = "SRRM2"),min.segment.length = unit(0, 'lines'),max.overlaps=Inf) + 
  ggrepel::geom_text_repel(data=proteins.df[grepl("SON",proteins.df$Gene.names),c("n","reported_phosphosites")],
                           ggplot2::aes(label = "SON"),min.segment.length = unit(0, 'lines'),max.overlaps=Inf) + 
  ggplot2::theme_classic(base_size = 16) 

# Create and export graph
Cairo(file=paste0("phosphosites_detected_vs_expected_scatter.pdf"), type='pdf', 
      width=400, height=320, pointsize=10, dpi=75,bg="transparent")
print(g)
dev.off()
