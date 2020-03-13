###### Stability Code #######
#### Loading Functions ####

check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)>0){
    print(paste0("Installing ", paste(new.pkg,collapse = ", ")))
    install.packages(new.pkg, dependencies = TRUE)
  }
  res <- lapply(pkg,load.packages)
}

check.Bioconductor.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)>0){
    print(paste0("Installing ", paste(new.pkg,collapse = ", ")))
    BiocManager::install(new.pkg, dependencies = TRUE,ask=FALSE)
  }
  res <- lapply(pkg,load.packages)
}

load.packages <- function(pkg){
  print(paste0("Loading ",pkg))
  suppressMessages(require(pkg,character.only = T))
}

consensus.dataframe.res <- function(list.of.dfs,method="fdr",pval=0.01,quantile.value=0.25){
  res_pval <- p.adjust(unlist(lapply(rownames(list.of.dfs[[1]]),doFisher,list.of.dfs)),
                       method = method)
  res_fc <- unlist(mclapply(names(res_pval),check_fc,list.of.dfs,mc.cores = 5))
  df <- data.frame(RE=names(res_pval),
                   p.adjusted = res_pval,
                   fc=res_fc,
                   stringsAsFactors = F)
  df <- na.omit(df)
  df <- df[df$p.adjusted<pval,]
  fcs_p <- quantile(abs(df$fc),quantile.value)
  df <- df[abs(df$fc)>fcs_p,]
  return(df)
}

check_fc <- function(gene,list.of.dfs){
  fcvalues <- c()
  res <- c()
  for (i in 1:length(list.of.dfs)){
    df <- list.of.dfs[[i]]
    fcvalues <- c(fcvalues,sign(df[gene,"logFC"]))
  }
  fcvalues <- unique(fcvalues)
  if (length(fcvalues)==1){
    for (i in 1:length(list.of.dfs)){
      df <- list.of.dfs[[i]]
      res <- c(res,df[gene,"logFC"])
    }
    res <- mean(res)
  } else{
    res <- NA
  }
  names(res) <- gene
  return(res)
}

doFisher <- function(gene,list.of.dfs){
  pvalues <- c()
  for (i in 1:length(list.of.dfs)){
    df <- list.of.dfs[[i]]
    pvalues <- c(pvalues,df[gene,"P.Value"])
  }
  res <- maximump(pvalues)$p
  names(res) <- gene
  return(res)
}


extract_fgsea <- function(tf,fgseaRes){
  fgseaRes <- fgseaRes[fgseaRes$pathway == tf,]
  fgseaRes <- unlist(fgseaRes$leadingEdge)
  df <- data.frame(TF=rep(tf,length(fgseaRes)),
                   Target=fgseaRes,stringsAsFactors = F)
  return(df)
}

dolimma = function(data, input1, input2) {
  data = data[,c(input1, input2)]
  input1Name =  deparse(substitute(input1)) 
  input2Name =  deparse(substitute(input2)) 
  design = matrix(data=0,ncol=2,nrow=ncol(data),
                  dimnames = list(c(input1, input2), c(input1Name, input2Name)))
  design[1:length(input1),1] = 1
  design[(length(input1)+1):(length(input1) + length(input2)),2] = 1
  fit <- lmFit(data, design)
  contrast = paste0(input1Name, "-", input2Name)
  cont.matrix <- makeContrasts(contrasts = contrast, levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2, 0.05)
  limma <- topTable(fit2, adjust="holm", sort.by="B", resort.by = "p", number=dim(data)[1])
  return(limma)
}

run_fgsea <- function(expressionMatrix,cases,control,examplePathways){
  library(fgsea)
  diffExpression <- dolimma(data = expressionMatrix,input1 = cases,input2 = control)
  diffExpression <- diffExpression[order(diffExpression$t,decreasing = T),]
  exampleRanks <- diffExpression$t
  names(exampleRanks) <- rownames(diffExpression)
  set.seed(123)
  fgseaRes <- fgsea(pathways = examplePathways, 
                    stats = exampleRanks,
                    minSize=1,
                    maxSize=500,
                    nperm=100000)
  
  res <- list(fgseaRes=fgseaRes,diffExpression=diffExpression)
  return(res)
}

obtain_tf_gene_association <- function(tf,viper_regulon){
  regulon <- viper_regulon[[tf]]
  df = data.frame(TF=rep(tf,length(regulon$tfmode)),
                  Target=names(regulon$tfmode),
                  TFMode=as.numeric(regulon$tfmode),
                  stringsAsFactors = F)
  return(df)
}

dorothea_to_table <- function(viper_regulon){
  names(viper_regulon) = sapply(strsplit(names(viper_regulon), split = ' - '), head, 1)
  results <- lapply(names(viper_regulon),obtain_tf_gene_association,viper_regulon)
  df = results[[1]]
  for (i in 2:length(results)){
    df <- unique(rbind(df,results[[i]]))
  }
  return(df)
}

table_to_regulon <- function(dorothea_table){
  dorothea_regulon = list()
  for (tf in unique(dorothea_table$TF)){
    regulon <- dorothea_table[dorothea_table$TF==tf,]
    insert_list = list('tfmode'=regulon$TFMode,
                       'likelihood'=rep(1,nrow(regulon)))
    names(insert_list$tfmode) <- regulon$Target
    dorothea_regulon[[tf]] <- insert_list
  }
  return(dorothea_regulon)
}

extract_genes <- function(term,ego_res){
  ego_res <- ego_res[ego_res$Description == term,]
  genes <- list(strsplit(ego_res$geneID,"/"))
  return(genes)
}

extract_clue <- function(id,clue){
  targets <- clue[clue$ID==id,"Target"]
  targets <- unlist(strsplit(targets,", "))
  return(targets)
}

clue_table <- function(tfs){
  clue <- read.delim("clue.txt",stringsAsFactors = F)
  clue <- clue[,-c(1,2,3)]
  clue <- clue[clue$Target != "",]
  res <- lapply(clue$ID,extract_clue,clue)
  names(res) <- clue$ID
  df_clue <- clue[0,]
  for (id in names(res)){
    insert_df <- clue[clue$ID==id,]
    duptimes <- length(res[[id]])
    idx <- rep(1:nrow(insert_df), duptimes)
    insert_df <- insert_df[idx,]
    insert_df$Target <- res[[id]]
    df_clue <- rbind(df_clue,insert_df)
  }
  rownames(df_clue) <- NULL
  df_clue <- df_clue[df_clue$Target %in% tfs,]
  return(df_clue)
}

reduce_dimension <- function(tf,nodes,res){
  sub_res <- res[res$TF==tf,]
  sub_nodes_tf <- nodes[tf,]
  sub_nodes_gene <- nodes[nodes$id %in% sub_res$Target,]
  sub_nodes_gene <- inner_join(sub_nodes_gene,sub_res,
                               by = c("id"="Target"))
  
  if(sub_nodes_tf$diffExpression == 1){ # High Activity
    sub_nodes_gene <- sub_nodes_gene[sub_nodes_gene$diffExpression == sub_nodes_gene$TFMode,]
  } else{
    sub_nodes_gene <- sub_nodes_gene[sub_nodes_gene$diffExpression != sub_nodes_gene$TFMode,]
  }
  sub_res <- sub_res[sub_res$Target %in% sub_nodes_gene$id,]
  return(sub_res)
}

doWilcoxon <- function(TFMatrix,g1,g2){
  TFMatrix_ordered <- as.data.frame(t(TFMatrix[,c(g1,g2)]))
  TFMatrix_ordered$Group <- c(rep("g1",length(g1)),rep("g2",length(g2)))
  res <- lapply(rownames(TFMatrix),doWilcoxonByTF,TFMatrix_ordered)
  res <- do.call("rbind",res)
  return(res)
}

doWilcoxonByTF <- function(tf,TFMatrix){
  subTable <- TFMatrix[,c(tf,"Group")]
  colnames(subTable)[1] <- "TF"
  p <- wilcox.test(TF ~ Group,data = subTable)$p.value
  m <- mean(subTable[subTable$Group=="g1","TF"]) - mean(subTable[subTable$Group=="g2","TF"])
  res <- data.frame(P.Value=p,logFC=m,row.names = tf)
  return(res)
}


#### Loading packages ####

check.Bioconductor.packages(c("ConsensusClusterPlus","gprofiler2",'limma','viper',"clusterProfiler","rWikiPathways","org.Hs.eg.db","ReactomePA"))
check.packages(c("clusteval","scico","wesanderson","mclust",'pheatmap','ggpubr','reshape2','ggplot2',"RColorBrewer",'metap','grid','gridExtra','rlist','dplyr',"randomcoloR"))


load("expression_cluster/toload.RData")

load("A_viperRegulon.rdata")
A_viper <- dorothea_to_table(viper_regulon)
dorothea_table <- unique(rbind(A_viper))
dorothea_table <- dorothea_table[dorothea_table$Target %in% rownames(expressionAdult),]
dorothea_regulon <- table_to_regulon(dorothea_table)
rownames(dorothea_table) <- paste0(dorothea_table$TF,"_",dorothea_table$Target)



TFActivitiesPediatric <- viper(expressionPediatric,dorothea_regulon,minsize = 1)
TFActivitiesAdult <- viper(expressionAdult,dorothea_regulon,minsize = 1)

diffActivity.pediatric.cl1 <- doWilcoxon(TFActivitiesPediatric,pediatric.cl1,pediatric.hc)
diffActivity.pediatric.cl2 <- doWilcoxon(TFActivitiesPediatric,pediatric.cl2,pediatric.hc)
diffActivity.adult.cl1 <- doWilcoxon(TFActivitiesAdult,adult.cl1,adult.hc)
diffActivity.adult.cl2 <- doWilcoxon(TFActivitiesAdult,adult.cl2,adult.hc)


####### All SLE #####
list.of.dfs <- list(diffActivity.adult.cl1,diffActivity.adult.cl2,
                    diffActivity.pediatric.cl1,diffActivity.pediatric.cl2)
res <- consensus.dataframe.res(list.of.dfs = list.of.dfs,method = "bonferroni",pval=0.05,quantile.value = 0)
res.meta.tfs.consensus <- res[order(res$fc),]
tfs.consensus <- rownames(res.meta.tfs.consensus)
tfs <- tfs.consensus

write.table(file = "expression_cluster/TFs.txt",
            tfs,row.names = F,col.names = F,quote = F)


outPediatric <- pheatmap(TFActivitiesPediatric[c(tfs.consensus),c(pediatric.hc,pediatric.cl1,pediatric.cl2)],
                         scale = "row",breaks = seq(-2,2,length.out = 100),
                         cluster_rows = F,cluster_cols = F,show_colnames = F,
                         show_rownames = F,legend = F,annotation_col = anncolPediatric[c(pediatric.hc,pediatric.cl1,pediatric.cl2),"cluster",drop=F],
                         annotation_colors = list("cluster"=c("Healthy"="#A5FF66",
                                                              "cluster1"="#FF7C4C",
                                                              "cluster2"="#293ECC")),
                         annotation_names_col = F,annotation_legend = F,
                         border_color = NA,fontsize = 5)
outAdult <- pheatmap(TFActivitiesAdult[c(tfs.consensus),c(adult.hc,adult.cl1,adult.cl2)],
                     scale = "row",breaks = seq(-2,2,length.out = 100),
                     cluster_rows = F,cluster_cols = F,show_colnames = F,
                     annotation_col = anncolAdult[c(adult.hc,adult.cl1,adult.cl2),"cluster",drop=F],
                     annotation_colors = list("cluster"=c("Healthy"="#A5FF66",
                                                          "cluster1"="#FF7C4C",
                                                          "cluster2"="#293ECC")),
                     annotation_names_col = F,
                     border_color = NA,fontsize = 5)
items=c("outPediatric","outAdult")
plot_list=list()
for (a in items){
  x=get(a)
  plot_list[[a]] = x[[4]]     ##to save each plot into a list. note the [[4]]
}
F1A<-grid.arrange(arrangeGrob(grobs= plot_list,ncol=2))
ggsave("expression_cluster/Figure_2A.pdf",F1A,width = 6,height = 6)
dev.off()
ggsave("expression_cluster/Figure_2A.png",F1A,width = 8,height = 5)


df_clue <- clue_table(tfs.consensus)


examplePathways <- list()
for (i in 1:length(tfs)){
  tf <- tfs[i]
  examplePathways <- list.append(examplePathways,
                                 tf=names(dorothea_regulon[[tf]]$tfmode))
  names(examplePathways)[i] <- tf
}

fgseaRes.cl1.pediatric <- run_fgsea(expressionPediatric,pediatric.cl1,pediatric.hc,examplePathways)
diffexpression_pediatric.cl1 <- fgseaRes.cl1.pediatric$diffExpression
fgseaRes.cl1.pediatric <- fgseaRes.cl1.pediatric$fgseaRes
fgseaRes.cl2.pediatric <- run_fgsea(expressionPediatric,pediatric.cl2,pediatric.hc,examplePathways)
diffexpression_pediatric.cl2 <- fgseaRes.cl2.pediatric$diffExpression
fgseaRes.cl2.pediatric <- fgseaRes.cl2.pediatric$fgseaRes

fgseaRes.cl1.adult <- run_fgsea(expressionAdult,adult.cl1,adult.hc,examplePathways)
diffexpression_adult.cl1 <- fgseaRes.cl1.adult$diffExpression
fgseaRes.cl1.adult <- fgseaRes.cl1.adult$fgseaRes
fgseaRes.cl2.adult <- run_fgsea(expressionAdult,adult.cl2,adult.hc,examplePathways)
diffexpression_adult.cl2 <- fgseaRes.cl2.adult$diffExpression
fgseaRes.cl2.adult <- fgseaRes.cl2.adult$fgseaRes



res <- lapply(tfs[tfs %in% fgseaRes.cl1.pediatric$pathway],extract_fgsea,fgseaRes=fgseaRes.cl1.pediatric)
res_pediatric.cl1 <- do.call("rbind", res)
rownames(res_pediatric.cl1) <- paste0(res_pediatric.cl1$TF,"_",res_pediatric.cl1$Target)
res <- lapply(tfs[tfs %in% fgseaRes.cl2.pediatric$pathway],extract_fgsea,fgseaRes=fgseaRes.cl2.pediatric)
res_pediatric.cl2 <- do.call("rbind", res)
rownames(res_pediatric.cl2) <- paste0(res_pediatric.cl2$TF,"_",res_pediatric.cl2$Target)
res <- lapply(tfs[tfs %in% fgseaRes.cl1.adult$pathway],extract_fgsea,fgseaRes=fgseaRes.cl1.adult)
res_adult.cl1 <- do.call("rbind", res)
rownames(res_adult.cl1) <- paste0(res_adult.cl1$TF,"_",res_adult.cl1$Target)
res <- lapply(tfs[tfs %in% fgseaRes.cl2.adult$pathway],extract_fgsea,fgseaRes=fgseaRes.cl2.adult)
res_adult.cl2 <- do.call("rbind", res)
rownames(res_adult.cl2) <- paste0(res_adult.cl2$TF,"_",res_adult.cl2$Target)

keep_int <- c()
for (int in unique(c(rownames(res_pediatric.cl1),rownames(res_pediatric.cl2),
                     rownames(res_adult.cl1),rownames(res_adult.cl2)))){
  a <- res_pediatric.cl1[int,]
  b <- res_pediatric.cl2[int,]
  c <- res_adult.cl1[int,]
  d <- res_adult.cl2[int,]
  total <- na.omit(rbind(a,b,c,d))
  if (nrow(total)>=3){
    keep_int <- c(keep_int,int)
  }
}

res <- dorothea_table[rownames(dorothea_table) %in% keep_int,]

diffRes <- consensus.dataframe.res(list.of.dfs = list(diffexpression_adult.cl1,diffexpression_adult.cl2,diffexpression_pediatric.cl1,diffexpression_pediatric.cl2),
                                   method = "fdr",pval = 1.1,quantile.value = 0)

diffRes <- diffRes[order(diffRes$fc),]
res <- res[res$Target %in% rownames(diffRes),]
tfs.toplot <- tfs.consensus[tfs.consensus %in% res$TF]

nodes <- unique(c(res$TF,res$Target))
nodes <- data.frame(id=nodes,
                    type=ifelse(nodes %in% tfs,'TF','Gene'),
                    diffExpression=rep(NA,length(nodes)),
                    stringsAsFactors = F,row.names = nodes)

for (tf in tfs.toplot){
  nodes[tf,"diffExpression"] <- sign(res.meta.tfs.consensus[rownames(res.meta.tfs.consensus) == tf,"fc"])
}

for (gene in unique(res$Target)){
  if (!(gene %in% tfs.toplot)){
    direction <- sign(diffRes[gene,"fc"])
    nodes[gene,"diffExpression"] <- as.numeric(direction)
  }
}


res_network <- lapply(tfs.toplot,reduce_dimension,nodes,res)
res_network <- do.call("rbind", res_network)
nodes <- nodes[nodes$id %in% c(res_network$TF,res_network$Target),]

nodes$group <- ifelse(nodes[,"type"]=="TF",
                      ifelse(nodes[,"diffExpression"]==1,"tf_up","tf_down"),
                      ifelse(nodes[,"diffExpression"]==1,"gene_up","gene_down"))




res_network$source <- res_network$TF
res_network$target <- res_network$Target
res_network$interaction <- ifelse(res_network$TFMode==1,
                                  "activation","inhibition")

length(unique(res_network$Target))
length(unique(res_network$TF))

res_network <- res_network[res_network$Target !="T",]

tfs.toplot <- unique(res_network$TF)

calculate_mean_fold_change <- function(split.genes,diffRes){
  split.genes <- unlist(split.genes)
  fc <- mean(diffRes[split.genes,"fc"])
  return(fc)
}



annotate_function <- function(tf,res_network,diffRes){
  targets <- c(res_network[res_network$TF==tf,"Target"])
  targets_paste <- paste(c(res_network[res_network$TF==tf,"Target"]),collapse = " ")
  system(paste0("python3 gc4pro.py --genes ",targets_paste," --outdir ",tf))
  if (file.exists(paste0("./gc4Results/",tf,"/input1-GO_BP_Results-GeneCodis4.tsv"))){
    ego_res <- read.delim(paste0("./gc4Results/",tf,"/input1-GO_BP_Results-GeneCodis4.tsv"),
                                stringsAsFactors = F)
    genes <- strsplit(ego_res$Genes,",")
    fc <- unlist(lapply(genes,calculate_mean_fold_change,diffRes))
    ego_res$fcs <- fc
    return(ego_res)
  }
}




res <- mclapply(tfs.toplot,
                annotate_function,
                res_network = res_network,
                diffRes = diffRes,
                mc.cores = 6)

names(res) <- tfs.toplot
res <- res[which(!sapply(res, is.null))]

terms <- c()
for (i in 1:length(res)){
  subRes <- res[i]
  terms <- c(terms,subRes[[1]]$Terms)
}

terms <- unique(terms)

annTable <- data.frame(matrix(0,nrow = length(res),ncol = length(terms)))
annGenes <- data.frame(matrix(0,nrow = length(res),ncol = length(terms)))
rownames(annTable) <- rownames(annGenes) <- names(res)
colnames(annTable) <- colnames(annGenes) <- terms


for (i in 1:length(res)){
  subRes <- res[i]
  tf <- names(subRes)
  pvalues <- subRes[[1]]$Hyp.pValAdj
  fcs <- subRes[[1]]$fcs
  terms <- subRes[[1]]$Terms
  
  annTable[tf,terms] <- pvalues
  annGenes[tf,terms] <- fcs
}

annotation_table <- as.data.frame(t(annTable[rowSums(annTable != 0) > 0,]))
annotation_genes <- as.data.frame(t(annGenes[colnames(annotation_table),rownames(annotation_table)]))

annotation_df <- data.frame(REs=0,
                            Ann=0,
                            logFC=0,
                            pval=0)

annotation_df <- annotation_df[-1,]


for (term in rownames(annotation_table)){
  for (tf in colnames(annotation_table)){
    pvalues <- annotation_table[term,tf]
    pvalues <- pvalues[pvalues != 0]
    
    fcs <- annotation_genes[term,tf]
    if (length(pvalues)>0){
      n_row <-  data.frame(REs=tf,
                           Ann=term,
                           logFC=fcs,
                           pval=-log(pvalues))
      annotation_df <- rbind(annotation_df,n_row)
    }
  }
}


annotation_df <- annotation_df[annotation_df$pval > -log(0.001),]


out <- pheatmap(annotation_genes[unique(annotation_df$Ann),],clustering_method = "ward.D2")
term_order <- rownames(annotation_genes[unique(annotation_df$Ann),][out$tree_row[["order"]],])
gene_order <- colnames(annotation_genes[unique(annotation_df$Ann),][,out$tree_col[["order"]]])
annotation_table <- annotation_table[term_order,gene_order]

annotation_df$REs <- factor(annotation_df$REs,levels = colnames(annotation_table))
limit <- c(-max(abs(annotation_df$logFC)),max(abs(annotation_df$logFC)))


tfs <- ggplot(annotation_df, aes(REs, Ann,colour=logFC))+
  geom_point(aes(size=pval))+
  scale_colour_scico(palette = "vik", limit = limit)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size = 5),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "gray60"),
        panel.background = element_rect(fill = "white"))+
  xlab(label = "")+
  ylab(label = "")+
  labs(size = "-log(pval)")


consensus_tfs <- F1A
save(consensus_tfs,tfs,file = "expression_cluster/TFs_consensus.RData")
tfs
ggsave(filename = "expression_cluster/tfs.tiff",
       width = 10,height = 10)
dev.off()


tfs_consensus <- tfs.consensus


####### Cluster1 #####
list.of.dfs <- list(diffActivity.adult.cl1[ !(row.names(diffActivity.adult.cl1) %in% tfs_consensus), ],
                    diffActivity.pediatric.cl1[ !(row.names(diffActivity.pediatric.cl1) %in% tfs_consensus), ])
res <- consensus.dataframe.res(list.of.dfs = list.of.dfs,method = "bonferroni",pval=0.05,quantile.value = 0)
res.meta.tfs.consensus <- res[order(res$fc),]
tfs.consensus <- rownames(res.meta.tfs.consensus)
tfs <- tfs.consensus

write.table(file = "expression_cluster/TFs.txt",
            tfs,row.names = F,col.names = F,quote = F)


outPediatric <- pheatmap(TFActivitiesPediatric[c(tfs.consensus),c(pediatric.hc,pediatric.cl1,pediatric.cl2)],
                         scale = "row",breaks = seq(-2,2,length.out = 100),
                         cluster_rows = F,cluster_cols = F,show_colnames = F,
                         show_rownames = F,legend = F,annotation_col = anncolPediatric[c(pediatric.hc,pediatric.cl1,pediatric.cl2),"cluster",drop=F],
                         annotation_colors = list("cluster"=c("Healthy"="#A5FF66",
                                                              "cluster1"="#FF7C4C",
                                                              "cluster2"="#293ECC")),
                         annotation_names_col = F,annotation_legend = F,
                         border_color = NA,fontsize = 5)
outAdult <- pheatmap(TFActivitiesAdult[c(tfs.consensus),c(adult.hc,adult.cl1,adult.cl2)],
                     scale = "row",breaks = seq(-2,2,length.out = 100),
                     cluster_rows = F,cluster_cols = F,show_colnames = F,
                     annotation_col = anncolAdult[c(adult.hc,adult.cl1,adult.cl2),"cluster",drop=F],
                     annotation_colors = list("cluster"=c("Healthy"="#A5FF66",
                                                          "cluster1"="#FF7C4C",
                                                          "cluster2"="#293ECC")),
                     annotation_names_col = F,
                     border_color = NA,fontsize = 5)
items=c("outPediatric","outAdult")
plot_list=list()
for (a in items){
  x=get(a)
  plot_list[[a]] = x[[4]]     ##to save each plot into a list. note the [[4]]
}
F1A<-grid.arrange(arrangeGrob(grobs= plot_list,ncol=2))
ggsave("expression_cluster/Figure_2A.pdf",F1A,width = 6,height = 6)
dev.off()
ggsave("expression_cluster/Figure_2A.png",F1A,width = 8,height = 5)


df_clue <- clue_table(tfs.consensus)


examplePathways <- list()
for (i in 1:length(tfs)){
  tf <- tfs[i]
  examplePathways <- list.append(examplePathways,
                                 tf=names(dorothea_regulon[[tf]]$tfmode))
  names(examplePathways)[i] <- tf
}

fgseaRes.cl1.pediatric <- run_fgsea(expressionPediatric,pediatric.cl1,pediatric.hc,examplePathways)
diffexpression_pediatric.cl1 <- fgseaRes.cl1.pediatric$diffExpression
fgseaRes.cl1.pediatric <- fgseaRes.cl1.pediatric$fgseaRes

fgseaRes.cl1.adult <- run_fgsea(expressionAdult,adult.cl1,adult.hc,examplePathways)
diffexpression_adult.cl1 <- fgseaRes.cl1.adult$diffExpression
fgseaRes.cl1.adult <- fgseaRes.cl1.adult$fgseaRes


res <- lapply(tfs[tfs %in% fgseaRes.cl1.pediatric$pathway],extract_fgsea,fgseaRes=fgseaRes.cl1.pediatric)
res_pediatric.cl1 <- do.call("rbind", res)
rownames(res_pediatric.cl1) <- paste0(res_pediatric.cl1$TF,"_",res_pediatric.cl1$Target)
res <- lapply(tfs[tfs %in% fgseaRes.cl1.adult$pathway],extract_fgsea,fgseaRes=fgseaRes.cl1.adult)
res_adult.cl1 <- do.call("rbind", res)
rownames(res_adult.cl1) <- paste0(res_adult.cl1$TF,"_",res_adult.cl1$Target)

keep_int <- c()
for (int in unique(c(rownames(res_pediatric.cl1),
                     rownames(res_adult.cl1)))){
  a <- res_pediatric.cl1[int,]
  c <- res_adult.cl1[int,]
  total <- na.omit(rbind(a,c))
  if (nrow(total)>=2){
    keep_int <- c(keep_int,int)
  }
}

res <- dorothea_table[rownames(dorothea_table) %in% keep_int,]

diffRes <- consensus.dataframe.res(list.of.dfs = list(diffexpression_adult.cl1,diffexpression_pediatric.cl1),
                                   method = "fdr",pval = 1.1,quantile.value = 0)

diffRes <- diffRes[order(diffRes$fc),]
res <- res[res$Target %in% rownames(diffRes),]
tfs.toplot <- tfs.consensus[tfs.consensus %in% res$TF]

nodes <- unique(c(res$TF,res$Target))
nodes <- data.frame(id=nodes,
                    type=ifelse(nodes %in% tfs,'TF','Gene'),
                    diffExpression=rep(NA,length(nodes)),
                    stringsAsFactors = F,row.names = nodes)

for (tf in tfs.toplot){
  nodes[tf,"diffExpression"] <- sign(res.meta.tfs.consensus[rownames(res.meta.tfs.consensus) == tf,"fc"])
}

for (gene in unique(res$Target)){
  if (!(gene %in% tfs.toplot)){
    direction <- sign(diffRes[gene,"fc"])
    nodes[gene,"diffExpression"] <- as.numeric(direction)
  }
}


res_network <- lapply(tfs.toplot,reduce_dimension,nodes,res)
res_network <- do.call("rbind", res_network)
nodes <- nodes[nodes$id %in% c(res_network$TF,res_network$Target),]

nodes$group <- ifelse(nodes[,"type"]=="TF",
                      ifelse(nodes[,"diffExpression"]==1,"tf_up","tf_down"),
                      ifelse(nodes[,"diffExpression"]==1,"gene_up","gene_down"))




res_network$source <- res_network$TF
res_network$target <- res_network$Target
res_network$interaction <- ifelse(res_network$TFMode==1,
                                  "activation","inhibition")

length(unique(res_network$Target))
length(unique(res_network$TF))

res_network <- res_network[res_network$Target !="T",]

tfs.toplot <- unique(res_network$TF)

calculate_mean_fold_change <- function(split.genes,diffRes){
  split.genes <- unlist(split.genes)
  fc <- mean(diffRes[split.genes,"fc"])
  return(fc)
}

wikipathways_table <- read.delim("GO_BP_table.tsv",stringsAsFactors = F)
wikipathways_table <- wikipathways_table[grep("GC-9606",wikipathways_table$id),]
genes_table <- read.delim("gene_table.tsv",stringsAsFactors = F)
wikipathways_table <- merge(genes_table,wikipathways_table,on="id")[,c(2,5)]
annotation_table <- read.delim("annotation_info_table.tsv",stringsAsFactors = F)
wikipathways_table <- merge(wikipathways_table,annotation_table,on="annotation_id")


annotate_function <- function(tf,res_network,wikipathways_table,diffRes){
  targets <- c(res_network[res_network$TF==tf,"Target"])
  ewp.up <-enricher(
    targets,
    universe = unique(wikipathways_table$symbol),
    pAdjustMethod = "BH",
    pvalueCutoff = 1.1,
    qvalueCutoff = 1.1,
    TERM2GENE = wikipathways_table[,c(3,2)])
  if (!is.null(ewp.up)){
    ego_res <- ewp.up@result
    if (nrow(ego_res)>0){
      n_terms <- ego_res$Description
      pvalues <- ego_res$p.adjust
      genes <- strsplit(ego_res$geneID,"/")
      fc <- unlist(lapply(genes,calculate_mean_fold_change,diffRes))
      res <- list(pvalues=pvalues,
                  fcs = fc,
                  terms = n_terms)
      return(res)
    }
  }
}

res <- mclapply(tfs.toplot,
                annotate_function,
                res_network=res_network,
                wikipathways_table = wikipathways_table,
                diffRes = diffRes,
                mc.cores = 6)

names(res) <- tfs.toplot
res <- res[which(!sapply(res, is.null))]

terms <- c()
for (i in 1:length(res)){
  subRes <- res[i]
  terms <- c(terms,subRes[[1]]$terms)
}

terms <- unique(terms)

annTable <- data.frame(matrix(0,nrow = length(res),ncol = length(terms)))
annGenes <- data.frame(matrix(0,nrow = length(res),ncol = length(terms)))
rownames(annTable) <- rownames(annGenes) <- names(res)
colnames(annTable) <- colnames(annGenes) <- terms


for (i in 1:length(res)){
  subRes <- res[i]
  tf <- names(subRes)
  pvalues <- subRes[[1]]$pvalues
  fcs <- subRes[[1]]$fcs
  terms <- subRes[[1]]$terms
  
  annTable[tf,terms] <- pvalues
  annGenes[tf,terms] <- fcs
}

annotation_table <- as.data.frame(t(annTable[rowSums(annTable != 0) > 0,]))
annotation_genes <- as.data.frame(t(annGenes[colnames(annotation_table),rownames(annotation_table)]))

annotation_df <- data.frame(REs=0,
                            Ann=0,
                            logFC=0,
                            pval=0)

annotation_df <- annotation_df[-1,]


for (term in rownames(annotation_table)){
  for (tf in colnames(annotation_table)){
    pvalues <- annotation_table[term,tf]
    pvalues <- pvalues[pvalues != 0]
    
    fcs <- annotation_genes[term,tf]
    if (length(pvalues)>0){
      n_row <-  data.frame(REs=tf,
                           Ann=term,
                           logFC=fcs,
                           pval=-log(pvalues))
      annotation_df <- rbind(annotation_df,n_row)
    }
  }
}


annotation_df <- annotation_df[annotation_df$pval > -log(0.01) & abs(annotation_df$logFC)>0.5,]


out <- pheatmap(annotation_genes[unique(annotation_df$Ann),],clustering_method = "ward.D2")
term_order <- rownames(annotation_genes[unique(annotation_df$Ann),][out$tree_row[["order"]],])
gene_order <- colnames(annotation_genes[unique(annotation_df$Ann),][,out$tree_col[["order"]]])
annotation_table <- annotation_table[term_order,gene_order]

annotation_df$REs <- factor(annotation_df$REs,levels = colnames(annotation_table))
limit <- c(-max(abs(annotation_df$logFC)),max(abs(annotation_df$logFC)))


tfs <- ggplot(annotation_df, aes(REs, Ann,colour=logFC))+
  geom_point(aes(size=pval))+
  scale_colour_scico(palette = "vik", limit = limit)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size = 5),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "gray60"),
        panel.background = element_rect(fill = "white"))+
  xlab(label = "")+
  ylab(label = "")+
  labs(size = "-log(pval)")


consensus_tfs <- F1A
save(consensus_tfs,tfs,file = "expression_cluster/TFs_cl1_consensus.RData")
tfs
ggsave(filename = "expression_cluster/tfs.tiff",
       width = 10,height = 10)
dev.off()


####### Cluster2 #####
list.of.dfs <- list(diffActivity.adult.cl2[ !(row.names(diffActivity.adult.cl2) %in% tfs_consensus), ],
                    diffActivity.pediatric.cl2[ !(row.names(diffActivity.pediatric.cl2) %in% tfs_consensus), ])
res <- consensus.dataframe.res(list.of.dfs = list.of.dfs,method = "bonferroni",pval=0.05,quantile.value = 0)
res.meta.tfs.consensus <- res[order(res$fc),]
tfs.consensus <- rownames(res.meta.tfs.consensus)
tfs <- tfs.consensus

write.table(file = "expression_cluster/TFs.txt",
            tfs,row.names = F,col.names = F,quote = F)


outPediatric <- pheatmap(TFActivitiesPediatric[c(tfs.consensus),c(pediatric.hc,pediatric.cl2,pediatric.cl2)],
                         scale = "row",breaks = seq(-2,2,length.out = 100),
                         cluster_rows = F,cluster_cols = F,show_colnames = F,
                         show_rownames = F,legend = F,annotation_col = anncolPediatric[c(pediatric.hc,pediatric.cl2,pediatric.cl2),"cluster",drop=F],
                         annotation_colors = list("cluster"=c("Healthy"="#A5FF66",
                                                              "cluster1"="#FF7C4C",
                                                              "cluster2"="#293ECC")),
                         annotation_names_col = F,annotation_legend = F,
                         border_color = NA,fontsize = 5)
outAdult <- pheatmap(TFActivitiesAdult[c(tfs.consensus),c(adult.hc,adult.cl2,adult.cl2)],
                     scale = "row",breaks = seq(-2,2,length.out = 100),
                     cluster_rows = F,cluster_cols = F,show_colnames = F,
                     annotation_col = anncolAdult[c(adult.hc,adult.cl2,adult.cl2),"cluster",drop=F],
                     annotation_colors = list("cluster"=c("Healthy"="#A5FF66",
                                                          "cluster1"="#FF7C4C",
                                                          "cluster2"="#293ECC")),
                     annotation_names_col = F,
                     border_color = NA,fontsize = 5)
items=c("outPediatric","outAdult")
plot_list=list()
for (a in items){
  x=get(a)
  plot_list[[a]] = x[[4]]     ##to save each plot into a list. note the [[4]]
}
F1A<-grid.arrange(arrangeGrob(grobs= plot_list,ncol=2))
ggsave("expression_cluster/Figure_2A.pdf",F1A,width = 6,height = 6)
dev.off()
ggsave("expression_cluster/Figure_2A.png",F1A,width = 8,height = 5)


df_clue <- clue_table(tfs.consensus)


examplePathways <- list()
for (i in 1:length(tfs)){
  tf <- tfs[i]
  examplePathways <- list.append(examplePathways,
                                 tf=names(dorothea_regulon[[tf]]$tfmode))
  names(examplePathways)[i] <- tf
}

fgseaRes.cl2.pediatric <- run_fgsea(expressionPediatric,pediatric.cl2,pediatric.hc,examplePathways)
diffexpression_pediatric.cl2 <- fgseaRes.cl2.pediatric$diffExpression
fgseaRes.cl2.pediatric <- fgseaRes.cl2.pediatric$fgseaRes

fgseaRes.cl2.adult <- run_fgsea(expressionAdult,adult.cl2,adult.hc,examplePathways)
diffexpression_adult.cl2 <- fgseaRes.cl2.adult$diffExpression
fgseaRes.cl2.adult <- fgseaRes.cl2.adult$fgseaRes


res <- lapply(tfs[tfs %in% fgseaRes.cl2.pediatric$pathway],extract_fgsea,fgseaRes=fgseaRes.cl2.pediatric)
res_pediatric.cl2 <- do.call("rbind", res)
rownames(res_pediatric.cl2) <- paste0(res_pediatric.cl2$TF,"_",res_pediatric.cl2$Target)
res <- lapply(tfs[tfs %in% fgseaRes.cl2.adult$pathway],extract_fgsea,fgseaRes=fgseaRes.cl2.adult)
res_adult.cl2 <- do.call("rbind", res)
rownames(res_adult.cl2) <- paste0(res_adult.cl2$TF,"_",res_adult.cl2$Target)

keep_int <- c()
for (int in unique(c(rownames(res_pediatric.cl2),
                     rownames(res_adult.cl2)))){
  a <- res_pediatric.cl2[int,]
  c <- res_adult.cl2[int,]
  total <- na.omit(rbind(a,c))
  if (nrow(total)>=2){
    keep_int <- c(keep_int,int)
  }
}

res <- dorothea_table[rownames(dorothea_table) %in% keep_int,]

diffRes <- consensus.dataframe.res(list.of.dfs = list(diffexpression_adult.cl2,diffexpression_pediatric.cl2),
                                   method = "fdr",pval = 1.1,quantile.value = 0)

diffRes <- diffRes[order(diffRes$fc),]
res <- res[res$Target %in% rownames(diffRes),]
tfs.toplot <- tfs.consensus[tfs.consensus %in% res$TF]

nodes <- unique(c(res$TF,res$Target))
nodes <- data.frame(id=nodes,
                    type=ifelse(nodes %in% tfs,'TF','Gene'),
                    diffExpression=rep(NA,length(nodes)),
                    stringsAsFactors = F,row.names = nodes)

for (tf in tfs.toplot){
  nodes[tf,"diffExpression"] <- sign(res.meta.tfs.consensus[rownames(res.meta.tfs.consensus) == tf,"fc"])
}

for (gene in unique(res$Target)){
  if (!(gene %in% tfs.toplot)){
    direction <- sign(diffRes[gene,"fc"])
    nodes[gene,"diffExpression"] <- as.numeric(direction)
  }
}


res_network <- lapply(tfs.toplot,reduce_dimension,nodes,res)
res_network <- do.call("rbind", res_network)
nodes <- nodes[nodes$id %in% c(res_network$TF,res_network$Target),]

nodes$group <- ifelse(nodes[,"type"]=="TF",
                      ifelse(nodes[,"diffExpression"]==1,"tf_up","tf_down"),
                      ifelse(nodes[,"diffExpression"]==1,"gene_up","gene_down"))




res_network$source <- res_network$TF
res_network$target <- res_network$Target
res_network$interaction <- ifelse(res_network$TFMode==1,
                                  "activation","inhibition")

length(unique(res_network$Target))
length(unique(res_network$TF))

res_network <- res_network[res_network$Target !="T",]

tfs.toplot <- unique(res_network$TF)

calculate_mean_fold_change <- function(split.genes,diffRes){
  split.genes <- unlist(split.genes)
  fc <- mean(diffRes[split.genes,"fc"])
  return(fc)
}

wikipathways_table <- read.delim("GO_BP_table.tsv",stringsAsFactors = F)
wikipathways_table <- wikipathways_table[grep("GC-9606",wikipathways_table$id),]
genes_table <- read.delim("gene_table.tsv",stringsAsFactors = F)
wikipathways_table <- merge(genes_table,wikipathways_table,on="id")[,c(2,5)]
annotation_table <- read.delim("annotation_info_table.tsv",stringsAsFactors = F)
wikipathways_table <- merge(wikipathways_table,annotation_table,on="annotation_id")

res <- mclapply(tfs.toplot,
                annotate_function,
                res_network=res_network,
                wikipathways_table = wikipathways_table,
                diffRes = diffRes,
                mc.cores = 6)

names(res) <- tfs.toplot
res <- res[which(!sapply(res, is.null))]

terms <- c()
for (i in 1:length(res)){
  subRes <- res[i]
  terms <- c(terms,subRes[[1]]$terms)
}

terms <- unique(terms)

annTable <- data.frame(matrix(0,nrow = length(res),ncol = length(terms)))
annGenes <- data.frame(matrix(0,nrow = length(res),ncol = length(terms)))
rownames(annTable) <- rownames(annGenes) <- names(res)
colnames(annTable) <- colnames(annGenes) <- terms


for (i in 1:length(res)){
  subRes <- res[i]
  tf <- names(subRes)
  pvalues <- subRes[[1]]$pvalues
  fcs <- subRes[[1]]$fcs
  terms <- subRes[[1]]$terms
  
  annTable[tf,terms] <- pvalues
  annGenes[tf,terms] <- fcs
}

annotation_table <- as.data.frame(t(annTable[rowSums(annTable != 0) > 0,]))
annotation_genes <- as.data.frame(t(annGenes[colnames(annotation_table),rownames(annotation_table)]))

annotation_df <- data.frame(REs=0,
                            Ann=0,
                            logFC=0,
                            pval=0)

annotation_df <- annotation_df[-1,]


for (term in rownames(annotation_table)){
  for (tf in colnames(annotation_table)){
    pvalues <- annotation_table[term,tf]
    pvalues <- pvalues[pvalues != 0]
    
    fcs <- annotation_genes[term,tf]
    if (length(pvalues)>0){
      n_row <-  data.frame(REs=tf,
                           Ann=term,
                           logFC=fcs,
                           pval=-log(pvalues))
      annotation_df <- rbind(annotation_df,n_row)
    }
  }
}


annotation_df <- annotation_df[annotation_df$pval > -log(0.05) & abs(annotation_df$logFC)>0.5,]


out <- pheatmap(annotation_genes[unique(annotation_df$Ann),],clustering_method = "ward.D2")
term_order <- rownames(annotation_genes[unique(annotation_df$Ann),][out$tree_row[["order"]],])
gene_order <- colnames(annotation_genes[unique(annotation_df$Ann),][,out$tree_col[["order"]]])
annotation_table <- annotation_table[term_order,gene_order]

annotation_df$REs <- factor(annotation_df$REs,levels = colnames(annotation_table))
limit <- c(-max(abs(annotation_df$logFC)),max(abs(annotation_df$logFC)))


tfs <- ggplot(annotation_df, aes(REs, Ann,colour=logFC))+
  geom_point(aes(size=pval))+
  scale_colour_scico(palette = "vik", limit = limit)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size = 5),
        panel.grid.major = element_line(size = 0.15, linetype = 'solid',
                                        colour = "gray60"),
        panel.background = element_rect(fill = "white"))+
  xlab(label = "")+
  ylab(label = "")+
  labs(size = "-log(pval)")


consensus_tfs <- F1A
save(consensus_tfs,tfs,file = "expression_cluster/TFs_cl2_consensus.RData")
tfs
ggsave(filename = "expression_cluster/tfs.tiff",
       width = 10,height = 10)
dev.off()