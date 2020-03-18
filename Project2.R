library(reshape2)
library(ggplot2)
library(dplyr)
library(plyr)

fpkm_matrix = read.table("/project/bf528/project_2/data/fpkm_matrix.csv",header=TRUE)
fpkm_matrix

gene_exp = read.table("/project/bf528/project_2/data/P4_vs_P7_cuffdiff_out/gene_exp.diff", header=TRUE)
gene_exp

gene_exp1 = read.table("/projectnb/bf528/users/group6/project_2/analysis/cuffdiff_out/gene_exp.diff", header=TRUE)
gene_exp1
       
fpkm = read.table("/projectnb/bf528/users/group6/project_2/analysis/cuffdiff_out/genes.fpkm_tracking", header=TRUE)
fpkm


df1 = data.frame(fpkm$gene_short_name, fpkm$P0_FPKM, fpkm$Ad_FPKM)
target = c("Mpc1", "Prdx3", "Acat1", "Echs1", "Slc25a11", "Phyh")
index=df1$fpkm.gene_short_name %in% target

df1_mit=df1[index,]
df1_mit=melt(df1_mit)

df1_mit$variable = as.character(df1_mit$variable)
df1_mit$variable[df1_mit$variable == "fpkm.P0_FPKM"] <-"P0"
df1_mit$variable[df1_mit$variable == "fpkm.Ad_FPKM"] <-"Ad"


#Getting 
fpkm_matrix$P4 = (fpkm_matrix$P4_1_FPKM + fpkm_matrix$P4_2_FPKM)/2
fpkm_matrix$P7 = (fpkm_matrix$P7_1_FPKM + fpkm_matrix$P7_2_FPKM)/2

drops <- c("Ad_1_FPKM","Ad_2_FPKM", "P0_2_FPKM", "P4_1_FPKM", "P4_2_FPKM", "P7_1_FPKM", "P7_2_FPKM", "fpkm_mean_P4", "fpkm_mean_P7")
updated_fpkm = fpkm_matrix[ , !(names(fpkm_matrix) %in% drops)]


#Mitochondria
get1 = c("ENSMUSG00000024997", "ENSMUSG00000032047", "ENSMUSG00000025465", "ENSMUSG00000014606", "ENSMUSG00000026664")  
index1=updated_fpkm$tracking_id %in% get1
p47_mit=updated_fpkm[index1,]
p47_mit=melt(p47_mit)   
p47_mit$tracking_id = ifelse(p47_mit$tracking_id=="ENSMUSG00000024997", "Prdx3", ifelse(p47_mit$tracking_id=="ENSMUSG00000032047", "Acat1", ifelse(p47_mit$tracking_id=="ENSMUSG00000025465", "Echs1", ifelse(p47_mit$tracking_id=="ENSMUSG00000014606", "Slc25a11", ifelse(p47_mit$tracking_id=="ENSMUSG00000026664", "Phyh", NA))) ))

colnames(p47_mit)[1] <- "fpkm.gene_short_name"

x=rbind(p47_mit, df1_mit)
x
x$new=factor(x$variable, labels=c("P0", "P4", "P7", "Ad"), levels=c("P0", "P4", "P7", "Ad"), ordered = TRUE)


ggplot(x, aes(x=x$new, y=x$value, group=fpkm.gene_short_name)) +
  geom_line(aes(color=fpkm.gene_short_name)) +
  geom_point(aes(color=fpkm.gene_short_name)) +
  labs(title="Mitochondria",x="", y = "FPKM") +
  theme(legend.title = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5))

#Sarcomere
df2 = data.frame(fpkm$gene_short_name, fpkm$P0_FPKM, fpkm$Ad_FPKM)
target = c("Pdlim5", "Pygm", "Myoz2", "Des", "Csrp3", "Tcap","Cryab")
index=df2$fpkm.gene_short_name %in% target

df_sac=df2[index,]
df_sac=melt(df_sac)

df_sac$variable = as.character(df_sac$variable)
df_sac$variable[df_sac$variable == "fpkm.P0_FPKM"] <-"P0"
df_sac$variable[df_sac$variable == "fpkm.Ad_FPKM"] <-"Ad"

get1 = c("ENSMUSG00000028273", "ENSMUSG00000032648", "ENSMUSG00000028116", "ENSMUSG00000026208", "ENSMUSG00000030470", "ENSMUSG00000007877", "ENSMUSG00000032060")  
index1=updated_fpkm$tracking_id %in% get1
p47_sac=updated_fpkm[index1,]
p47_sac=melt(p47_sac)   
p47_sac$tracking_id = ifelse(p47_sac$tracking_id=="ENSMUSG00000028273", "Pdlim5", ifelse(p47_sac$tracking_id=="ENSMUSG00000032648", "Pygm", ifelse(p47_sac$tracking_id=="ENSMUSG00000028116", "Myoz2", ifelse(p47_sac$tracking_id=="ENSMUSG00000026208", "Des", ifelse(p47_sac$tracking_id=="ENSMUSG00000030470", "Csrp3", ifelse(p47_sac$tracking_id=="ENSMUSG00000007877", "Tcap",  ifelse(p47_sac$tracking_id=="ENSMUSG00000032060", "Cryab", NA)))))))

colnames(p47_sac)[1] <- "fpkm.gene_short_name"

x1=rbind(p47_sac, df_sac)
x1
x1$new=factor(x1$variable, labels=c("P0", "P4", "P7", "Ad"), levels=c("P0", "P4", "P7", "Ad"), ordered = TRUE)


ggplot(x1, aes(x=x1$new, y=x1$value, group=fpkm.gene_short_name)) +
  geom_line(aes(color=fpkm.gene_short_name)) +
  geom_point(aes(color=fpkm.gene_short_name)) +
  labs(title="Sarcomere",x="", y = "FPKM") +
  theme(legend.title = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5))

#Cell cycle
df3 = data.frame(fpkm$gene_short_name, fpkm$P0_FPKM, fpkm$Ad_FPKM)
target = c("Cdc7", "E2f8", "Cdk7", "Cdc26", "Cdc6", "E2f1","Cdc27", "Cdc45", "Rad51", "Aurkb", "Cdc23")
index=df3$fpkm.gene_short_name %in% target

df_cc=df2[index,]
df_cc=melt(df_cc)

df_cc$variable = as.character(df_cc$variable)
df_cc$variable[df_cc$variable == "fpkm.P0_FPKM"] <-"P0"
df_cc$variable[df_cc$variable == "fpkm.Ad_FPKM"] <-"Ad"

get1 = c("ENSMUSG00000029283", "ENSMUSG00000046179", "ENSMUSG00000069089", "ENSMUSG00000066149", "ENSMUSG00000017499", "ENSMUSG00000027490", "ENSMUSG00000020687", "ENSMUSG00000000028", "ENSMUSG00000027323", "ENSMUSG00000020897", "ENSMUSG00000024370")  
index1=updated_fpkm$tracking_id %in% get1
p47_cc=updated_fpkm[index1,]
p47_cc=melt(p47_cc)   
p47_cc$tracking_id = ifelse(p47_cc$tracking_id=="ENSMUSG00000029283", "Cdc7", ifelse(p47_cc$tracking_id=="ENSMUSG00000046179", "E2f8", ifelse(p47_cc$tracking_id=="ENSMUSG00000069089", "Cdk7", ifelse(p47_cc$tracking_id=="ENSMUSG00000066149", "Cdc26", ifelse(p47_cc$tracking_id=="ENSMUSG00000017499", "Cdc6", ifelse(p47_cc$tracking_id=="ENSMUSG00000027490", "E2f1",  ifelse(p47_cc$tracking_id=="ENSMUSG00000020687", "Cdc27", ifelse(p47_cc$tracking_id=="ENSMUSG00000000028", "Cdc45",  ifelse(p47_cc$tracking_id== "ENSMUSG00000027323", "Rad51", ifelse(p47_cc$tracking_id== "ENSMUSG00000020897", "Aurkb", ifelse(p47_cc$tracking_id== "ENSMUSG00000024370", "Cdc23", NA)))))))))))

colnames(p47_cc)[1] <- "fpkm.gene_short_name"

x2=rbind(p47_cc, df_cc)
x2
x2$new=factor(x2$variable, labels=c("P0", "P4", "P7", "Ad"), levels=c("P0", "P4", "P7", "Ad"), ordered = TRUE)

ggplot(x2, aes(x=x2$new, y=x2$value, group=fpkm.gene_short_name)) +
  geom_line(aes(color=fpkm.gene_short_name)) +
  geom_point(aes(color=fpkm.gene_short_name)) +
  labs(title="Cell cycle",x="", y = "FPKM")+
  theme(legend.title = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5))

#Creating heatmap of top 1000 genes
inf = c("Inf","-Inf")
gene_exp1_rank = gene_exp1[!(gene_exp1$log2.fold_change. %in% inf), ]
gene_exp1_rank
gene_exp1_rank=gene_exp1_rank[order(gene_exp1_rank$log2.fold_change.),]

newdata <- mtcars[order(mpg),]

top= head(gene_exp1_rank, 500)
top
  
bot= tail(gene_exp1_rank, 500)
bot
top_genes = rbind(top, bot)
top_genes=subset(top_genes, select=gene)

fpkm_matrix = subset(fpkm_matrix, select=-c(P4, P7))
fpkm = subset(fpkm, select=c(gene_short_name, P0_FPKM))
names(fpkm)[names(fpkm)=="P0_FPKM"]="P0_1_FPKM"

supplement<-read.table("/usr4/bf527/smit2/supplement2.csv",sep = ",", header=TRUE)
supplement2 = supplement[which(supplement$GeneID %in% top_genes$gene),]
supplement2 = subset(supplement2, select=c(Gene, GeneID))
names(supplement2)[names(supplement2)=="GeneID"]="gene_short_name"

merge1 = merge(supplement2, fpkm, by="gene_short_name")
names(merge1)[names(merge1)=="Gene"]="tracking_id"

merge2= merge(merge1, fpkm_matrix, by="tracking_id", na.rm = T)

merge2 = subset(merge2, select=-(tracking_id))

row.names(merge2)=make.names(merge2[,1], TRUE)
merge2=merge2[,-1]
merge2 <- merge2[c(1,4,5,6,7,8,2,3)]
names(merge2)[names(merge2)=="P0_1_FPKM"] <- "P0_1"
names(merge2)[names(merge2)=="P0_2_FPKM"] <- "P0_2"
names(merge2)[names(merge2)=="P1_1_FPKM"] <- "P1_1"
names(merge2)[names(merge2)=="P1_2_FPKM"] <- "P1_2"
names(merge2)[names(merge2)=="P4_1_FPKM"] <- "P4_1"
names(merge2)[names(merge2)=="P4_2_FPKM"] <- "P4_2"
names(merge2)[names(merge2)=="P7_1_FPKM"] <- "P7_1"
names(merge2)[names(merge2)=="P7_2_FPKM"] <- "P7_2"
names(merge2)[names(merge2)=="Ad_1_FPKM"] <- "Ad_1"
names(merge2)[names(merge2)=="Ad_2_FPKM"] <- "Ad_2"

merge_matrix = data.matrix(merge2)

mode(merge_matrix)='numeric'
heatmap(merge_matrix, scale="row", name)


