
data = read.table("gene_exp.diff", head= TRUE, sep='\t' )
head(data)

#sorting data based on q-values(lowest to highest)
data <- data[order(data$q_value),]

#dim(example_data)

#taking top ten values
diff_data <- data[1:10,]


diff_data <- diff_data[ -c(1:2, 4:7, 11,14) ]

#creating histogram of differentially expressed genes
hist2<- hist(data$log2.fold_change., breaks = 15, labels = TRUE,
    
             main = paste("Histogram of log2 fold change"),
xlab = "log2 fold change values")  

#selecting signicant genes
sig_data <- subset(example_data, example_data$significant=="yes")

#histogram of signicant differentailly expressed genes
hist4<- hist(sig_data$log2.fold_change., breaks = 30,labels = TRUE,
             main = paste("Histogram of log2 fold change"),
             xlab = "log2 fold change values for significant genes")  


#choosing up regulated genes based on log2 fold change values
up_data <- subset(sig_data$gene ,sig_data$log2.fold_change.>1)

#choosing down regulated genes based on log2 fold change values
down_data <- subset(sig_data$gene,sig_data$log2.fold_change.<1)

length(up_data)
length(down_data)

#creating csv files for both
write.csv(up_data,"up_genes_project2.csv")
write.csv(down_data,"down_genes_project2.csv")
