library(Seurat)
library(Scissor)
library(Matrix)
library(preprocessCore)
#---------教程源代码--------


location <- "https://xialab.s3-us-west-2.amazonaws.com/Duanchen/Scissor_data/"
load("E:\\PostGraduate\\phenotype-spatial\\Scissors\\data from article\\scRNA-seq.RData")
sc_dataset <- Seurat_preprocessing(sc_dataset, verbose = F)

DimPlot(sc_dataset, reduction = 'umap', label = T, label.size = 10)
load("E:\\PostGraduate\\phenotype-spatial\\Scissors\\data from article\\TCGA_LUAD_exp1.RData")
load("E:\\PostGraduate\\phenotype-spatial\\Scissors\\data from article\\TCGA_LUAD_survival.RData")

phenotype <- bulk_survival[,2:3]
colnames(phenotype) <- c("time", "status")
head(phenotype)

infos1 <- Scissor(bulk_dataset, sc_dataset, phenotype, alpha = 0.05, 
                  family = "cox", Save_file = 'Scissor_LUAD_survival.RData')


Scissor_select <- rep(0, ncol(sc_dataset))
names(Scissor_select) <- colnames(sc_dataset)
Scissor_select[infos1$Scissor_pos] <- 1
Scissor_select[infos1$Scissor_neg] <- 2
sc_dataset <- AddMetaData(sc_dataset, metadata = Scissor_select, col.name = "scissor")
DimPlot(sc_dataset, reduction = 'umap', group.by = 'scissor', cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1))



















#-------------做修改后的代码----------
load("E:\\PostGraduate\\phenotype-spatial\\Scissors\\data from article\\scRNA-seq.RData")
load("E:\\PostGraduate\\phenotype-spatial\\Scissors\\data from article\\TCGA_LUAD_exp1.RData")
load("E:\\PostGraduate\\phenotype-spatial\\Scissors\\data from article\\TCGA_LUAD_survival.RData")
sc_count <- CreateSeuratObject(counts = sc_dataset, project = "Scissor_Single_Cell", 
                             min.cells = 400, min.features = 0)

sc_count[["percent.mt"]] <- PercentageFeatureSet(sc_count, pattern = "^MT[-]")
VlnPlot(sc_count,features = c("nCount_Spatial","nFeature_Spatial","percent.mt"))


sc_count <- NormalizeData(object = sc_count, normalization.method = "LogNormalize", 
                      scale.factor = 10000, verbose = TRUE)
# 先进行PCA降维，再选择前30个维度进行聚类和umap降维
sc_count <- FindVariableFeatures(object = sc_count, selection.method = "vst", 
                             verbose = TRUE)
sc_count <- ScaleData(object = sc_count, verbose = TRUE)
sc_count <- RunPCA(object = sc_count, features = VariableFeatures(sc_count), 
               verbose = TRUE)
sc_count <- FindNeighbors(object = sc_count, dims = 1:10, 
                      verbose = TRUE)
sc_count <- FindClusters(object = sc_count, resolution = 0.6, 
                     verbose = TRUE)
sc_count <- RunTSNE(object = sc_count, dims = 1:10)
sc_count <- RunUMAP(object = sc_count, dims = 1:10, verbose = TRUE)

phenotype <- bulk_survival[,2:3]
colnames(phenotype) <- c("time", "status")
head(phenotype)























#Scissor 源代码
common <- intersect(rownames(bulk_dataset), rownames(sc_count))
if (length(common) == 0) {
  stop("There is no common genes between the given single-cell and bulk samples.")
}
if (class(sc_count) == "Seurat") {
  sc_exprs <- as.matrix(sc_count@assays$RNA@data)
  network <- as.matrix(sc_count@graphs$RNA_snn)
}

family<- "cox"
diag(network) <- 0
network[which(network != 0)] <- 1
dataset0 <- cbind(bulk_dataset[common, ], sc_exprs[common, 
])
dataset1 <- normalize.quantiles(dataset0)
rownames(dataset1) <- rownames(dataset0)
colnames(dataset1) <- colnames(dataset0)
Expression_bulk <- dataset1[, 1:ncol(bulk_dataset)]
Expression_cell <- dataset1[, (ncol(bulk_dataset) + 
                                 1):ncol(dataset1)]
X <- cor(Expression_bulk, Expression_cell)
quality_check <- quantile(X)
print("|**************************************************|")
print("Performing quality-check for the correlations")
print("The five-number summary of correlations:")
print(quality_check)
print("|**************************************************|")
if (quality_check[3] < 0.01) {
  warning("The median correlation between the single-cell and bulk samples is relatively low.")
}
if (family == "binomial") {
  Y <- as.numeric(phenotype)
  z <- table(Y)
  if (length(z) != length(tag)) {
    stop("The length differs between tags and phenotypes. Please check Scissor inputs and selected regression type.")
  }
  else {
    print(sprintf("Current phenotype contains %d %s and %d %s samples.", 
                  z[1], tag[1], z[2], tag[2]))
    print("Perform logistic regression on the given phenotypes:")
  }
}
if (family == "gaussian") {
  Y <- as.numeric(phenotype)
  z <- table(Y)
  if (length(z) != length(tag)) {
    stop("The length differs between tags and phenotypes. Please check Scissor inputs and selected regression type.")
  }
  else {
    tmp <- paste(z, tag)
    print(paste0("Current phenotype contains ", 
                 paste(tmp[1:(length(z) - 1)], collapse = ", "), 
                 ", and ", tmp[length(z)], " samples."))
    print("Perform linear regression on the given phenotypes:")
  }
}
if (family == "cox") {
  Y <- as.matrix(phenotype)
  if (ncol(Y) != 2) {
    stop("The size of survival data is wrong. Please check Scissor inputs and selected regression type.")
  }
  else {
    print("Perform cox regression on the given clinical outcomes:")
  }
}
# save(X, Y, network, Expression_bulk, Expression_cell, 
#      file = Save_file)

alpha <- c()
if (is.null(alpha)) {
  alpha <- c(0.005, 0.05, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 
             0.6, 0.7, 0.8, 0.9)
}
cutoff <- 0.2

alpha <- c(0.05)

nrow(X)
length(Y[,"status"])
for (i in 1:length(alpha)) {
  set.seed(123)
  fit0 <- APML1(X, Y, family = family, penalty = "Net", 
                alpha = alpha[i], Omega = network, nlambda = 100, 
                nfolds = min(10, nrow(X)))
  fit1 <- APML1(X, Y, family = family, penalty = "Net", 
                alpha = alpha[i], Omega = network, lambda = fit0$lambda.min)
  if (family == "binomial") {
    Coefs <- as.numeric(fit1$Beta[2:(ncol(X) + 1)])
  }
  else {
    Coefs <- as.numeric(fit1$Beta)
  }
  Cell1 <- colnames(X)[which(Coefs > 0)]
  Cell2 <- colnames(X)[which(Coefs < 0)]
  percentage <- (length(Cell1) + length(Cell2))/ncol(X)
  
  print(sprintf("alpha = %s", alpha[i]))
  print(sprintf("Scissor identified %d Scissor+ cells and %d Scissor- cells.", 
                length(Cell1), length(Cell2)))
  print(sprintf("The percentage of selected cell is: %s%%", 
                formatC(percentage * 100, format = "f", digits = 3)))
  if (percentage < cutoff) {
    break
  }
  cat("\n")
}

infos1 <- list(para = list(alpha = alpha[i], lambda = fit0$lambda.min, 
                           family = family), Coefs = Coefs, Scissor_pos = Cell1, 
               Scissor_neg = Cell2)

#查看结果
names(infos1)
length(infos1$Scissor_pos)
infos1$Scissor_pos[1:4]
length(infos1$Scissor_neg)
infos1$Scissor_neg
#将注释的和表型相关的注释信息加在Seurat object上面
Scissor_select <- rep(0, ncol(sc_count))
names(Scissor_select) <- colnames(sc_count)
Scissor_select[infos1$Scissor_pos] <- 1
Scissor_select[infos1$Scissor_neg] <- 2
sc_count <- AddMetaData(sc_count, metadata = Scissor_select, col.name = "scissor")
DimPlot(sc_count, reduction = 'umap', group.by = 'scissor', cols = c('grey','indianred1','royalblue'), pt.size = 2, order = c(2,1))













