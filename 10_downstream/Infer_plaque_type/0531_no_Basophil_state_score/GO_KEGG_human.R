library(openxlsx)#读取.xlsx文件
library(ggplot2)#柱状图和点状图
library(stringr)#基因ID转换
library(ggnewscale)
library(circlize)#绘制富集分析圈图
library(ComplexHeatmap)#绘制图例

library(enrichplot)#GO,KEGG,GSEA
library(clusterProfiler)#GO,KEGG,GSEA
library(GOplot)#弦图，弦表图，系统聚类图
library(DOSE)
library(topGO)#绘制通路网络图

###############################################unstable########################################################3
#载入差异表达数据，只需基因ID(GO,KEGG,GSEA需要)和Log2FoldChange(GSEA需要)即可
info <- read.xlsx( "/home/lixiangyu/zr/Annotate/ANNOTATE_new/10_downstream/Infer_plaque_type/0531_no_Basophil_state_score/output_DESeq2/human/unstable_up.xlsx", rowNames = F,colNames = T)
#指定富集分析的物种库
GO_database <- 'org.Hs.eg.db' #GO分析指定物种，物种缩写索引表详见http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
KEGG_database <- 'hsa' #KEGG分析指定物种，物种缩写索引表详见http://www.genome.jp/kegg/catalog/org_list.html
#gene ID转换
gene <- bitr(info$gene_symbol,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
head(gene)

##GO分析
GO<-enrichGO( gene$ENTREZID,#GO富集分析
              OrgDb = GO_database,
              keyType = "ENTREZID",#设定读取的gene ID类型
              ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
              pvalueCutoff = 0.05,#设定p值阈值
              qvalueCutoff = 0.05,#设定q值阈值
              readable = T
)
GO
##KEGG分析
KEGG<-enrichKEGG(gene$ENTREZID,#KEGG富集分析
                 organism = KEGG_database,
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
KEGG
###GO/KEGG富集柱状图+点状图
outdir <- "/home/lixiangyu/zr/Annotate/ANNOTATE_new/10_downstream/Infer_plaque_type/0531_no_Basophil_state_score/output_GO_KEGG/output_human/unstable"
p1 <- barplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")#柱状图
p2 <- barplot(KEGG,showCategory = 40,title = 'KEGG Pathway')
ggsave(filename = file.path(outdir, "GO_barplot.png"),plot = p1,width = 10,height = 8,dpi = 300)
ggsave(filename = file.path(outdir, "KEGG_barplot.png"),plot = p2,width = 10,height = 8,dpi = 300)
p3 <- dotplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")#点状图
p4 <- dotplot(KEGG)
ggsave(filename = file.path(outdir, "GO_dotplot.png"),plot = p3,width = 10,height = 8,dpi = 300)
ggsave(filename = file.path(outdir, "KEGG_dotplot.png"),plot = p4,width = 10,height = 8,dpi = 300)

###############################################stable#########################################################
#载入差异表达数据，只需基因ID(GO,KEGG,GSEA需要)和Log2FoldChange(GSEA需要)即可
info <- read.xlsx( "/home/lixiangyu/zr/Annotate/ANNOTATE_new/10_downstream/Infer_plaque_type/0531_no_Basophil_state_score/output_DESeq2/human/stable_up.xlsx", rowNames = F,colNames = T)
#指定富集分析的物种库
GO_database <- 'org.Hs.eg.db' #GO分析指定物种，物种缩写索引表详见http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
KEGG_database <- 'hsa' #KEGG分析指定物种，物种缩写索引表详见http://www.genome.jp/kegg/catalog/org_list.html
#gene ID转换
gene <- bitr(info$gene_symbol,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
head(gene)

##GO分析
GO<-enrichGO( gene$ENTREZID,#GO富集分析
              OrgDb = GO_database,
              keyType = "ENTREZID",#设定读取的gene ID类型
              ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
              pvalueCutoff = 0.05,#设定p值阈值
              qvalueCutoff = 0.05,#设定q值阈值
              readable = T
)
GO
##KEGG分析
KEGG<-enrichKEGG(gene$ENTREZID,#KEGG富集分析
                 organism = KEGG_database,
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
KEGG
###GO/KEGG富集柱状图+点状图
outdir <- "/home/lixiangyu/zr/Annotate/ANNOTATE_new/10_downstream/Infer_plaque_type/0531_no_Basophil_state_score/output_GO_KEGG/output_human/stable"
p1 <- barplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")#柱状图
p2 <- barplot(KEGG,showCategory = 40,title = 'KEGG Pathway')
ggsave(filename = file.path(outdir, "GO_barplot.png"),plot = p1,width = 10,height = 8,dpi = 300)
ggsave(filename = file.path(outdir, "KEGG_barplot.png"),plot = p2,width = 10,height = 8,dpi = 300)
p3 <- dotplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")#点状图
p4 <- dotplot(KEGG)
ggsave(filename = file.path(outdir, "GO_dotplot.png"),plot = p3,width = 10,height = 8,dpi = 300)
ggsave(filename = file.path(outdir, "KEGG_dotplot.png"),plot = p4,width = 10,height = 8,dpi = 300)
