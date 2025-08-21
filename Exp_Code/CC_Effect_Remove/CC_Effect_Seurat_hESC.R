library(Seurat)
library(anndata)
library(reticulate)

# 读取数据
res <- anndata::read_h5ad('/home/wenbaole/Cycle/regression/ref_dataset_Leng.h5ad')
temp_matrix <- t(res$X)

# 创建Seurat对象
seurat_obj <- CreateSeuratObject(counts = temp_matrix)

# 数据预处理和归一化
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst")
seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))

# 细胞周期评分
# 加载细胞周期基因集（Seurat内置）
s.genes <- cc.genes$s.genes  # S期基因
g2m.genes <- cc.genes$g2m.genes  # G2/M期基因

# 计算细胞周期评分
seurat_obj <- CellCycleScoring(seurat_obj,
                               s.features = s.genes,
                               g2m.features = g2m.genes,
                               set.ident = TRUE)

# 查看细胞周期分布
print(table(seurat_obj@meta.data$Phase))

# 去除细胞周期效应
seurat_obj <- ScaleData(seurat_obj,
                        vars.to.regress = c("S.Score", "G2M.Score"),
                        features = rownames(seurat_obj))

# 获取去除细胞周期效应后的矩阵
corrected_matrix <- GetAssayData(seurat_obj, assay = "RNA", slot = "scale.data")

# 查看结果
print(dim(corrected_matrix))
print("原始矩阵维度:")
print(dim(temp_matrix))
print("去除细胞周期效应后矩阵维度:")
print(dim(corrected_matrix))


# 可选：保存结果
# write.csv(corrected_matrix, "cell_cycle_corrected_matrix.csv")

corrected_seurat <- CreateSeuratObject(counts = corrected_matrix)
# 设置corrected矩阵到scale.data槽位
corrected_seurat <- SetAssayData(corrected_seurat, 
                                 slot = "scale.data", 
                                 new.data = corrected_matrix)

corrected_seurat <- RunPCA(corrected_seurat, 
                           features = rownames(corrected_matrix),
                           npcs = 2)

# 获取PC2坐标
pc2_coordinates <- Embeddings(corrected_seurat, reduction = "pca")[, "PC_2"]
pc1_coordinates <- Embeddings(corrected_seurat, reduction = "pca")[, "PC_1"]
# 从原始adata中获取stage标签
stage_labels <- res$obs$stage

# 确保标签与细胞顺序一致
if (length(stage_labels) != length(pc2_coordinates)) {
  stop("标签数量与细胞数量不匹配")
}

library(cluster)

pca_2d_df <- data.frame(
  PC1 = pc1_coordinates,
  PC2 = pc2_coordinates,
  row.names = names(pc1_coordinates)
)

# 计算二维空间中的距离矩阵
pca_2d_dist <- dist(pca_2d_df)

# 计算轮廓系数（基于PC1和PC2二维空间）
silhouette_score_2d <- silhouette(as.numeric(as.factor(stage_labels)), 
                                  dist = as.matrix(pca_2d_dist))

# 输出结果
cat("=== 基于PC1和PC2二维空间的轮廓系数 ===\n")
cat("平均轮廓系数:", mean(silhouette_score_2d[, "sil_width"]), "\n")
print(summary(silhouette_score_2d))

result_df <- data.frame(
  cell_id = rownames(pca_2d_df),
  PC1 = pca_2d_df$PC1,
  PC2 = pca_2d_df$PC2,
  stage = stage_labels,
  silhouette_width = silhouette_score_2d[, "sil_width"]
)

library(ggplot2)

# PC1 vs PC2散点图着色按stage，点大小表示轮廓系数
p2 <- ggplot(result_df, aes(x = PC1, y = PC2, color = as.factor(stage), size = silhouette_width)) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(range = c(1, 4)) +
  scale_color_discrete(name = "Stage") +  # 使用离散色彩标度
  labs(title = "PC1 vs PC2 by Stage (点大小表示轮廓系数)",
       x = "PC1", y = "PC2") +
  theme_minimal() +
  guides(size = guide_legend(title = "Silhouette Width"),
         color = guide_legend(title = "Stage"))

p2

save(seurat_obj, corrected_seurat, file = "/home/wenbaole/Cycle/res//CC_hESC_seurat_objects.rda")

