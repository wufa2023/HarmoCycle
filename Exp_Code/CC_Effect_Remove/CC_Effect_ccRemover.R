# ===================================================================
# 步骤 1: 安装和加载 ccRemover 包
# ===================================================================

res <- anndata::read_h5ad('/root/Cycle/HarmoCycle_v1/Fig5/ref1-6.h5ad')
temp_matrix <- t(res$X)

# 检查 ccRemover 包是否已经安装，如果没有，则自动从 CRAN 安装
if (!requireNamespace("ccRemover", quietly = TRUE)) {
  install.packages("ccRemover")
}

# 加载包
library(ccRemover)

cat("--- ccRemover 包已成功加载 ---\n\n")


# ===================================================================
# 步骤 2: 加载示例数据并进行预处理
# ===================================================================

# 为保证结果可复现，设置随机种子
set.seed(10)

# # 加载包内自带的示例数据 t.cell_data
# data(t.cell_data)

# t_cell_data_subset <- t.cell_data

# # 1. 转置为 基因 x 细胞
data_matrix_genes_by_cells <- temp_matrix
cat("\n转置后的数据维度 (基因 x 细胞):\n")
print(dim(data_matrix_genes_by_cells))

# 2. 对每一行（基因）进行中心化
#    我们使用 t(scale(t(...))) 的技巧，这样可以方便地在行上操作
#    注意：scale 函数默认在列上操作，所以需要两次转置
log_data <- log1p(data_matrix_genes_by_cells)
scaled_data_transposed <- scale(t(log_data), center = TRUE, scale = FALSE)
data_centered <- t(scaled_data_transposed)

# 检查一下某一行（基因）的均值是否接近于0，以验证中心化是否成功
cat("\n检查第一个基因中心化后的均值 (应接近0):\n")
print(mean(data_centered[1, ]))
cat("\n预处理完成！\n\n")


gene_means <- attr(scaled_data_transposed, "scaled:center")

cat("--- 已保存每个基因的平均表达量 ---\n")
cat("前5个基因的平均值:\n")
print(head(gene_means, 5))
cat("\n")

gene_names <- rownames(data_centered)
cat("--- 正在使用 gene_indexer 识别细胞周期基因 ---\n")
cell_cycle_gene_indices <- gene_indexer(
  gene_names = gene_names,
  species = "human",
  name_type = "symbol"
)
if_cc <- rep(FALSE, nrow(data_centered))
if_cc[cell_cycle_gene_indices] <- TRUE
cat(paste("在", length(gene_names), "个基因中，识别出", sum(if_cc), "个细胞周期基因。\n\n"))

dat_input <- list(x = data_centered, if_cc = if_cc)

cat("--- 正在运行 ccRemover，请稍候... ---\n")
xhat_centered <- ccRemover(dat_input, bar = TRUE)
cat("--- ccRemover 运行完成！ ---\n\n")

cat("校正后但仍为中心化的矩阵 `xhat_centered` 的维度 (基因 x 细胞):\n")
print(dim(xhat_centered))

cat("\n中心化校正数据预览 (注意值围绕0分布，有负值):\n")
print(xhat_centered[1:5, 1:5])
cat("\n")

xhat_final <- sweep(xhat_centered, MARGIN = 1, STATS = gene_means, FUN = "+")

cat("--- 已将基因平均值加回 ---\n")
cat("最终校正后矩阵 `xhat_final` 的维度 (应与之前相同):\n")
print(dim(xhat_final))

# 查看最终校正后的数据，现在的值恢复了其原始的量级，且大部分为非负数
cat("\n最终校正后数据预览 (值已恢复正常量级):\n")
print(xhat_final[1:5, 1:5])
cat("\n--- 所有步骤完成！`xhat_final` 已可用于下游分析。---\n")

write.csv(xhat_centered, '/root/Cycle/HarmoCycle_v1/Fig5/after_log_6datasets_xhat.csv')
write.csv(xhat_final, '/root/Cycle/HarmoCycle_v1/Fig5/after_log_6datasets_final.csv')
