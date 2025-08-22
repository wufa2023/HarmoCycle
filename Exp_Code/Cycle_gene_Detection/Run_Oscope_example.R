library(Oscope)
library(anndata)


# 读取h5ad文件（假设数据存储在.X中，行为基因，列为细胞）
adata <- read_h5ad("/remote-home/share/data2/wbl/wbl_cellcycle/ground-truth/ref_dataset_Leng.h5ad")  # 替换为实际文件路径
expr <- as.matrix(adata$X)     # 提取表达矩阵（G×S）
expr_matrix <- t(expr)

# 检查数据维度：确保行为基因，列为样本
stopifnot(nrow(expr_matrix) > 0 & ncol(expr_matrix) > 0)
#rownames(expr_matrix) 
#colnames(expr_matrix)


# 使用中位数归一化（MedianNorm，类似DESeq方法）
Sizes <- MedianNorm(expr_matrix)

# 或使用分位数归一化（如75%分位数归一化）
# Sizes <- QuantileNorm(expr_matrix, 0.75)

# 获取归一化后的表达矩阵
DataNorm <- GetNormalizedMat(expr_matrix, Sizes)


# 计算均值和方差，筛选基因
MV <- CalcMV(Data = DataNorm, Sizes = NULL, NormData = TRUE)  # NormData=TRUE表示输入已归一化数据
DataSubset <- DataNorm[MV$GeneToUse, ]  # 提取筛选后的基因子集

#sum_genes <- length(MV$GeneToUse)
("最终保留的基因数:", length(MV$GeneToUse))

write.table(MV$GeneToUse, file = "/root/Cycle/HarmoCycle_v1/Result/Oscope/Leng_GeneToUse.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# %% ---- 开始分析流程并计时 ----

# 加载必要的库 (请确保 oscope 及其依赖已安装)

# 记录开始时间
start_time <- Sys.time()
print(paste("分析开始时间:", start_time))


# --- 1. 数据预处理 ---
# 将表达值缩放到[-1, 1]，并处理离群值
print("步骤 1: 数据标准化...")
DataInput <- NormForSine(DataSubset)


# --- 2. 识别共振荡基因 ---
# 运行配对正弦模型
print("步骤 2: 运行 OscopeSine 识别共振荡基因...")
SineRes <- OscopeSine(DataInput, parallel = FALSE)  # 可并行计算加速
end_time <- Sys.time()
print(paste("分析结束时间:", end_time))

# 计算总运行时间
total_runtime <- end_time - start_time

# 打印总运行时间
print("-----------------------------------------")
print("分析完成！")
print(paste("总运行时间:"))
print(total_runtime)
print("-----------------------------------------")

save(SineRes, file = "/root/Cycle/HarmoCycle_v1/Result/Oscope/Leng_SineRes.rda")


print("步骤 3: 使用 OscopeKM 聚类基因群...")
KMRes <- OscopeKM(SineRes, maxK = 3)

# 查看聚类结果
print("聚类结果预览:")
print(KMRes)

save(KMRes, file = "/root/Cycle/HarmoCycle_v1/Result/Oscope/Leng_KMRes.rda")


unique_genes <- union(KMRes$cluster1, KMRes$cluster2)



# 初始化 unique_genes 为一个空集合
unique_genes <- c()

# 遍历 KMRes 中的所有聚类
for (cluster_name in names(KMRes)) {
    # 将当前聚类的基因合并到 unique_genes 中
    unique_genes <- union(unique_genes, KMRes[[cluster_name]])
}

# 打印 unique_genes 的长度
cat("Unique genes count:", length(unique_genes), "\n")


