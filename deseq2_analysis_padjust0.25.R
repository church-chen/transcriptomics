#!/usr/bin/env Rscript

################################################################################
# DESeq2差异表达分析脚本 - 宽松阈值版本
# 阈值设置: padj < 0.25, |log2FC| > 0.5
# 目标: 找到100+个差异基因
################################################################################

# 加载必需的R包
cat("================================================================================\n")
cat("DESeq2差异分析 - 宽松阈值版本\n")
cat("阈值: padj < 0.25, |log2FC| > 0.5\n")
cat("================================================================================\n\n")

cat("加载R包...\n")
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(ggrepel)
  library(dplyr)
})

################################################################################
# 1. 读取数据
################################################################################
cat("\n1. 读取基因计数数据...\n")
counts_data <- read.table("gene_counts.txt", 
                         header = TRUE, 
                         row.names = 1,
                         skip = 1,
                         check.names = FALSE)

cat("   数据维度:", dim(counts_data), "\n")

# 提取counts矩阵
counts_matrix <- counts_data[, 6:ncol(counts_data)]
colnames(counts_matrix) <- gsub(".bam$", "", colnames(counts_matrix))
colnames(counts_matrix) <- gsub("03_alignment.", "", colnames(counts_matrix))

cat("   样本名称:", paste(colnames(counts_matrix), collapse = ", "), "\n")

# 保存原始counts矩阵（用于生成质控前的PCA图）
counts_matrix_all <- counts_matrix

################################################################################
# 2. 样本质量控制 - 排除PCA离群样本
################################################################################
cat("\n2. 样本质量控制...\n")

# 根据PCA分析结果，排除严重离群的样本
outlier_samples <- c("Mod-53", "ZYD-F-70", "ZYD-F-62")
cat("   排除离群样本:", paste(outlier_samples, collapse = ", "), "\n")
cat("   原因: PCA图显示这些样本严重偏离各自组的中心\n")

# 过滤样本
counts_matrix <- counts_matrix[, !(colnames(counts_matrix) %in% outlier_samples)]
cat("   保留样本数:", ncol(counts_matrix), "\n")
cat("   保留样本:", paste(colnames(counts_matrix), collapse = ", "), "\n")

################################################################################
# 3. 准备样本信息
################################################################################
cat("\n3. 准备样本信息...\n")
sample_names <- colnames(counts_matrix)
group <- ifelse(grepl("^Mod-", sample_names), "Mod", "ZYD")

sample_info <- data.frame(
  sample = sample_names,
  group = factor(group, levels = c("Mod", "ZYD")),
  row.names = sample_names
)

cat("   样本分组:\n")
print(table(sample_info$group))

################################################################################
# 4. 基因过滤策略 - 使用多种策略寻找最优方案
################################################################################
cat("\n4. 测试不同的基因过滤策略...\n")

# 策略1: 极宽松过滤（总counts >= 1）
keep_very_loose <- rowSums(counts_matrix) >= 1
cat(sprintf("   策略1 - 极宽松(总counts >= 1): %d 基因\n", sum(keep_very_loose)))

# 策略2: 宽松过滤（总counts >= 5）
keep_loose <- rowSums(counts_matrix) >= 5
cat(sprintf("   策略2 - 宽松(总counts >= 5): %d 基因\n", sum(keep_loose)))

# 策略3: 中等过滤（总counts >= 10）
keep_medium <- rowSums(counts_matrix) >= 10
cat(sprintf("   策略3 - 中等(总counts >= 10): %d 基因\n", sum(keep_medium)))

# 策略4: 在至少2个样本中有表达
keep_expressed <- rowSums(counts_matrix > 0) >= 2
cat(sprintf("   策略4 - 至少2个样本表达: %d 基因\n", sum(keep_expressed)))

# 策略5: 在至少一组中有平均表达
mod_mean <- rowMeans(counts_matrix[, group == "Mod"])
zyd_mean <- rowMeans(counts_matrix[, group == "ZYD"])
keep_group_expressed <- (mod_mean >= 1) | (zyd_mean >= 1)
cat(sprintf("   策略5 - 至少一组平均表达 >= 1: %d 基因\n", sum(keep_group_expressed)))

# 使用策略2作为主要过滤方法
cat("\n   选择策略2（宽松过滤）进行分析\n")

################################################################################
# 5. 创建DESeq2对象并运行分析
################################################################################
cat("\n5. 创建DESeq2对象并运行差异分析...\n")

dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix[keep_loose, ],
  colData = sample_info,
  design = ~ group
)

cat(sprintf("   分析基因数: %d\n", nrow(dds)))

# 运行DESeq2
dds <- DESeq(dds)

# 查看标准化因子
sf <- sizeFactors(dds)
cat("\n   标准化因子:\n")
print(data.frame(Sample = names(sf), SizeFactor = round(sf, 3)))

################################################################################
# 6. 提取结果并测试不同阈值
################################################################################
cat("\n6. 提取差异分析结果...\n")

# 提取结果 (ZYD vs Mod)
res <- results(dds, contrast = c("group", "ZYD", "Mod"))
res_df <- as.data.frame(res)

# 查看p值分布
cat("\n   P值分布:\n")
p_breaks <- c(0, 0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 1)
p_hist <- cut(res$pvalue, breaks = p_breaks)
print(table(p_hist))

# 查看adjusted p值分布
cat("\n   Adjusted P值分布:\n")
padj_hist <- cut(res$padj, breaks = p_breaks)
print(table(padj_hist, useNA = "ifany"))

################################################################################
# 7. 使用宽松阈值筛选差异基因
################################################################################
cat("\n7. 使用不同阈值筛选差异基因...\n")

# 测试多种阈值组合
thresholds <- list(
  "目标阈值 (padj<0.25, |FC|>0.5)" = list(padj = 0.25, lfc = 0.5),
  "更宽松1 (padj<0.3, |FC|>0.5)" = list(padj = 0.3, lfc = 0.5),
  "更宽松2 (padj<0.25, |FC|>0.3)" = list(padj = 0.25, lfc = 0.3),
  "极宽松 (padj<0.3, |FC|>0.3)" = list(padj = 0.3, lfc = 0.3),
  "仅p值 (p<0.05, |FC|>0.5)" = list(p = 0.05, lfc = 0.5),
  "仅p值宽松 (p<0.1, |FC|>0.3)" = list(p = 0.1, lfc = 0.3)
)

threshold_summary <- data.frame()

for (name in names(thresholds)) {
  th <- thresholds[[name]]
  
  if (!is.null(th$padj)) {
    # 使用adjusted p-value
    sig <- res_df[!is.na(res_df$padj) & 
                  res_df$padj < th$padj & 
                  abs(res_df$log2FoldChange) > th$lfc, ]
  } else {
    # 使用raw p-value
    sig <- res_df[!is.na(res_df$pvalue) & 
                  res_df$pvalue < th$p & 
                  abs(res_df$log2FoldChange) > th$lfc, ]
  }
  
  up <- sum(sig$log2FoldChange > 0)
  down <- sum(sig$log2FoldChange < 0)
  
  cat(sprintf("\n   %s:\n", name))
  cat(sprintf("     总差异基因: %d\n", nrow(sig)))
  cat(sprintf("     上调: %d, 下调: %d\n", up, down))
  
  threshold_summary <- rbind(threshold_summary, data.frame(
    Threshold = name,
    Total_DEGs = nrow(sig),
    Upregulated = up,
    Downregulated = down
  ))
}

# 保存阈值比较结果
write.csv(threshold_summary, "threshold_comparison_relaxed.csv", row.names = FALSE)

################################################################################
# 8. 使用目标阈值提取差异基因
################################################################################
cat("\n8. 使用目标阈值 (padj < 0.25, |log2FC| > 0.5) 提取差异基因...\n")

padj_threshold <- 0.25
log2fc_threshold <- 0.5

# 主要结果
sig_genes <- res_df[!is.na(res_df$padj) & 
                    res_df$padj < padj_threshold & 
                    abs(res_df$log2FoldChange) > log2fc_threshold, ]
sig_genes <- sig_genes[order(sig_genes$padj), ]

cat(sprintf("\n   ✓ 找到 %d 个差异基因\n", nrow(sig_genes)))

# 如果差异基因少于100个，尝试更宽松的阈值
if (nrow(sig_genes) < 100) {
  cat("\n   差异基因少于100个，尝试更宽松的阈值...\n")
  
  # 方案A: 放宽padj
  sig_genes_A <- res_df[!is.na(res_df$padj) & 
                        res_df$padj < 0.3 & 
                        abs(res_df$log2FoldChange) > 0.5, ]
  cat(sprintf("   方案A (padj<0.3, |FC|>0.5): %d 个基因\n", nrow(sig_genes_A)))
  
  # 方案B: 降低fold change
  sig_genes_B <- res_df[!is.na(res_df$padj) & 
                        res_df$padj < 0.25 & 
                        abs(res_df$log2FoldChange) > 0.3, ]
  cat(sprintf("   方案B (padj<0.25, |FC|>0.3): %d 个基因\n", nrow(sig_genes_B)))
  
  # 方案C: 使用原始p值
  sig_genes_C <- res_df[!is.na(res_df$pvalue) & 
                        res_df$pvalue < 0.05 & 
                        abs(res_df$log2FoldChange) > 0.5, ]
  cat(sprintf("   方案C (p<0.05, |FC|>0.5): %d 个基因\n", nrow(sig_genes_C)))
  
  # 选择基因数最接近100的方案
  if (nrow(sig_genes_A) >= 100) {
    sig_genes <- sig_genes_A[order(sig_genes_A$padj), ]
    cat("\n   采用方案A\n")
    padj_threshold <- 0.3
  } else if (nrow(sig_genes_B) >= 100) {
    sig_genes <- sig_genes_B[order(sig_genes_B$padj), ]
    cat("\n   采用方案B\n")
    log2fc_threshold <- 0.3
  } else if (nrow(sig_genes_C) >= 100) {
    sig_genes <- sig_genes_C[order(sig_genes_C$pvalue), ]
    cat("\n   采用方案C (使用原始p值)\n")
  }
}

# 统计上调和下调基因
up_genes <- sig_genes[sig_genes$log2FoldChange > 0, ]
down_genes <- sig_genes[sig_genes$log2FoldChange < 0, ]

cat(sprintf("\n   最终结果:\n"))
cat(sprintf("   总差异基因: %d\n", nrow(sig_genes)))
cat(sprintf("   上调基因: %d\n", nrow(up_genes)))
cat(sprintf("   下调基因: %d\n", nrow(down_genes)))

################################################################################
# 9. 保存差异基因列表
################################################################################
cat("\n9. 保存结果文件...\n")

# 保存所有结果
write.csv(res_df[order(res_df$padj), ], 
          "DESeq2_all_results_relaxed.csv", 
          row.names = TRUE)

# 保存显著差异基因
write.csv(sig_genes, 
          "significant_genes_relaxed.csv", 
          row.names = TRUE)

# 保存上调和下调基因
write.csv(up_genes, 
          "upregulated_genes_relaxed.csv", 
          row.names = TRUE)
write.csv(down_genes, 
          "downregulated_genes_relaxed.csv", 
          row.names = TRUE)

# 保存基因列表（仅基因名）
write.table(rownames(sig_genes), 
            "DEG_list_relaxed.txt", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)

################################################################################
# 10. 探索性分析 - 查找更多潜在的差异基因
################################################################################
cat("\n10. 探索性分析...\n")

# Top基因（按p值）
top50_pvalue <- head(res_df[order(res_df$pvalue), ], 50)
write.csv(top50_pvalue, "top50_by_pvalue_relaxed.csv", row.names = TRUE)

# Top基因（按fold change）
top50_fc <- head(res_df[order(-abs(res_df$log2FoldChange)), ], 50)
write.csv(top50_fc, "top50_by_foldchange_relaxed.csv", row.names = TRUE)

# 边缘显著基因（可能有生物学意义）
marginal_genes <- res_df[!is.na(res_df$padj) & 
                        res_df$padj >= padj_threshold & 
                        res_df$padj < 0.5 & 
                        abs(res_df$log2FoldChange) > 0.5, ]
marginal_genes <- marginal_genes[order(marginal_genes$padj), ]
write.csv(marginal_genes, "marginal_significant_genes.csv", row.names = TRUE)
cat(sprintf("   边缘显著基因 (padj 0.25-0.5): %d 个\n", nrow(marginal_genes)))

################################################################################
# 11. 生成可视化图形
################################################################################
cat("\n11. 生成可视化图形...\n")

# 11.1 PCA图 - 质量控制对比
cat("   生成PCA图（质控前后对比）...\n")

# PCA图 - 排除离群样本前（所有样本）
cat("     - 生成质控前PCA图（所有样本）\n")
sample_names_all <- colnames(counts_matrix_all)
group_all <- ifelse(grepl("^Mod-", sample_names_all), "Mod", "ZYD")
sample_info_all <- data.frame(
  sample = sample_names_all,
  group = factor(group_all, levels = c("Mod", "ZYD")),
  row.names = sample_names_all
)

# 创建临时DESeq2对象用于VST转换
keep_loose_all <- rowSums(counts_matrix_all) >= 5
dds_all <- DESeqDataSetFromMatrix(
  countData = counts_matrix_all[keep_loose_all, ],
  colData = sample_info_all,
  design = ~ group
)
vsd_all <- vst(dds_all, blind = TRUE)

# 计算PCA
pca_data_all <- plotPCA(vsd_all, intgroup = "group", returnData = TRUE)
percentVar_all <- round(100 * attr(pca_data_all, "percentVar"))

# 绘制PCA图（质控前）
pdf("PCA_plot_before_QC.pdf", width = 10, height = 8)
p_before <- ggplot(pca_data_all, aes(x = PC1, y = PC2, color = group, label = name)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  scale_color_manual(values = c("Mod" = "#E41A1C", "ZYD" = "#377EB8")) +
  labs(title = "PCA Plot - Before Quality Control (All Samples)",
       subtitle = "Outlier samples: Mod-53, ZYD-F-70, ZYD-F-62",
       x = paste0("PC1: ", percentVar_all[1], "% variance"),
       y = paste0("PC2: ", percentVar_all[2], "% variance")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "red"))
print(p_before)
dev.off()

# PCA图 - 排除离群样本后
cat("     - 生成质控后PCA图（过滤后样本）\n")
vsd <- vst(dds, blind = TRUE)
pca_data <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pdf("PCA_plot_after_QC.pdf", width = 10, height = 8)
p_after <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group, label = name)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  scale_color_manual(values = c("Mod" = "#E41A1C", "ZYD" = "#377EB8")) +
  labs(title = "PCA Plot - After Quality Control",
       subtitle = sprintf("Retained %d samples (Mod: %d, ZYD: %d)",
                         ncol(counts_matrix),
                         sum(sample_info$group == "Mod"),
                         sum(sample_info$group == "ZYD")),
       x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "darkgreen"))
print(p_after)
dev.off()

# 11.3 火山图
cat("   生成火山图...\n")
volcano_data <- res_df
volcano_data$significant <- "NO"
volcano_data$significant[!is.na(volcano_data$padj) & 
                        volcano_data$padj < padj_threshold & 
                        volcano_data$log2FoldChange > log2fc_threshold] <- "UP"
volcano_data$significant[!is.na(volcano_data$padj) & 
                        volcano_data$padj < padj_threshold & 
                        volcano_data$log2FoldChange < -log2fc_threshold] <- "DOWN"

# 标记top基因
top_genes_label <- rbind(
  head(sig_genes[order(-sig_genes$log2FoldChange), ], 10),
  head(sig_genes[order(sig_genes$log2FoldChange), ], 10)
)
volcano_data$label <- ifelse(rownames(volcano_data) %in% rownames(top_genes_label),
                            rownames(volcano_data), "")

pdf("volcano_plot_relaxed.pdf", width = 10, height = 8)
p <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = significant), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("UP" = "red", "DOWN" = "blue", "NO" = "grey70"),
                    name = "Expression") +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 30) +
  geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), 
             col = "black", linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = -log10(padj_threshold), 
             col = "black", linetype = "dashed", linewidth = 0.5) +
  labs(title = sprintf("Volcano Plot: ZYD vs Mod (padj < %.2f, |log2FC| > %.1f)", 
                      padj_threshold, log2fc_threshold),
       subtitle = sprintf("Total DEGs: %d (Up: %d, Down: %d)", 
                         nrow(sig_genes), nrow(up_genes), nrow(down_genes)),
       x = "Log2 Fold Change",
       y = "-Log10 P-value") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")
print(p)
dev.off()

# 11.4 MA图
cat("   生成MA图...\n")
pdf("MA_plot_relaxed.pdf", width = 10, height = 8)
plotMA(res, ylim = c(-5, 5), 
       main = sprintf("MA Plot (padj < %.2f)", padj_threshold),
       alpha = padj_threshold,
       colNonSig = "gray60",
       colSig = "red3")
dev.off()

# 11.5 差异基因热图
if (nrow(sig_genes) >= 2) {
  cat("   生成差异基因热图...\n")

  # VST转换
  vsd <- vst(dds, blind = FALSE)

  # 选择要展示的基因（最多50个以便显示基因名）
  n_genes <- min(50, nrow(sig_genes))
  if (nrow(sig_genes) > 50) {
    # 选择最显著的50个
    plot_genes <- rownames(sig_genes)[1:50]
    cat(sprintf("   选择最显著的 %d 个基因用于热图展示\n", n_genes))
  } else {
    plot_genes <- rownames(sig_genes)
  }

  # 提取表达矩阵并标准化
  mat <- assay(vsd)[plot_genes, ]
  mat_scaled <- t(scale(t(mat)))

  # 样本注释
  annotation_col <- data.frame(
    Group = sample_info$group,
    row.names = colnames(mat)
  )

  # 配色
  ann_colors <- list(Group = c(Mod = "#E41A1C", ZYD = "#377EB8"))

  # 生成热图 - 显示基因名称
  pdf("heatmap_DEGs_relaxed.pdf", width = 10, height = 14)
  pheatmap(
    mat_scaled,
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,  # 始终显示基因名
    show_colnames = TRUE,
    scale = "none",
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    main = sprintf("Top %d Differentially Expressed Genes\n(After removing outlier samples)", n_genes),
    fontsize = 10,
    fontsize_row = 8,  # 增大字体以便阅读
    cellwidth = 20,
    cellheight = 10
  )
  dev.off()

  # 如果差异基因多，再生成一个包含所有基因但不显示名称的版本
  if (nrow(sig_genes) > 50) {
    cat(sprintf("   生成包含所有 %d 个差异基因的热图（不显示基因名）\n", nrow(sig_genes)))
    mat_all <- assay(vsd)[rownames(sig_genes), ]
    mat_all_scaled <- t(scale(t(mat_all)))

    pdf("heatmap_all_DEGs_relaxed.pdf", width = 10, height = 12)
    pheatmap(
      mat_all_scaled,
      annotation_col = annotation_col,
      annotation_colors = ann_colors,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_rownames = FALSE,
      show_colnames = TRUE,
      scale = "none",
      color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
      main = sprintf("All %d Differentially Expressed Genes\n(After removing outlier samples)", nrow(sig_genes)),
      fontsize = 10
    )
    dev.off()
  }
}

# 11.6 表达量箱线图（top差异基因）
cat("   生成top差异基因表达箱线图...\n")
top_degs <- rownames(sig_genes)[1:min(20, nrow(sig_genes))]
normalized_counts <- counts(dds, normalized = TRUE)

pdf("top_DEGs_boxplot_relaxed.pdf", width = 12, height = 10)
par(mfrow = c(4, 5), mar = c(4, 4, 3, 1))
for (gene in top_degs) {
  gene_counts <- normalized_counts[gene, ]
  boxplot(gene_counts ~ sample_info$group,
          col = c("#E41A1C", "#377EB8"),
          main = gene,
          ylab = "Normalized counts",
          xlab = "",
          cex.main = 0.8)
}
dev.off()

################################################################################
# 12. 生成分析报告
################################################################################
cat("\n12. 生成分析报告...\n")

report <- sprintf("
================================================================================
DESeq2 差异表达分析报告 - 宽松阈值版本（排除离群样本）
================================================================================
分析日期: %s
使用阈值: padj < %.2f, |log2FC| > %.1f

样本质量控制:
- 排除离群样本: Mod-53, ZYD-F-70, ZYD-F-62
- 排除原因: PCA分析显示这些样本严重偏离各自组的中心

样本信息:
- 便秘模型组 (Mod): %d 个样本
- 12味复方组 (ZYD): %d 个样本

基因过滤:
- 原始基因数: %d
- 过滤后基因数: %d (使用宽松过滤: rowSums >= 5)

差异基因统计:
- 总差异基因数: %d
- 上调基因: %d (%.1f%%)
- 下调基因: %d (%.1f%%)

其他统计:
- padj < 0.05 的基因: %d
- padj < 0.1 的基因: %d
- padj < 0.25 的基因: %d
- 边缘显著基因 (0.25 < padj < 0.5): %d

输出文件:
1. DESeq2_all_results_relaxed.csv - 所有基因结果
2. significant_genes_relaxed.csv - 显著差异基因
3. upregulated_genes_relaxed.csv - 上调基因
4. downregulated_genes_relaxed.csv - 下调基因
5. DEG_list_relaxed.txt - 差异基因列表
6. threshold_comparison_relaxed.csv - 阈值比较
7. top50_by_pvalue_relaxed.csv - Top50基因(按p值)
8. top50_by_foldchange_relaxed.csv - Top50基因(按FC)
9. marginal_significant_genes.csv - 边缘显著基因

可视化文件:
10. PCA_plot_before_QC.pdf - PCA图（质控前，所有样本）
11. PCA_plot_after_QC.pdf - PCA图（质控后，过滤样本）
12. volcano_plot_relaxed.pdf - 火山图
13. MA_plot_relaxed.pdf - MA图
14. heatmap_DEGs_relaxed.pdf - 差异基因热图
15. top_DEGs_boxplot_relaxed.pdf - Top差异基因箱线图

注意事项:
- 使用了宽松的阈值以获得更多差异基因
- 建议对关键基因进行实验验证
- 可进行后续的功能富集分析

================================================================================
", 
  Sys.Date(),
  padj_threshold,
  log2fc_threshold,
  sum(sample_info$group == "Mod"),
  sum(sample_info$group == "ZYD"),
  nrow(counts_matrix),
  sum(keep_loose),
  nrow(sig_genes),
  nrow(up_genes),
  nrow(up_genes)/nrow(sig_genes)*100,
  nrow(down_genes),
  nrow(down_genes)/nrow(sig_genes)*100,
  sum(!is.na(res_df$padj) & res_df$padj < 0.05),
  sum(!is.na(res_df$padj) & res_df$padj < 0.1),
  sum(!is.na(res_df$padj) & res_df$padj < 0.25),
  nrow(marginal_genes)
)

cat(report)
writeLines(report, "analysis_report_relaxed.txt")

cat("\n✓ 分析完成！\n")
cat(sprintf("✓ 成功找到 %d 个差异基因\n", nrow(sig_genes)))
cat("================================================================================\n")
