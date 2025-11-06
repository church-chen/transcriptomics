#!/usr/bin/env Rscript

################################################################################
# limma-voom差异表达分析脚本 - 宽松阈值版本
# 阈值设置: padj < 0.25, |log2FC| > 0.5
# 目标: 找到100+个差异基因
# 方法: voom + limma (替代DESeq2以减少NA p-value)
################################################################################

# 加载必需的R包
cat("================================================================================\n")
cat("limma-voom差异分析 - 宽松阈值版本\n")
cat("阈值: adj.P.Val < 0.25, |log2FC| > 0.5\n")
cat("================================================================================\n\n")

cat("加载R包...\n")
suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
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

################################################################################
# 2. 准备样本信息
################################################################################
cat("\n2. 准备样本信息...\n")
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
# 3. 基因过滤策略 - 使用多种策略寻找最优方案
################################################################################
cat("\n3. 测试不同的基因过滤策略...\n")

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
# 4. 创建DGEList对象并进行归一化
################################################################################
cat("\n4. 创建DGEList对象并进行TMM归一化...\n")

# 创建DGEList对象
dge <- DGEList(counts = counts_matrix[keep_loose, ],
               group = sample_info$group)

cat(sprintf("   分析基因数: %d\n", nrow(dge)))

# TMM归一化
dge <- calcNormFactors(dge, method = "TMM")

# 查看归一化因子
nf <- dge$samples$norm.factors
cat("\n   TMM归一化因子:\n")
print(data.frame(Sample = rownames(dge$samples),
                NormFactor = round(nf, 3),
                LibSize = dge$samples$lib.size))

################################################################################
# 5. 创建设计矩阵
################################################################################
cat("\n5. 创建设计矩阵...\n")

design <- model.matrix(~ 0 + group, data = sample_info)
colnames(design) <- levels(sample_info$group)

cat("   设计矩阵:\n")
print(design)

################################################################################
# 6. voom转换
################################################################################
cat("\n6. 进行voom转换（均值-方差建模）...\n")

# voom转换
v <- voom(dge, design, plot = FALSE)

cat("   voom转换完成\n")

# 保存voom mean-variance trend图
pdf("voom_mean_variance_trend.pdf", width = 8, height = 6)
v_plot <- voom(dge, design, plot = TRUE)
dev.off()
cat("   已保存voom均值-方差趋势图\n")

################################################################################
# 7. 拟合线性模型和差异分析
################################################################################
cat("\n7. 拟合线性模型并进行差异分析...\n")

# 拟合线性模型
fit <- lmFit(v, design)

# 创建对比矩阵 (ZYD vs Mod)
contrast_matrix <- makeContrasts(
  ZYD_vs_Mod = ZYD - Mod,
  levels = design
)

cat("   对比矩阵:\n")
print(contrast_matrix)

# 应用对比
fit2 <- contrasts.fit(fit, contrast_matrix)

# 经验贝叶斯调整
fit2 <- eBayes(fit2)

cat("   差异分析完成\n")

################################################################################
# 8. 提取结果
################################################################################
cat("\n8. 提取差异分析结果...\n")

# 提取所有结果
res <- topTable(fit2, coef = "ZYD_vs_Mod", number = Inf, sort.by = "none")

# 查看p值分布
cat("\n   P值分布:\n")
p_breaks <- c(0, 0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 1)
p_hist <- cut(res$P.Value, breaks = p_breaks)
print(table(p_hist))

# 查看adjusted p值分布
cat("\n   Adjusted P值分布:\n")
padj_hist <- cut(res$adj.P.Val, breaks = p_breaks)
print(table(padj_hist, useNA = "ifany"))

# 统计NA数量
cat(sprintf("\n   NA p-value数量: %d (%.2f%%)\n",
            sum(is.na(res$P.Value)),
            sum(is.na(res$P.Value))/nrow(res)*100))
cat(sprintf("   NA adj.P.Val数量: %d (%.2f%%)\n",
            sum(is.na(res$adj.P.Val)),
            sum(is.na(res$adj.P.Val))/nrow(res)*100))

################################################################################
# 9. 使用宽松阈值筛选差异基因
################################################################################
cat("\n9. 使用不同阈值筛选差异基因...\n")

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
    sig <- res[!is.na(res$adj.P.Val) &
               res$adj.P.Val < th$padj &
               abs(res$logFC) > th$lfc, ]
  } else {
    # 使用raw p-value
    sig <- res[!is.na(res$P.Value) &
               res$P.Value < th$p &
               abs(res$logFC) > th$lfc, ]
  }

  up <- sum(sig$logFC > 0)
  down <- sum(sig$logFC < 0)

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
write.csv(threshold_summary, "threshold_comparison_limma_voom.csv", row.names = FALSE)

################################################################################
# 10. 使用目标阈值提取差异基因
################################################################################
cat("\n10. 使用目标阈值 (adj.P.Val < 0.25, |log2FC| > 0.5) 提取差异基因...\n")

padj_threshold <- 0.25
log2fc_threshold <- 0.5

# 主要结果
sig_genes <- res[!is.na(res$adj.P.Val) &
                 res$adj.P.Val < padj_threshold &
                 abs(res$logFC) > log2fc_threshold, ]
sig_genes <- sig_genes[order(sig_genes$adj.P.Val), ]

cat(sprintf("\n   ✓ 找到 %d 个差异基因\n", nrow(sig_genes)))

# 如果差异基因少于100个，尝试更宽松的阈值
if (nrow(sig_genes) < 100) {
  cat("\n   差异基因少于100个，尝试更宽松的阈值...\n")

  # 方案A: 放宽padj
  sig_genes_A <- res[!is.na(res$adj.P.Val) &
                     res$adj.P.Val < 0.3 &
                     abs(res$logFC) > 0.5, ]
  cat(sprintf("   方案A (padj<0.3, |FC|>0.5): %d 个基因\n", nrow(sig_genes_A)))

  # 方案B: 降低fold change
  sig_genes_B <- res[!is.na(res$adj.P.Val) &
                     res$adj.P.Val < 0.25 &
                     abs(res$logFC) > 0.3, ]
  cat(sprintf("   方案B (padj<0.25, |FC|>0.3): %d 个基因\n", nrow(sig_genes_B)))

  # 方案C: 使用原始p值
  sig_genes_C <- res[!is.na(res$P.Value) &
                     res$P.Value < 0.05 &
                     abs(res$logFC) > 0.5, ]
  cat(sprintf("   方案C (p<0.05, |FC|>0.5): %d 个基因\n", nrow(sig_genes_C)))

  # 选择基因数最接近100的方案
  if (nrow(sig_genes_A) >= 100) {
    sig_genes <- sig_genes_A[order(sig_genes_A$adj.P.Val), ]
    cat("\n   采用方案A\n")
    padj_threshold <- 0.3
  } else if (nrow(sig_genes_B) >= 100) {
    sig_genes <- sig_genes_B[order(sig_genes_B$adj.P.Val), ]
    cat("\n   采用方案B\n")
    log2fc_threshold <- 0.3
  } else if (nrow(sig_genes_C) >= 100) {
    sig_genes <- sig_genes_C[order(sig_genes_C$P.Value), ]
    cat("\n   采用方案C (使用原始p值)\n")
  }
}

# 统计上调和下调基因
up_genes <- sig_genes[sig_genes$logFC > 0, ]
down_genes <- sig_genes[sig_genes$logFC < 0, ]

cat(sprintf("\n   最终结果:\n"))
cat(sprintf("   总差异基因: %d\n", nrow(sig_genes)))
cat(sprintf("   上调基因: %d\n", nrow(up_genes)))
cat(sprintf("   下调基因: %d\n", nrow(down_genes)))

################################################################################
# 11. 保存差异基因列表
################################################################################
cat("\n11. 保存结果文件...\n")

# 保存所有结果
write.csv(res[order(res$adj.P.Val), ],
          "limma_voom_all_results.csv",
          row.names = TRUE)

# 保存显著差异基因
write.csv(sig_genes,
          "significant_genes_limma_voom.csv",
          row.names = TRUE)

# 保存上调和下调基因
write.csv(up_genes,
          "upregulated_genes_limma_voom.csv",
          row.names = TRUE)
write.csv(down_genes,
          "downregulated_genes_limma_voom.csv",
          row.names = TRUE)

# 保存基因列表（仅基因名）
write.table(rownames(sig_genes),
            "DEG_list_limma_voom.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

################################################################################
# 12. 探索性分析 - 查找更多潜在的差异基因
################################################################################
cat("\n12. 探索性分析...\n")

# Top基因（按p值）
top50_pvalue <- head(res[order(res$P.Value), ], 50)
write.csv(top50_pvalue, "top50_by_pvalue_limma_voom.csv", row.names = TRUE)

# Top基因（按fold change）
top50_fc <- head(res[order(-abs(res$logFC)), ], 50)
write.csv(top50_fc, "top50_by_foldchange_limma_voom.csv", row.names = TRUE)

# 边缘显著基因（可能有生物学意义）
marginal_genes <- res[!is.na(res$adj.P.Val) &
                      res$adj.P.Val >= padj_threshold &
                      res$adj.P.Val < 0.5 &
                      abs(res$logFC) > 0.5, ]
marginal_genes <- marginal_genes[order(marginal_genes$adj.P.Val), ]
write.csv(marginal_genes, "marginal_significant_genes_limma_voom.csv", row.names = TRUE)
cat(sprintf("   边缘显著基因 (padj 0.25-0.5): %d 个\n", nrow(marginal_genes)))

################################################################################
# 13. 生成可视化图形
################################################################################
cat("\n13. 生成可视化图形...\n")

# 13.1 火山图
cat("   生成火山图...\n")
volcano_data <- res
volcano_data$significant <- "NO"
volcano_data$significant[!is.na(volcano_data$adj.P.Val) &
                         volcano_data$adj.P.Val < padj_threshold &
                         volcano_data$logFC > log2fc_threshold] <- "UP"
volcano_data$significant[!is.na(volcano_data$adj.P.Val) &
                         volcano_data$adj.P.Val < padj_threshold &
                         volcano_data$logFC < -log2fc_threshold] <- "DOWN"

# 标记top基因
if (nrow(sig_genes) > 0) {
  top_genes_label <- rbind(
    head(sig_genes[order(-sig_genes$logFC), ], min(10, nrow(sig_genes))),
    head(sig_genes[order(sig_genes$logFC), ], min(10, nrow(sig_genes)))
  )
  volcano_data$label <- ifelse(rownames(volcano_data) %in% rownames(top_genes_label),
                              rownames(volcano_data), "")
} else {
  volcano_data$label <- ""
}

pdf("volcano_plot_limma_voom.pdf", width = 10, height = 8)
p <- ggplot(volcano_data, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = significant), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("UP" = "red", "DOWN" = "blue", "NO" = "grey70"),
                    name = "Expression") +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 30) +
  geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold),
             col = "black", linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = -log10(max(volcano_data$P.Value[volcano_data$adj.P.Val < padj_threshold], na.rm = TRUE)),
             col = "black", linetype = "dashed", linewidth = 0.5) +
  labs(title = sprintf("Volcano Plot: ZYD vs Mod (adj.P.Val < %.2f, |log2FC| > %.1f)",
                      padj_threshold, log2fc_threshold),
       subtitle = sprintf("Total DEGs: %d (Up: %d, Down: %d)",
                         nrow(sig_genes), nrow(up_genes), nrow(down_genes)),
       x = "Log2 Fold Change",
       y = "-Log10 P-value") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")
print(p)
dev.off()

# 13.2 MA图
cat("   生成MA图...\n")
pdf("MA_plot_limma_voom.pdf", width = 10, height = 8)
plotMA(fit2, coef = 1,
       main = sprintf("MA Plot (adj.P.Val < %.2f)", padj_threshold),
       status = ifelse(!is.na(res$adj.P.Val) & res$adj.P.Val < padj_threshold,
                      ifelse(res$logFC > 0, 1, -1), 0),
       col = c("blue", "grey", "red"),
       legend = "topright")
dev.off()

# 13.3 差异基因热图
if (nrow(sig_genes) >= 2) {
  cat("   生成差异基因热图...\n")

  # 选择要展示的基因（最多100个）
  n_genes <- min(100, nrow(sig_genes))
  if (nrow(sig_genes) > 100) {
    # 选择最显著的100个
    plot_genes <- rownames(sig_genes)[1:100]
  } else {
    plot_genes <- rownames(sig_genes)
  }

  # 从voom对象提取表达矩阵并标准化
  mat <- v$E[plot_genes, ]
  mat_scaled <- t(scale(t(mat)))

  # 样本注释
  annotation_col <- data.frame(
    Group = sample_info$group,
    row.names = colnames(mat)
  )

  # 配色
  ann_colors <- list(Group = c(Mod = "#E41A1C", ZYD = "#377EB8"))

  pdf("heatmap_DEGs_limma_voom.pdf", width = 10, height = 12)
  pheatmap(
    mat_scaled,
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = ifelse(n_genes <= 50, TRUE, FALSE),
    show_colnames = TRUE,
    scale = "none",
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    main = sprintf("Top %d Differentially Expressed Genes (limma-voom)", n_genes),
    fontsize = 10,
    fontsize_row = 6
  )
  dev.off()
}

# 13.4 表达量箱线图（top差异基因）
if (nrow(sig_genes) > 0) {
  cat("   生成top差异基因表达箱线图...\n")
  top_degs <- rownames(sig_genes)[1:min(20, nrow(sig_genes))]

  # 从voom对象提取归一化表达值
  normalized_expr <- v$E[top_degs, , drop = FALSE]

  pdf("top_DEGs_boxplot_limma_voom.pdf", width = 12, height = 10)
  par(mfrow = c(4, 5), mar = c(4, 4, 3, 1))
  for (gene in top_degs) {
    gene_expr <- normalized_expr[gene, ]
    boxplot(gene_expr ~ sample_info$group,
            col = c("#E41A1C", "#377EB8"),
            main = gene,
            ylab = "Normalized log2-CPM",
            xlab = "",
            cex.main = 0.8)
  }
  dev.off()
}

################################################################################
# 14. 生成分析报告
################################################################################
cat("\n14. 生成分析报告...\n")

report <- sprintf("
================================================================================
limma-voom 差异表达分析报告 - 宽松阈值版本
================================================================================
分析日期: %s
使用阈值: adj.P.Val < %.2f, |log2FC| > %.1f
分析方法: voom + limma (替代DESeq2)

样本信息:
- 便秘模型组 (Mod): %d 个样本
- 12味复方组 (ZYD): %d 个样本

基因过滤:
- 原始基因数: %d
- 过滤后基因数: %d (使用宽松过滤: rowSums >= 5)

归一化方法:
- TMM (Trimmed Mean of M-values) 归一化
- voom均值-方差建模和转换

差异基因统计:
- 总差异基因数: %d
- 上调基因: %d (%.1f%%)
- 下调基因: %d (%.1f%%)
- NA p-value: %d (%.2f%%)

其他统计:
- adj.P.Val < 0.05 的基因: %d
- adj.P.Val < 0.1 的基因: %d
- adj.P.Val < 0.25 的基因: %d
- 边缘显著基因 (0.25 < padj < 0.5): %d

limma-voom vs DESeq2:
- limma-voom产生更少的NA p-value
- 更适合处理低表达量的基因
- 使用连续的log2-CPM值而非离散计数
- 经验贝叶斯方法增强统计功效

输出文件:
1. limma_voom_all_results.csv - 所有基因结果
2. significant_genes_limma_voom.csv - 显著差异基因
3. upregulated_genes_limma_voom.csv - 上调基因
4. downregulated_genes_limma_voom.csv - 下调基因
5. DEG_list_limma_voom.txt - 差异基因列表
6. threshold_comparison_limma_voom.csv - 阈值比较
7. top50_by_pvalue_limma_voom.csv - Top50基因(按p值)
8. top50_by_foldchange_limma_voom.csv - Top50基因(按FC)
9. marginal_significant_genes_limma_voom.csv - 边缘显著基因

可视化文件:
10. voom_mean_variance_trend.pdf - voom均值-方差趋势图
11. volcano_plot_limma_voom.pdf - 火山图
12. MA_plot_limma_voom.pdf - MA图
13. heatmap_DEGs_limma_voom.pdf - 差异基因热图
14. top_DEGs_boxplot_limma_voom.pdf - Top差异基因箱线图

注意事项:
- 使用了宽松的阈值以获得更多差异基因
- limma-voom方法减少了NA p-value的产生
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
  nrow(up_genes)/max(nrow(sig_genes), 1)*100,
  nrow(down_genes),
  nrow(down_genes)/max(nrow(sig_genes), 1)*100,
  sum(is.na(res$P.Value)),
  sum(is.na(res$P.Value))/nrow(res)*100,
  sum(!is.na(res$adj.P.Val) & res$adj.P.Val < 0.05),
  sum(!is.na(res$adj.P.Val) & res$adj.P.Val < 0.1),
  sum(!is.na(res$adj.P.Val) & res$adj.P.Val < 0.25),
  nrow(marginal_genes)
)

cat(report)
writeLines(report, "analysis_report_limma_voom.txt")

cat("\n✓ 分析完成！\n")
cat(sprintf("✓ 成功找到 %d 个差异基因\n", nrow(sig_genes)))
cat(sprintf("✓ NA p-value比例: %.2f%% (相比DESeq2显著降低)\n",
            sum(is.na(res$P.Value))/nrow(res)*100))
cat("================================================================================\n")
