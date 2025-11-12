#!/usr/bin/env Rscript

################################################################################
# GO/KEGG/GSEA 富集分析脚本 - 小鼠（改进版）
# 输入：DESeq2差异分析结果
# 输出：ORA + GSEA富集分析结果和可视化图表
# 改进要点：
# 1. 增加GSEA分析（使用全部基因）
# 2. 同时支持严格阈值（padj<0.05）和宽松阈值（padj<0.25）
# 3. 统一使用ENTREZID进行富集分析
# 4. 增加更多可视化选项
################################################################################

cat("================================================================================\n")
cat("GO/KEGG/GSEA 富集分析 - 小鼠（改进版）\n")
cat("================================================================================\n\n")

# 加载必需的R包
cat("加载R包...\n")
suppressPackageStartupMessages({
  library(DESeq2)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(DOSE)
})

################################################################################
# 1. 读取DESeq2完整结果
################################################################################
cat("\n1. 读取DESeq2结果...\n")

# 读取DESeq2完整结果（包含所有基因）
# 假设你的DESeq2结果保存为RDS文件
# 如果是CSV文件，用 read.csv() 读取
tryCatch({
  # 优先读取RDS文件
  deseq_results <- readRDS("deseq2_results.rds")
  cat("   成功读取RDS文件\n")
}, error = function(e) {
  # 备选：读取CSV文件
  cat("   RDS文件不存在，尝试读取CSV文件...\n")
  deseq_results <- read.csv("deseq2_results.csv", row.names = 1)
})

cat(sprintf("   总基因数: %d\n", nrow(deseq_results)))

# 移除NA值
deseq_results <- deseq_results[!is.na(deseq_results$padj), ]
cat(sprintf("   去除NA后基因数: %d\n", nrow(deseq_results)))

################################################################################
# 2. 基因筛选 - 使用两种阈值
################################################################################
cat("\n2. 基因筛选...\n")

# 严格阈值（推荐用于文章）
sig_genes_strict <- deseq_results[deseq_results$padj < 0.05 & 
                                   abs(deseq_results$log2FoldChange) > 1, ]
cat(sprintf("   严格阈值 (padj<0.05, |log2FC|>1): %d 个DEGs\n", nrow(sig_genes_strict)))
cat(sprintf("     上调: %d, 下调: %d\n", 
            sum(sig_genes_strict$log2FoldChange > 0),
            sum(sig_genes_strict$log2FoldChange < 0)))

# 宽松阈值（探索性分析）
sig_genes_relaxed <- deseq_results[deseq_results$padj < 0.25 & 
                                    abs(deseq_results$log2FoldChange) > 0.3, ]
cat(sprintf("   宽松阈值 (padj<0.25, |log2FC|>0.3): %d 个DEGs\n", nrow(sig_genes_relaxed)))
cat(sprintf("     上调: %d, 下调: %d\n", 
            sum(sig_genes_relaxed$log2FoldChange > 0),
            sum(sig_genes_relaxed$log2FoldChange < 0)))

# 根据基因数量选择使用哪个阈值进行ORA分析
if (nrow(sig_genes_strict) < 10) {
  cat("\n   ⚠️  严格阈值基因数太少，ORA分析使用宽松阈值\n")
  sig_genes <- sig_genes_relaxed
  threshold_label <- "relaxed"
} else {
  cat("\n   ✓ ORA分析使用严格阈值\n")
  sig_genes <- sig_genes_strict
  threshold_label <- "strict"
}

# 保存筛选后的DEGs
write.csv(sig_genes, paste0("significant_genes_", threshold_label, ".csv"))

################################################################################
# 3. 基因ID转换（ENSEMBL -> ENTREZID + SYMBOL）
################################################################################
cat("\n3. 基因ID转换...\n")

# 提取ENSEMBL ID（去除版本号）
all_ensembl <- gsub("\\..*", "", rownames(deseq_results))
sig_ensembl <- gsub("\\..*", "", rownames(sig_genes))

# 分离上调和下调基因
up_genes <- sig_genes[sig_genes$log2FoldChange > 0, ]
down_genes <- sig_genes[sig_genes$log2FoldChange < 0, ]
up_ensembl <- gsub("\\..*", "", rownames(up_genes))
down_ensembl <- gsub("\\..*", "", rownames(down_genes))

# ENSEMBL转ENTREZID和SYMBOL
cat("   转换所有基因ID...\n")
id_conversion <- bitr(all_ensembl,
                     fromType = "ENSEMBL",
                     toType = c("ENTREZID", "SYMBOL"),
                     OrgDb = org.Mm.eg.db)

cat(sprintf("   转换成功: %d/%d (%.1f%%)\n",
            nrow(id_conversion), length(all_ensembl),
            nrow(id_conversion)/length(all_ensembl)*100))

# 转换差异基因
up_entrez <- bitr(up_ensembl,
                  fromType = "ENSEMBL",
                  toType = c("ENTREZID", "SYMBOL"),
                  OrgDb = org.Mm.eg.db)

down_entrez <- bitr(down_ensembl,
                    fromType = "ENSEMBL",
                    toType = c("ENTREZID", "SYMBOL"),
                    OrgDb = org.Mm.eg.db)

cat(sprintf("   上调基因转换成功: %d/%d\n", nrow(up_entrez), length(up_ensembl)))
cat(sprintf("   下调基因转换成功: %d/%d\n", nrow(down_entrez), length(down_ensembl)))

# 保存ID对应表
write.csv(id_conversion, "gene_id_conversion.csv", row.names = FALSE)

################################################################################
# 4. 准备GSEA基因列表（全部基因）
################################################################################
cat("\n4. 准备GSEA基因列表...\n")

# 合并log2FoldChange和ENTREZID
gsea_data <- data.frame(
  ensembl = all_ensembl,
  log2FC = deseq_results$log2FoldChange,
  stringsAsFactors = FALSE
)

# 合并ENTREZID
gsea_data <- merge(gsea_data, id_conversion[, c("ENSEMBL", "ENTREZID")], 
                   by.x = "ensembl", by.y = "ENSEMBL", all.x = TRUE)

# 去除NA和重复
gsea_data <- gsea_data[!is.na(gsea_data$ENTREZID) & !is.na(gsea_data$log2FC), ]
gsea_data <- gsea_data[!duplicated(gsea_data$ENTREZID), ]

# 创建排序的基因列表
gene_list <- gsea_data$log2FC
names(gene_list) <- gsea_data$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

cat(sprintf("   GSEA基因列表: %d 个基因\n", length(gene_list)))
cat(sprintf("   log2FC范围: %.2f 到 %.2f\n", min(gene_list), max(gene_list)))

################################################################################
# 5. GO富集分析 (ORA - Over-Representation Analysis)
################################################################################
cat("\n5. GO富集分析 (ORA)...\n")

# 5.1 上调基因 GO 富集
cat("   5.1 上调基因 GO 富集...\n")
go_up <- NULL
if (nrow(up_entrez) > 0) {
  go_up <- enrichGO(gene = up_entrez$ENTREZID,
                    OrgDb = org.Mm.eg.db,
                    keyType = "ENTREZID",
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2,
                    readable = TRUE)
  
  if (!is.null(go_up) && nrow(go_up) > 0) {
    cat(sprintf("     富集到 %d 个GO term\n", nrow(go_up)))
    write.csv(as.data.frame(go_up), 
              paste0("GO_ORA_upregulated_", threshold_label, ".csv"), 
              row.names = FALSE)
  } else {
    cat("     未找到显著富集的GO term\n")
  }
}

# 5.2 下调基因 GO 富集
cat("   5.2 下调基因 GO 富集...\n")
go_down <- NULL
if (nrow(down_entrez) > 0) {
  go_down <- enrichGO(gene = down_entrez$ENTREZID,
                      OrgDb = org.Mm.eg.db,
                      keyType = "ENTREZID",
                      ont = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2,
                      readable = TRUE)
  
  if (!is.null(go_down) && nrow(go_down) > 0) {
    cat(sprintf("     富集到 %d 个GO term\n", nrow(go_down)))
    write.csv(as.data.frame(go_down), 
              paste0("GO_ORA_downregulated_", threshold_label, ".csv"), 
              row.names = FALSE)
  } else {
    cat("     未找到显著富集的GO term\n")
  }
}

################################################################################
# 6. KEGG通路富集分析 (ORA)
################################################################################
cat("\n6. KEGG通路富集分析 (ORA)...\n")

# 6.1 上调基因 KEGG 富集
cat("   6.1 上调基因 KEGG 富集...\n")
kegg_up <- NULL
if (nrow(up_entrez) > 0) {
  kegg_up <- enrichKEGG(gene = up_entrez$ENTREZID,
                        organism = "mmu",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2)
  
  if (!is.null(kegg_up) && nrow(kegg_up) > 0) {
    kegg_up <- setReadable(kegg_up, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
    cat(sprintf("     富集到 %d 个KEGG通路\n", nrow(kegg_up)))
    write.csv(as.data.frame(kegg_up), 
              paste0("KEGG_ORA_upregulated_", threshold_label, ".csv"), 
              row.names = FALSE)
  } else {
    cat("     未找到显著富集的KEGG通路\n")
  }
}

# 6.2 下调基因 KEGG 富集
cat("   6.2 下调基因 KEGG 富集...\n")
kegg_down <- NULL
if (nrow(down_entrez) > 0) {
  kegg_down <- enrichKEGG(gene = down_entrez$ENTREZID,
                          organism = "mmu",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)
  
  if (!is.null(kegg_down) && nrow(kegg_down) > 0) {
    kegg_down <- setReadable(kegg_down, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
    cat(sprintf("     富集到 %d 个KEGG通路\n", nrow(kegg_down)))
    write.csv(as.data.frame(kegg_down), 
              paste0("KEGG_ORA_downregulated_", threshold_label, ".csv"), 
              row.names = FALSE)
  } else {
    cat("     未找到显著富集的KEGG通路\n")
  }
}

################################################################################
# 7. GSEA分析 - GO (Gene Set Enrichment Analysis) ⭐ 新增
################################################################################
cat("\n7. GSEA分析 - GO...\n")

gsea_go <- NULL
tryCatch({
  gsea_go <- gseGO(geneList = gene_list,
                   OrgDb = org.Mm.eg.db,
                   ont = "ALL",
                   keyType = "ENTREZID",
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   verbose = FALSE)
  
  if (!is.null(gsea_go) && nrow(gsea_go) > 0) {
    cat(sprintf("   ✓ GSEA-GO: 富集到 %d 个GO term\n", nrow(gsea_go)))
    write.csv(as.data.frame(gsea_go), "GSEA_GO_results.csv", row.names = FALSE)
  } else {
    cat("   未找到显著富集的GO term (GSEA)\n")
  }
}, error = function(e) {
  cat(sprintf("   ⚠️  GSEA-GO分析出错: %s\n", e$message))
})

################################################################################
# 8. GSEA分析 - KEGG ⭐ 新增
################################################################################
cat("\n8. GSEA分析 - KEGG...\n")

gsea_kegg <- NULL
tryCatch({
  gsea_kegg <- gseKEGG(geneList = gene_list,
                       organism = "mmu",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH",
                       verbose = FALSE)
  
  if (!is.null(gsea_kegg) && nrow(gsea_kegg) > 0) {
    cat(sprintf("   ✓ GSEA-KEGG: 富集到 %d 个KEGG通路\n", nrow(gsea_kegg)))
    # 添加基因symbol
    gsea_kegg <- setReadable(gsea_kegg, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
    write.csv(as.data.frame(gsea_kegg), "GSEA_KEGG_results.csv", row.names = FALSE)
  } else {
    cat("   未找到显著富集的KEGG通路 (GSEA)\n")
  }
}, error = function(e) {
  cat(sprintf("   ⚠️  GSEA-KEGG分析出错: %s\n", e$message))
})

################################################################################
# 9. 可视化 - ORA结果
################################################################################
cat("\n9. 生成ORA可视化图表...\n")

# 9.1 GO ORA可视化
if (!is.null(go_up) && nrow(go_up) > 0) {
  cat("   9.1 上调基因 GO 可视化...\n")
  
  pdf(paste0("GO_ORA_dotplot_upregulated_", threshold_label, ".pdf"), 
      width = 12, height = 10)
  print(dotplot(go_up, showCategory = 10, split = "ONTOLOGY", font.size = 10) +
          facet_grid(ONTOLOGY~., scales = "free"))
  dev.off()
  
  # 如果基因数量足够，画网络图
  if (nrow(go_up) >= 5 && sum(up_entrez$ENTREZID %in% unique(unlist(strsplit(go_up$geneID, "/")))) >= 3) {
    pdf(paste0("GO_ORA_cnetplot_upregulated_", threshold_label, ".pdf"), 
        width = 14, height = 12)
    print(cnetplot(go_up, showCategory = 5, node_label = "all"))
    dev.off()
  }
}

if (!is.null(go_down) && nrow(go_down) > 0) {
  cat("   9.2 下调基因 GO 可视化...\n")
  
  pdf(paste0("GO_ORA_dotplot_downregulated_", threshold_label, ".pdf"), 
      width = 12, height = 10)
  print(dotplot(go_down, showCategory = 10, split = "ONTOLOGY", font.size = 10) +
          facet_grid(ONTOLOGY~., scales = "free"))
  dev.off()
}

# 9.2 KEGG ORA可视化
if (!is.null(kegg_up) && nrow(kegg_up) > 0) {
  cat("   9.3 上调基因 KEGG 可视化...\n")
  
  pdf(paste0("KEGG_ORA_dotplot_upregulated_", threshold_label, ".pdf"), 
      width = 12, height = 8)
  print(dotplot(kegg_up, showCategory = 20, font.size = 10))
  dev.off()
  
  if (nrow(kegg_up) >= 5) {
    pdf(paste0("KEGG_ORA_cnetplot_upregulated_", threshold_label, ".pdf"), 
        width = 14, height = 12)
    print(cnetplot(kegg_up, showCategory = 5, node_label = "all"))
    dev.off()
  }
}

if (!is.null(kegg_down) && nrow(kegg_down) > 0) {
  cat("   9.4 下调基因 KEGG 可视化...\n")
  
  pdf(paste0("KEGG_ORA_dotplot_downregulated_", threshold_label, ".pdf"), 
      width = 12, height = 8)
  print(dotplot(kegg_down, showCategory = 20, font.size = 10))
  dev.off()
}

################################################################################
# 10. 可视化 - GSEA结果 ⭐ 新增
################################################################################
cat("\n10. 生成GSEA可视化图表...\n")

# 10.1 GSEA-GO可视化
if (!is.null(gsea_go) && nrow(gsea_go) > 0) {
  cat("   10.1 GSEA-GO 可视化...\n")
  
  # Dotplot
  pdf("GSEA_GO_dotplot.pdf", width = 12, height = 10)
  print(dotplot(gsea_go, showCategory = 20, split = "ONTOLOGY", font.size = 10) +
          facet_grid(ONTOLOGY~., scales = "free"))
  dev.off()
  
  # Ridge plot (展示富集分数分布)
  pdf("GSEA_GO_ridgeplot.pdf", width = 12, height = 10)
  print(ridgeplot(gsea_go, showCategory = 20))
  dev.off()
  
  # GSEA plot (展示top 3通路)
  pdf("GSEA_GO_enrichment_plot.pdf", width = 12, height = 8)
  print(gseaplot2(gsea_go, geneSetID = 1:min(3, nrow(gsea_go)), 
                  pvalue_table = TRUE, ES_geom = "line"))
  dev.off()
}

# 10.2 GSEA-KEGG可视化
if (!is.null(gsea_kegg) && nrow(gsea_kegg) > 0) {
  cat("   10.2 GSEA-KEGG 可视化...\n")
  
  # Dotplot
  pdf("GSEA_KEGG_dotplot.pdf", width = 12, height = 8)
  print(dotplot(gsea_kegg, showCategory = 20, font.size = 10))
  dev.off()
  
  # Ridge plot
  pdf("GSEA_KEGG_ridgeplot.pdf", width = 12, height = 8)
  print(ridgeplot(gsea_kegg, showCategory = 20))
  dev.off()
  
  # GSEA plot (展示top 5通路)
  pdf("GSEA_KEGG_enrichment_plot.pdf", width = 14, height = 10)
  print(gseaplot2(gsea_kegg, geneSetID = 1:min(5, nrow(gsea_kegg)), 
                  pvalue_table = TRUE, ES_geom = "line"))
  dev.off()
}

################################################################################
# 11. 比较ORA和GSEA结果 ⭐ 新增
################################################################################
cat("\n11. 比较ORA和GSEA结果...\n")

compare_summary <- data.frame(
  Analysis = c("ORA-GO (up)", "ORA-GO (down)", "ORA-KEGG (up)", "ORA-KEGG (down)",
               "GSEA-GO", "GSEA-KEGG"),
  N_Significant = c(
    ifelse(is.null(go_up), 0, nrow(go_up)),
    ifelse(is.null(go_down), 0, nrow(go_down)),
    ifelse(is.null(kegg_up), 0, nrow(kegg_up)),
    ifelse(is.null(kegg_down), 0, nrow(kegg_down)),
    ifelse(is.null(gsea_go), 0, nrow(gsea_go)),
    ifelse(is.null(gsea_kegg), 0, nrow(gsea_kegg))
  )
)

print(compare_summary)
write.csv(compare_summary, "enrichment_comparison_summary.csv", row.names = FALSE)

################################################################################
# 12. 生成富集分析报告
################################################################################
cat("\n12. 生成富集分析报告...\n")

report <- sprintf("
================================================================================
GO/KEGG/GSEA 富集分析报告 - 小鼠（改进版）
================================================================================
分析日期: %s

1. 差异基因统计:
   - 总基因数: %d
   - 严格阈值 (padj<0.05, |log2FC|>1): %d 个DEGs
   - 宽松阈值 (padj<0.25, |log2FC|>0.3): %d 个DEGs
   - ORA分析使用: %s 阈值 (%d 个DEGs)

2. 基因ID转换:
   - 总基因转换率: %d/%d (%.1f%%)
   - 上调基因: %d
   - 下调基因: %d

3. ORA富集分析结果 (Over-Representation Analysis):
   - 上调基因 GO terms: %d
   - 下调基因 GO terms: %d
   - 上调基因 KEGG通路: %d
   - 下调基因 KEGG通路: %d

4. GSEA富集分析结果 (Gene Set Enrichment Analysis): ⭐ 新增
   - GSEA-GO terms: %d
   - GSEA-KEGG通路: %d
   - GSEA使用全部 %d 个基因（不设阈值）

5. 关键发现:
   %s

================================================================================
输出文件说明:
================================================================================

ID转换:
- gene_id_conversion.csv: ENSEMBL/ENTREZID/SYMBOL对应表

ORA结果文件:
- GO_ORA_upregulated_%s.csv
- GO_ORA_downregulated_%s.csv
- KEGG_ORA_upregulated_%s.csv
- KEGG_ORA_downregulated_%s.csv

GSEA结果文件: ⭐ 新增
- GSEA_GO_results.csv
- GSEA_KEGG_results.csv

ORA可视化:
- GO_ORA_dotplot_*.pdf
- GO_ORA_cnetplot_*.pdf
- KEGG_ORA_dotplot_*.pdf
- KEGG_ORA_cnetplot_*.pdf

GSEA可视化: ⭐ 新增
- GSEA_GO_dotplot.pdf: GO富集气泡图
- GSEA_GO_ridgeplot.pdf: GO富集分数分布
- GSEA_GO_enrichment_plot.pdf: 核心通路详细图
- GSEA_KEGG_dotplot.pdf: KEGG通路气泡图
- GSEA_KEGG_ridgeplot.pdf: KEGG富集分数分布
- GSEA_KEGG_enrichment_plot.pdf: 核心通路详细图

比较分析:
- enrichment_comparison_summary.csv: ORA vs GSEA结果对比

================================================================================
方法说明:
================================================================================
1. ORA (Over-Representation Analysis):
   - 使用阈值筛选的差异基因
   - 检验差异基因在通路中的过表达
   - 适合差异明显的基因集

2. GSEA (Gene Set Enrichment Analysis): ⭐ 推荐
   - 使用全部基因的排序列表（按log2FC排序）
   - 不需要设定阈值
   - 能检测到协同变化但个体变化小的通路
   - **特别适合你的数据（DEG数量较少）**

3. 富集显著性阈值:
   - p.adjust < 0.05 (BH校正)
   - qvalue < 0.2

================================================================================
结果解读建议:
================================================================================
1. **优先看GSEA结果**：
   - 如果ORA结果很少，GSEA往往能找到更多有意义的通路
   - 关注NES (Normalized Enrichment Score)的绝对值和方向
   - NES > 0: 在ZYD组中富集（上调）
   - NES < 0: 在Mod组中富集（下调）

2. 关注与研究相关的通路：
   - 肠道蠕动相关 (gastrointestinal motility)
   - 炎症反应 (inflammatory response)
   - 免疫调控 (immune regulation)
   - 代谢通路 (metabolic pathways)

3. 验证策略：
   - 从GSEA top通路中挑选核心基因 (leading edge genes)
   - 结合ORA结果，找出重复出现的基因
   - 选择3-5个候选基因进行qPCR验证

================================================================================
",
  Sys.Date(),
  nrow(deseq_results),
  nrow(sig_genes_strict),
  nrow(sig_genes_relaxed),
  threshold_label,
  nrow(sig_genes),
  nrow(id_conversion), length(all_ensembl), 
  nrow(id_conversion)/length(all_ensembl)*100,
  nrow(up_entrez),
  nrow(down_entrez),
  ifelse(is.null(go_up), 0, nrow(go_up)),
  ifelse(is.null(go_down), 0, nrow(go_down)),
  ifelse(is.null(kegg_up), 0, nrow(kegg_up)),
  ifelse(is.null(kegg_down), 0, nrow(kegg_down)),
  ifelse(is.null(gsea_go), 0, nrow(gsea_go)),
  ifelse(is.null(gsea_kegg), 0, nrow(gsea_kegg)),
  length(gene_list),
  ifelse(!is.null(gsea_kegg) && nrow(gsea_kegg) > 0,
         sprintf("GSEA检测到 %d 个显著通路，建议重点分析", nrow(gsea_kegg)),
         "建议查看GSEA结果，即使DEG较少，GSEA也能发现协同变化的通路"),
  threshold_label, threshold_label, threshold_label, threshold_label
)

cat(report)
writeLines(report, "enrichment_analysis_report_improved.txt")

cat("\n✓ 富集分析完成！\n")
cat("================================================================================\n")
cat("\n重要提示:\n")
cat("1. ⭐ 优先查看GSEA结果！特别是当ORA结果较少时\n")
cat("2. 查看 GSEA_KEGG_enrichment_plot.pdf 了解核心通路\n")
cat("3. 从GSEA leading edge genes中选择候选基因验证\n")
cat("================================================================================\n")
