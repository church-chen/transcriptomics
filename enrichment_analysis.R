#!/usr/bin/env Rscript

################################################################################
# GO/KEGG 富集分析脚本 - 小鼠
# 输入：DESeq2差异分析结果
# 输出：GO/KEGG富集分析结果和可视化图表
################################################################################

cat("================================================================================\n")
cat("GO/KEGG 富集分析 - 小鼠\n")
cat("================================================================================\n\n")

# 加载必需的R包
cat("加载R包...\n")
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)  # 小鼠注释数据库
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
  library(stringr)
})

################################################################################
# 1. 读取差异基因数据
################################################################################
cat("\n1. 读取差异基因数据...\n")

# 读取差异基因结果
sig_genes <- read.csv("significant_genes_relaxed.csv", row.names = 1)
cat(sprintf("   总差异基因数: %d\n", nrow(sig_genes)))

# 分离上调和下调基因
up_genes <- sig_genes[sig_genes$log2FoldChange > 0, ]
down_genes <- sig_genes[sig_genes$log2FoldChange < 0, ]

cat(sprintf("   上调基因: %d\n", nrow(up_genes)))
cat(sprintf("   下调基因: %d\n", nrow(down_genes)))

################################################################################
# 2. 基因ID转换（ENSEMBL -> ENTREZID）
################################################################################
cat("\n2. 基因ID转换...\n")

# 提取ENSEMBL ID（去除版本号，如.1, .2等）
up_ensembl <- gsub("\\..*", "", rownames(up_genes))
down_ensembl <- gsub("\\..*", "", rownames(down_genes))
all_ensembl <- gsub("\\..*", "", rownames(sig_genes))

cat(sprintf("   上调基因 ENSEMBL IDs: %d\n", length(up_ensembl)))
cat(sprintf("   下调基因 ENSEMBL IDs: %d\n", length(down_ensembl)))

# ENSEMBL转ENTREZID（用于KEGG分析）
up_entrez <- bitr(up_ensembl,
                  fromType = "ENSEMBL",
                  toType = "ENTREZID",
                  OrgDb = org.Mm.eg.db)

down_entrez <- bitr(down_ensembl,
                    fromType = "ENSEMBL",
                    toType = "ENTREZID",
                    OrgDb = org.Mm.eg.db)

all_entrez <- bitr(all_ensembl,
                   fromType = "ENSEMBL",
                   toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)

cat(sprintf("   上调基因转换成功: %d/%d (%.1f%%)\n",
            nrow(up_entrez), length(up_ensembl),
            nrow(up_entrez)/length(up_ensembl)*100))
cat(sprintf("   下调基因转换成功: %d/%d (%.1f%%)\n",
            nrow(down_entrez), length(down_ensembl),
            nrow(down_entrez)/length(down_ensembl)*100))

################################################################################
# 3. GO富集分析
################################################################################
cat("\n3. GO富集分析...\n")

# 3.1 上调基因 GO 富集
cat("   3.1 上调基因 GO 富集分析...\n")
go_up <- enrichGO(gene = up_ensembl,
                  OrgDb = org.Mm.eg.db,
                  keyType = "ENSEMBL",
                  ont = "ALL",  # BP, MF, CC, ALL
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)

cat(sprintf("     富集到的GO term: %d\n", nrow(go_up)))

# 3.2 下调基因 GO 富集
cat("   3.2 下调基因 GO 富集分析...\n")
go_down <- enrichGO(gene = down_ensembl,
                    OrgDb = org.Mm.eg.db,
                    keyType = "ENSEMBL",
                    ont = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2,
                    readable = TRUE)

cat(sprintf("     富集到的GO term: %d\n", nrow(go_down)))

# 3.3 保存GO结果
if (!is.null(go_up) && nrow(go_up) > 0) {
  write.csv(as.data.frame(go_up), "GO_enrichment_upregulated.csv", row.names = FALSE)
}

if (!is.null(go_down) && nrow(go_down) > 0) {
  write.csv(as.data.frame(go_down), "GO_enrichment_downregulated.csv", row.names = FALSE)
}

################################################################################
# 4. KEGG通路富集分析
################################################################################
cat("\n4. KEGG通路富集分析...\n")

# 4.1 上调基因 KEGG 富集
cat("   4.1 上调基因 KEGG 富集分析...\n")
kegg_up <- enrichKEGG(gene = up_entrez$ENTREZID,
                      organism = "mmu",  # 小鼠
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2)

if (!is.null(kegg_up) && nrow(kegg_up) > 0) {
  # 添加基因symbol
  kegg_up <- setReadable(kegg_up, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  cat(sprintf("     富集到的KEGG通路: %d\n", nrow(kegg_up)))
} else {
  cat("     未找到显著富集的KEGG通路\n")
}

# 4.2 下调基因 KEGG 富集
cat("   4.2 下调基因 KEGG 富集分析...\n")
kegg_down <- enrichKEGG(gene = down_entrez$ENTREZID,
                        organism = "mmu",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2)

if (!is.null(kegg_down) && nrow(kegg_down) > 0) {
  kegg_down <- setReadable(kegg_down, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  cat(sprintf("     富集到的KEGG通路: %d\n", nrow(kegg_down)))
} else {
  cat("     未找到显著富集的KEGG通路\n")
}

# 4.3 保存KEGG结果
if (!is.null(kegg_up) && nrow(kegg_up) > 0) {
  write.csv(as.data.frame(kegg_up), "KEGG_enrichment_upregulated.csv", row.names = FALSE)
}

if (!is.null(kegg_down) && nrow(kegg_down) > 0) {
  write.csv(as.data.frame(kegg_down), "KEGG_enrichment_downregulated.csv", row.names = FALSE)
}

################################################################################
# 5. 可视化 - GO富集分析
################################################################################
cat("\n5. 生成GO富集分析可视化图表...\n")

# 5.1 上调基因 GO 可视化
if (!is.null(go_up) && nrow(go_up) > 0) {
  cat("   5.1 上调基因 GO 可视化...\n")

  # 柱状图
  pdf("GO_barplot_upregulated.pdf", width = 12, height = 10)
  p1 <- barplot(go_up,
                showCategory = 20,
                split = "ONTOLOGY",
                font.size = 10,
                title = "GO Enrichment - Upregulated Genes") +
    facet_grid(ONTOLOGY~., scales = "free")
  print(p1)
  dev.off()

  # 气泡图
  pdf("GO_dotplot_upregulated.pdf", width = 12, height = 10)
  p2 <- dotplot(go_up,
                showCategory = 10,
                split = "ONTOLOGY",
                font.size = 10,
                title = "GO Enrichment - Upregulated Genes") +
    facet_grid(ONTOLOGY~., scales = "free")
  print(p2)
  dev.off()

  # 基因-概念网络图（选择top10）
  if (nrow(go_up) >= 5) {
    pdf("GO_cnetplot_upregulated.pdf", width = 14, height = 12)
    p3 <- cnetplot(go_up,
                   showCategory = 5,
                   colorEdge = TRUE,
                   circular = FALSE,
                   node_label = "all")
    print(p3)
    dev.off()
  }
}

# 5.2 下调基因 GO 可视化
if (!is.null(go_down) && nrow(go_down) > 0) {
  cat("   5.2 下调基因 GO 可视化...\n")

  # 柱状图
  pdf("GO_barplot_downregulated.pdf", width = 12, height = 10)
  p4 <- barplot(go_down,
                showCategory = 20,
                split = "ONTOLOGY",
                font.size = 10,
                title = "GO Enrichment - Downregulated Genes") +
    facet_grid(ONTOLOGY~., scales = "free")
  print(p4)
  dev.off()

  # 气泡图
  pdf("GO_dotplot_downregulated.pdf", width = 12, height = 10)
  p5 <- dotplot(go_down,
                showCategory = 10,
                split = "ONTOLOGY",
                font.size = 10,
                title = "GO Enrichment - Downregulated Genes") +
    facet_grid(ONTOLOGY~., scales = "free")
  print(p5)
  dev.off()

  # 基因-概念网络图
  if (nrow(go_down) >= 5) {
    pdf("GO_cnetplot_downregulated.pdf", width = 14, height = 12)
    p6 <- cnetplot(go_down,
                   showCategory = 5,
                   colorEdge = TRUE,
                   circular = FALSE,
                   node_label = "all")
    print(p6)
    dev.off()
  }
}

################################################################################
# 6. 可视化 - KEGG通路富集分析
################################################################################
cat("\n6. 生成KEGG富集分析可视化图表...\n")

# 6.1 上调基因 KEGG 可视化
if (!is.null(kegg_up) && nrow(kegg_up) > 0) {
  cat("   6.1 上调基因 KEGG 可视化...\n")

  # 柱状图
  pdf("KEGG_barplot_upregulated.pdf", width = 12, height = 8)
  p7 <- barplot(kegg_up,
                showCategory = 20,
                font.size = 10,
                title = "KEGG Pathway Enrichment - Upregulated Genes")
  print(p7)
  dev.off()

  # 气泡图
  pdf("KEGG_dotplot_upregulated.pdf", width = 12, height = 8)
  p8 <- dotplot(kegg_up,
                showCategory = 20,
                font.size = 10,
                title = "KEGG Pathway Enrichment - Upregulated Genes")
  print(p8)
  dev.off()

  # 基因-概念网络图
  if (nrow(kegg_up) >= 5) {
    pdf("KEGG_cnetplot_upregulated.pdf", width = 14, height = 12)
    p9 <- cnetplot(kegg_up,
                   showCategory = 5,
                   colorEdge = TRUE,
                   circular = FALSE,
                   node_label = "all")
    print(p9)
    dev.off()
  }
}

# 6.2 下调基因 KEGG 可视化
if (!is.null(kegg_down) && nrow(kegg_down) > 0) {
  cat("   6.2 下调基因 KEGG 可视化...\n")

  # 柱状图
  pdf("KEGG_barplot_downregulated.pdf", width = 12, height = 8)
  p10 <- barplot(kegg_down,
                 showCategory = 20,
                 font.size = 10,
                 title = "KEGG Pathway Enrichment - Downregulated Genes")
  print(p10)
  dev.off()

  # 气泡图
  pdf("KEGG_dotplot_downregulated.pdf", width = 12, height = 8)
  p11 <- dotplot(kegg_down,
                 showCategory = 20,
                 font.size = 10,
                 title = "KEGG Pathway Enrichment - Downregulated Genes")
  print(p11)
  dev.off()

  # 基因-概念网络图
  if (nrow(kegg_down) >= 5) {
    pdf("KEGG_cnetplot_downregulated.pdf", width = 14, height = 12)
    p12 <- cnetplot(kegg_down,
                    showCategory = 5,
                    colorEdge = TRUE,
                    circular = FALSE,
                    node_label = "all")
    print(p12)
    dev.off()
  }
}

################################################################################
# 7. 生成富集分析报告
################################################################################
cat("\n7. 生成富集分析报告...\n")

report <- sprintf("
================================================================================
GO/KEGG 富集分析报告 - 小鼠
================================================================================
分析日期: %s

差异基因统计:
- 总差异基因: %d
- 上调基因: %d
- 下调基因: %d

基因ID转换:
- 上调基因 ENSEMBL->ENTREZ: %d/%d (%.1f%%)
- 下调基因 ENSEMBL->ENTREZ: %d/%d (%.1f%%)

GO富集分析结果:
- 上调基因富集GO term: %d
- 下调基因富集GO term: %d

KEGG通路富集分析结果:
- 上调基因富集KEGG通路: %d
- 下调基因富集KEGG通路: %d

输出文件:
================================================================================
CSV结果文件:
1. GO_enrichment_upregulated.csv - 上调基因GO富集结果
2. GO_enrichment_downregulated.csv - 下调基因GO富集结果
3. KEGG_enrichment_upregulated.csv - 上调基因KEGG富集结果
4. KEGG_enrichment_downregulated.csv - 下调基因KEGG富集结果

GO可视化文件:
5. GO_barplot_upregulated.pdf - 上调基因GO柱状图
6. GO_dotplot_upregulated.pdf - 上调基因GO气泡图
7. GO_cnetplot_upregulated.pdf - 上调基因GO网络图
8. GO_barplot_downregulated.pdf - 下调基因GO柱状图
9. GO_dotplot_downregulated.pdf - 下调基因GO气泡图
10. GO_cnetplot_downregulated.pdf - 下调基因GO网络图

KEGG可视化文件:
11. KEGG_barplot_upregulated.pdf - 上调基因KEGG柱状图
12. KEGG_dotplot_upregulated.pdf - 上调基因KEGG气泡图
13. KEGG_cnetplot_upregulated.pdf - 上调基因KEGG网络图
14. KEGG_barplot_downregulated.pdf - 下调基因KEGG柱状图
15. KEGG_dotplot_downregulated.pdf - 下调基因KEGG气泡图
16. KEGG_cnetplot_downregulated.pdf - 下调基因KEGG网络图

================================================================================
注意事项:
- 富集分析使用 p.adjust < 0.05, qvalue < 0.2 作为显著性阈值
- GO分析包含三个本体：BP (生物学过程), MF (分子功能), CC (细胞组分)
- KEGG通路分析使用小鼠数据库 (mmu)
- 基因symbol已添加到结果中，方便查看

建议:
1. 重点关注与便秘、肠道蠕动、炎症相关的GO term和KEGG通路
2. 查看network图，找出核心调控基因 (hub genes)
3. 结合文献验证关键通路的生物学意义

================================================================================
",
  Sys.Date(),
  nrow(sig_genes),
  nrow(up_genes),
  nrow(down_genes),
  nrow(up_entrez), length(up_ensembl), nrow(up_entrez)/length(up_ensembl)*100,
  nrow(down_entrez), length(down_ensembl), nrow(down_entrez)/length(down_ensembl)*100,
  ifelse(is.null(go_up), 0, nrow(go_up)),
  ifelse(is.null(go_down), 0, nrow(go_down)),
  ifelse(is.null(kegg_up), 0, nrow(kegg_up)),
  ifelse(is.null(kegg_down), 0, nrow(kegg_down))
)

cat(report)
writeLines(report, "enrichment_analysis_report.txt")

cat("\n✓ 富集分析完成！\n")
cat("================================================================================\n")
