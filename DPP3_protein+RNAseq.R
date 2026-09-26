#!/usr/bin/env Rscript

#### Fig S3 #####
# 蛋白组独立分析；对应RNA最终8张补充图，不做联合分析。RStudio可直接Source。
# 输入为已标准化蛋白丰度，不是RNA counts，不使用DESeq2，不重复标准化，不填补缺失。
suppressWarnings(invisible(Sys.setlocale("LC_ALL", "zh_CN.UTF-8")))
required_packages <- c("readxl", "ggplot2", "ggrepel", "limma", "statmod", "matrixStats", "msigdbr", "fgsea", "BiocParallel", "systemfonts", "png")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("请先安装所需R包：", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages(library(ggplot2))

# 输入、统一输出位置及必要参数。以后修改时改为-v2、-v3，PDF/PNG/结果表保持相同后缀。
INPUT_DIR <- "/Volumes/zzData/DPP3/young_old_转录组-蛋白组测序/Results-蛋白组"
OUTPUT_DIR <- "/Volumes/zzData/DPP3/Figure-原始数据整理/Figure 1/DPP3-Figure1-补充图蛋白组"
RNA_DIR <- "/Volumes/zzData/DPP3/Figure-原始数据整理/Figure 1/DPP3-Figure1-补充图转录组"
ABUNDANCE_FILE <- file.path(INPUT_DIR, "标准化蛋白丰度.xlsx")
SOURCE_DEP_FILE <- file.path(INPUT_DIR, "Young_vs_Old.DEP.xls")
BURDEN_DEFINITION_FILE <- file.path(RNA_DIR, "results", "RNA_cartilage_aging_burden_gene_sets.csv")
GSEA_DISPLAY_FILE <- file.path(RNA_DIR, "results", "RNA_GSEA_primary_display.csv")
VERSION <- ""
DPI <- 300L
SEED <- 20260828L
BOOTSTRAP_B <- 5000L
FDR_CUTOFF <- 0.05
LOG2FC_CUTOFF <- 1
stopifnot(dir.exists(OUTPUT_DIR), all(file.exists(c(ABUNDANCE_FILE, SOURCE_DEP_FILE, BURDEN_DEFINITION_FILE, GSEA_DISPLAY_FILE))))
output_file <- function(stem, extension) file.path(OUTPUT_DIR, paste0(stem, VERSION, ".", extension))
write_result <- function(x, stem) {
  f <- output_file(stem, "csv")
  if (file.exists(f)) stop("不覆盖已有结果，请递增VERSION：", f)
  write.csv(x, f, row.names = FALSE, fileEncoding = "UTF-8")
}
input_md5 <- tools::md5sum(c(ABUNDANCE_FILE, SOURCE_DEP_FILE, BURDEN_DEFINITION_FILE, GSEA_DISPLAY_FILE))

# 显式核验Arial，不替换成Helvetica；所有图形文字18pt，绘图区域和图例背景透明。
font_table <- systemfonts::system_fonts()
if (!any(font_table$family == "Arial") || !capabilities("cairo") || !capabilities("png")) stop("Arial或Cairo/PNG图形设备不可用，请先配置；本脚本不替换字体。")
base_family <- "Arial"
font_pt <- 18
text_size <- font_pt / ggplot2::.pt
group_colors <- c(Young = "#3C78A8", Old = "#C44E52")
change_colors <- c(Not_significant = "#BDBDBD", Down = "#3C78A8", Up = "#C44E52", Dpp3 = "#E69F00")
theme_pub <- function() {
  theme_classic(base_size = 18, base_family = "Arial") + theme(
    text = element_text(family = "Arial", size = 18, colour = "black"),
    axis.title = element_text(size = 18, colour = "black"), axis.text = element_text(size = 18, colour = "black"),
    axis.line = element_line(linewidth = 0.5, colour = "black"), axis.ticks = element_line(linewidth = 0.5, colour = "black"),
    legend.title = element_text(size = 18), legend.text = element_text(size = 18), strip.text = element_text(size = 18),
    plot.title = element_blank(), plot.subtitle = element_blank(), plot.caption = element_text(size = 18), plot.tag = element_text(size = 18),
    plot.background = element_rect(fill = "transparent", colour = NA), panel.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA), legend.box.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA), strip.background = element_rect(fill = "transparent", colour = NA),
    plot.margin = margin(12, 12, 10, 10))
}
p_label <- function(p) format.pval(p, digits = 3, eps = 0)
z_score <- function(x) as.numeric(scale(x))

# 显式read_excel读入线性和log2矩阵；DEP.xls实际为制表符文本，使用read.delim。
abundance <- as.data.frame(readxl::read_excel(ABUNDANCE_FILE, sheet = "Abundance"), check.names = FALSE)
log_abundance <- as.data.frame(readxl::read_excel(ABUNDANCE_FILE, sheet = "Log2_abundance"), check.names = FALSE)
normalization <- as.data.frame(readxl::read_excel(ABUNDANCE_FILE, sheet = "Normalization"), check.names = FALSE)
source_dep <- read.delim(SOURCE_DEP_FILE, check.names = FALSE, quote = "", comment.char = "", stringsAsFactors = FALSE)
sample_name <- c(paste0("Young_", 1:6), paste0("Old_", 1:6))
group <- factor(sub("_.*", "", sample_name), levels = c("Young", "Old"))
gene <- trimws(as.character(log_abundance$Protein_Name))
stopifnot(identical(names(abundance)[-1], sample_name), identical(names(log_abundance), names(abundance)),
          identical(abundance$Protein_Name, log_abundance$Protein_Name), !anyNA(gene), all(nzchar(gene)), !anyDuplicated(gene),
          all(table(group) == 6), !any(grepl(";|,|/", gene)))
linear_matrix <- as.matrix(abundance[, sample_name])
log_matrix <- as.matrix(log_abundance[, sample_name])
rownames(linear_matrix) <- rownames(log_matrix) <- gene
stopifnot(identical(is.na(linear_matrix), is.na(log_matrix)), all(linear_matrix[!is.na(linear_matrix)] > 0),
          all(is.finite(log_matrix[!is.na(log_matrix)])))
log_error <- max(abs(log2(linear_matrix) - log_matrix), na.rm = TRUE)
stopifnot(log_error < 1e-9)
observed_young <- rowSums(!is.na(log_matrix[, group == "Young"]))
observed_old <- rowSums(!is.na(log_matrix[, group == "Old"]))
testable <- observed_young >= 3 & observed_old >= 3
complete <- rowSums(is.na(log_matrix)) == 0
variable <- matrixStats::rowSds(log_matrix, na.rm = TRUE) > 0
complete_variable <- complete & variable
sample_qc <- data.frame(sample = sample_name, group = group, observed = colSums(!is.na(log_matrix)),
                       missing = colSums(is.na(log_matrix)), missing_fraction = colMeans(is.na(log_matrix)),
                       median_log2 = apply(log_matrix, 2, median, na.rm = TRUE))

# 差异分析：同一标准化输入、两组独立样本、Old-Young；仅对每组>=3实测值的蛋白检验。
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
fit <- limma::lmFit(log_matrix[testable & variable, , drop = FALSE], design)
fit <- limma::contrasts.fit(fit, limma::makeContrasts(Old - Young, levels = design))
fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
tested <- limma::topTable(fit, number = Inf, sort.by = "none", confint = TRUE)
tested$gene <- rownames(tested)
dep <- data.frame(gene = gene, n_Young_observed = observed_young, n_Old_observed = observed_old,
                  complete_12_samples = complete, Young_mean_log2 = rowMeans(log_matrix[, group == "Young"], na.rm = TRUE),
                  Old_mean_log2 = rowMeans(log_matrix[, group == "Old"], na.rm = TRUE))
hit <- match(gene, tested$gene)
dep$log2FC <- dep$Old_mean_log2 - dep$Young_mean_log2
dep$log2FC[!is.finite(dep$log2FC)] <- NA_real_
dep$FC_Old_over_Young <- 2^dep$log2FC
dep$moderated_t <- tested$t[hit]
dep$P_value <- tested$P.Value[hit]
dep$FDR_BH <- tested$adj.P.Val[hit]
dep$log2FC_CI_low <- tested$CI.L[hit]
dep$log2FC_CI_high <- tested$CI.R[hit]
dep$Regulation <- ifelse(is.na(dep$P_value), "Not_tested", "Not_significant")
dep$Regulation[which(dep$FDR_BH < FDR_CUTOFF & dep$log2FC >= LOG2FC_CUTOFF)] <- "Up"
dep$Regulation[which(dep$FDR_BH < FDR_CUTOFF & dep$log2FC <= -LOG2FC_CUTOFF)] <- "Down"
stopifnot(max(abs(dep$log2FC[!is.na(hit)] - tested$logFC[hit[!is.na(hit)]])) < 1e-9,
          max(abs(p.adjust(tested$P.Value, "BH") - tested$adj.P.Val)) < 1e-12)

# 来源DEP方向核验；重名条目不能确定代表蛋白时，不做强行的一对一匹配。
source_log2fc <- source_dep[["Young_vs_Old.Log2FC"]]
source_fc_error <- max(abs(source_log2fc - log2(source_dep$Old_Mean / source_dep$Young_Mean)), na.rm = TRUE)
stopifnot(source_fc_error < 1e-6)
source_unique <- !duplicated(source_dep$Gene_Name) & !duplicated(source_dep$Gene_Name, fromLast = TRUE)
source_match <- match(dep$gene, source_dep$Gene_Name[source_unique])
source_reference <- source_dep[source_unique, c("Gene_Name", "Young_vs_Old.Log2FC", "Young_vs_Old.p_value", "Young_vs_Old.p_adjust", "Young_vs_Old.Regulation")]
source_comparison <- cbind(dep[, c("gene", "log2FC", "P_value", "FDR_BH", "Regulation")],
                           source_reference[source_match, -1, drop = FALSE])
source_comparison$source_match_status <- ifelse(is.na(source_match), "Not uniquely matched", "Unique Gene_Name; descriptive comparison")

# 固定的10通路检验家族及7个展示词条；不根据蛋白组结果方向或P值删词条。
senmayo <- strsplit("Acvr1b Ang Angpt1 Angptl4 Areg Axl Bex3 Bmp2 Bmp6 C3 Ccl1 Ccl13 Ccl16 Ccl2 Ccl20 Ccl24 Ccl26 Ccl3 Ccl3l1 Ccl4 Ccl5 Ccl7 Ccl8 Cd55 Cd9 Csf1 Csf2 Csf2rb Cst4 Ctnnb1 Ctsb Cxcl1 Cxcl10 Cxcl12 Cxcl16 Cxcl2 Cxcl3 Cxcl8 Cxcr2 Dkk1 Edn1 Egf Egfr Ereg Esm1 Ets2 Fas Fgf1 Fgf2 Fgf7 Gdf15 Gem Gmfg Hgf Hmgb1 Icam1 Icam5 Igf1 Igfbp1 Igfbp2 Igfbp3 Igfbp4 Igfbp5 Igfbp6 Igfbp7 Il10 Il13 Il15 Il18 Il1a Il1b Il2 Il6 Il6st Il7 Inha Iqgap2 Itga2 Itpka Jun Kitl Lcp1 Mif Mmp10 Mmp12 Mmp13 Mmp14 Mmp2 Mmp3 Mmp9 Nap1l4 Nrg1 Pappa Pecam1 Pgf Plat Plau Plaur Ptbp1 Ptger2 Ptges Rps6ka5 Scamp4 Selplg Sema3f Serpinb3c Serpine1 Serpine2 Spp1 Spx Timp2 Tnf Tnfrsf10c Tnfrsf11b Tnfrsf1a Tnfrsf1b Tubgcp2 Vegfa Vegfc Vgf Wnt16 Wnt2", " ", fixed = TRUE)[[1]]
pathway_source <- c("Cellular senescence (GO)" = "GOBP_CELLULAR_SENESCENCE", "SASP (Reactome)" = "REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP",
  "Inflammatory response" = "HALLMARK_INFLAMMATORY_RESPONSE", "TNF-alpha/NF-kB signaling" = "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "IL-6/JAK/STAT3 signaling" = "HALLMARK_IL6_JAK_STAT3_SIGNALING", "Extracellular matrix organization" = "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION",
  "Cartilage development" = "GOBP_CARTILAGE_DEVELOPMENT", "Oxidative phosphorylation" = "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "Mitochondrial translation" = "REACTOME_MITOCHONDRIAL_TRANSLATION")
display_order <- c("IL-6/JAK/STAT3 signaling", "Inflammatory response", "SenMayo senescence / SASP", "Oxidative phosphorylation",
                   "Mitochondrial translation", "Extracellular matrix organization", "Cartilage development")
stopifnot(setequal(display_order, read.csv(GSEA_DISPLAY_FILE)$pathway), length(senmayo) == 122L)
pathway_snapshot_file <- file.path(OUTPUT_DIR, "PROT_fixed_pathway_definitions.csv")
if (file.exists(pathway_snapshot_file)) {
  pathway_definitions <- read.csv(pathway_snapshot_file, stringsAsFactors = FALSE)
  pathways_full <- split(pathway_definitions$gene, pathway_definitions$pathway)
} else {
  if (as.character(packageVersion("msigdbr")) != "7.5.1") stop("初次运行需与RNA相同的msigdbr 7.5.1，或提供已冻结的PROT_fixed_pathway_definitions.csv。")
  msig <- msigdbr::msigdbr(species = "Mus musculus")
  stopifnot(all(pathway_source %in% msig$gs_name))
  pathways_full <- c(list("SenMayo senescence / SASP" = senmayo), lapply(pathway_source, function(id) unique(msig$gene_symbol[msig$gs_name == id])))
  rm(msig)
  pathway_definitions <- data.frame(pathway = rep(names(pathways_full), lengths(pathways_full)), gene = unlist(pathways_full, use.names = FALSE),
                                    definition = rep(c("RNA script fixed SenMayo list", unname(pathway_source)), lengths(pathways_full)), msigdbr_version = "7.5.1")
}
stopifnot(setequal(names(pathways_full), c("SenMayo senescence / SASP", names(pathway_source))))
ranks <- sort(setNames(tested$t, tested$gene), decreasing = TRUE)
stopifnot(all(is.finite(ranks)), !anyDuplicated(names(ranks)), length(ranks) > 1000L)
pathways <- lapply(pathways_full, intersect, y = names(ranks))
stopifnot(all(lengths(pathways) >= 10L), all(lengths(pathways) <= 500L))
set.seed(SEED)
gsea <- as.data.frame(fgsea::fgseaMultilevel(pathways = pathways, stats = ranks, minSize = 10, maxSize = 500, eps = 0, BPPARAM = BiocParallel::SerialParam()))
gsea$leadingEdge <- vapply(gsea$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
gsea$direction <- ifelse(gsea$NES > 0, "Higher in Old", "Higher in Young")
gsea$primary_display <- gsea$pathway %in% display_order
gsea_display <- gsea[match(display_order, gsea$pathway), ]
stopifnot(nrow(gsea_display) == 7L, !anyNA(gsea_display$NES))

# Aging burden定义不变：固定RNA模块成员，Dpp3排除，蛋白Z→模块均值→模块Z→(+---)/4。
# 蛋白组缺失值不补零；仅用12样本均实测的固定成员子集，另存全部成员及未纳入原因。
burden_definition <- read.csv(BURDEN_DEFINITION_FILE, stringsAsFactors = FALSE)
module_ids <- c("Inflammation", "Cartilage_homeostasis", "ECM_homeostasis", "Mitochondrial_homeostasis")
module_sign <- c(1, -1, -1, -1)
module_labels <- c("Inflammation", "Loss of cartilage\ndevelopment", "Loss of ECM\norganization", "Loss of mitochondrial\nprograms")
stopifnot(setequal(unique(burden_definition$module), module_ids), !anyDuplicated(burden_definition[c("module", "gene")]),
          !any(tolower(burden_definition$gene) == "dpp3"))
burden_definition$measured <- burden_definition$gene %in% gene
burden_definition$complete_variable <- burden_definition$gene %in% gene[complete_variable]
burden_definition$included <- burden_definition$complete_variable & tolower(burden_definition$gene) != "dpp3"
burden_definition$exclusion_reason <- ifelse(burden_definition$included, "Included", ifelse(!burden_definition$measured, "Not measured", "Missing sample(s) or zero variance"))
score_gene_sets <- lapply(module_ids, function(m) burden_definition$gene[burden_definition$module == m & burden_definition$included])
names(score_gene_sets) <- module_ids
stopifnot(all(lengths(score_gene_sets) >= 50L), !any(tolower(unlist(score_gene_sets)) == "dpp3"))
module_score <- function(genes) {
  expression_subset <- log_matrix[genes, , drop = FALSE]
  expression_subset <- expression_subset[matrixStats::rowSds(expression_subset) > 0, , drop = FALSE]
  colMeans(t(scale(t(expression_subset))), na.rm = TRUE)
}
module_scores <- vapply(score_gene_sets, module_score, numeric(length(sample_name)))
rownames(module_scores) <- sample_name
module_z <- apply(module_scores, 2, z_score)
rownames(module_z) <- sample_name
oriented_modules <- sweep(module_z, 2, module_sign, "*")
aging_burden <- rowMeans(oriented_modules)
names(aging_burden) <- sample_name
stopifnot(max(abs(aging_burden - (module_z[, 1] - module_z[, 2] - module_z[, 3] - module_z[, 4])/4)) < 1e-12,
          all(is.finite(aging_burden)), abs(mean(aging_burden)) < 1e-12)
module_coverage <- do.call(rbind, lapply(module_ids, function(m) {
  d <- burden_definition[burden_definition$module == m, ]
  data.frame(module = m, RNA_frozen_members = nrow(d), measured = sum(d$measured), scored_complete = sum(d$included),
             scored_fraction_of_definition = mean(d$included), sign = module_sign[match(m, module_ids)])
}))
burden_data <- data.frame(sample = sample_name, group = group,
  setNames(as.data.frame(module_scores), paste0(module_ids, "_mean_gene_z")),
  setNames(as.data.frame(module_z), paste0(module_ids, "_z")),
  setNames(as.data.frame(oriented_modules), paste0(module_ids, "_aging_oriented")), burden = aging_burden, check.names = FALSE)

# 与RNA过程分析相同的Hedges g精确J校正、5000次组内bootstrap及全部924种标签分配。
young_idx <- which(group == "Young"); old_idx <- which(group == "Old")
hedges_df <- length(group) - 2L
hedges_J <- exp(lgamma(hedges_df / 2) - log(hedges_df / 2) / 2 - lgamma((hedges_df - 1) / 2))
set.seed(SEED)
bootstrap_young <- matrix(sample.int(6, 6 * BOOTSTRAP_B, replace = TRUE), nrow = 6)
bootstrap_old <- matrix(sample.int(6, 6 * BOOTSTRAP_B, replace = TRUE), nrow = 6)
old_assignments <- combn(seq_along(group), 6)
effect_statistics <- function(values) {
  y <- values[young_idx]; o <- values[old_idx]
  delta <- mean(o) - mean(y)
  pooled_sd <- sqrt((5 * var(y) + 5 * var(o)) / 10)
  yb <- matrix(y[bootstrap_young], nrow = 6); ob <- matrix(o[bootstrap_old], nrow = 6)
  gb <- hedges_J * (colMeans(ob) - colMeans(yb)) / sqrt((5 * matrixStats::colVars(yb) + 5 * matrixStats::colVars(ob)) / 10)
  finite_g <- is.finite(gb)
  stopifnot(all(is.finite(values)), pooled_sd > 0, sum(finite_g) >= 0.99 * BOOTSTRAP_B)
  ci <- quantile(gb[finite_g], c(0.025, 0.975), names = FALSE)
  perm_old <- colMeans(matrix(values[old_assignments], nrow = 6))
  perm_delta <- perm_old - (sum(values) - 6 * perm_old)/6
  n_extreme <- sum(abs(perm_delta) >= abs(delta) - 1e-12)
  welch <- t.test(o, y, var.equal = FALSE)
  data.frame(n_Young = 6, n_Old = 6, Young_mean = mean(y), Old_mean = mean(o), Old_minus_Young = delta,
    Hedges_g = hedges_J * delta/pooled_sd, g_CI_low = ci[1], g_CI_high = ci[2], bootstrap_valid = sum(finite_g),
    Welch_t = unname(welch$statistic), Welch_df = unname(welch$parameter), Welch_P = welch$p.value,
    mean_difference_CI_low = welch$conf.int[1], mean_difference_CI_high = welch$conf.int[2],
    exact_permutation_P = n_extreme/ncol(old_assignments), extreme_allocations = n_extreme, total_allocations = ncol(old_assignments))
}
components <- cbind(oriented_modules, Composite = aging_burden)
component_effects <- do.call(rbind, lapply(seq_len(ncol(components)), function(i) cbind(component = colnames(components)[i], effect_statistics(components[, i]))))
component_effects$Welch_FDR <- p.adjust(component_effects$Welch_P, "BH")
component_effects$permutation_FDR <- p.adjust(component_effects$exact_permutation_P, "BH")
full_effect <- component_effects[component_effects$component == "Composite", ]
lomo <- vapply(seq_along(module_ids), function(i) rowMeans(oriented_modules[, -i, drop = FALSE]), numeric(12))
colnames(lomo) <- paste0("Without_", module_ids)
lomo_effects <- do.call(rbind, lapply(seq_len(ncol(lomo)), function(i) cbind(definition = colnames(lomo)[i], effect_statistics(lomo[, i]))))
lomo_effects$permutation_FDR <- p.adjust(lomo_effects$exact_permutation_P, "BH")

# Dpp3是预先指定靶点；不筛选其他高相关蛋白来宣称burden有效。年龄组校正df=12-1-2=9。
stopifnot(sum(tolower(gene) == "dpp3") == 1L)
dpp3_data <- data.frame(sample = sample_name, group = group, log2_abundance = as.numeric(log_matrix[tolower(gene) == "dpp3", ]), burden = aging_burden)
stopifnot(!anyNA(dpp3_data))
dpp3_correlations <- do.call(rbind, lapply(c("All", "Young", "Old"), function(gr) {
  take <- if (gr == "All") rep(TRUE, 12) else group == gr
  ct <- cor.test(dpp3_data$log2_abundance[take], aging_burden[take], method = "spearman", exact = FALSE)
  data.frame(analysis = gr, n = sum(take), rho = unname(ct$estimate), P = ct$p.value, df = sum(take) - 2)
}))
rx <- residuals(lm(rank(log2_abundance) ~ group, data = dpp3_data))
ry <- residuals(lm(rank(burden) ~ group, data = dpp3_data))
adjusted_rho <- cor(rx, ry)
partial_df <- nrow(dpp3_data) - 3L
adjusted_P <- 2 * pt(-abs(adjusted_rho * sqrt(partial_df/(1 - adjusted_rho^2))), df = partial_df)
stopifnot(abs(adjusted_P - summary(lm(rank(burden) ~ rank(log2_abundance) + group, data = dpp3_data))$coefficients[2, 4]) < 1e-9)
dpp3_correlations <- rbind(dpp3_correlations, data.frame(analysis = "Age_group_adjusted", n = 12, rho = adjusted_rho, P = adjusted_P, df = partial_df))
dpp3_correlations$FDR_4_tests <- p.adjust(dpp3_correlations$P, "BH")

print(table(dep$Regulation))
print(gsea_display[, c("pathway", "NES", "padj", "size")])
print(module_coverage)
print(component_effects[, c("component", "Hedges_g", "Welch_P", "exact_permutation_P")])
print(dpp3_correlations)

# 先检查全部输出，保证旧文件不覆盖；基因集快照为固定定义，存在时仅只读复用。
figure_stems <- c("Fig_S3a_蛋白组_PCA", "Fig_S3b_差异蛋白_火山图", "Fig_S3c_Top差异蛋白_热图", "Fig_S3d_衰老通路_GSEA",
                  "Fig_S3e_衰老方向统一_模块热图", "Fig_S3f_模块与综合评分_效应量", "Fig_S3g_蛋白组_aging_burden", "Fig_S3h_Dpp3与aging_burden_相关")
figure_sizes <- data.frame(panel = letters[1:8], stem = figure_stems, width_in = c(7.5, 8, 10, 10, 10.8, 11, 7, 8.5),
                          height_in = c(6.3, 6.5, 10, 10, 6.6, 6.6, 6.5, 7.3), dpi = DPI, font = "Arial", font_pt = 18)
figure_files <- as.vector(outer(figure_stems, c("pdf", "png"), Vectorize(output_file)))
result_stems <- c("PROT_DE_results", "PROT_source_DEP_comparison", "PROT_sample_QC", "PROT_log2_abundance", "PROT_input_normalization",
  "PROT_GSEA_results", "PROT_GSEA_primary_display", "PROT_GSEA_ranked_proteins", "PROT_burden_gene_coverage", "PROT_burden_module_coverage",
  "PROT_burden_sample_scores", "PROT_module_effects", "PROT_leave_one_module_out", "PROT_Dpp3_sample_values", "PROT_Dpp3_burden_correlations",
  "PROT_Dpp3_DE_result", "PROT_PCA_scores", "PROT_PCA_proteins", "PROT_sample_correlations", "PROT_top_DEP_heatmap_values",
  "PROT_module_heatmap_values", "PROT_analysis_summary", "PROT_figure_manifest", "PROT_export_QA", "PROT_input_integrity", "PROT_source_DEP_audit")
planned_files <- c(figure_files, output_file(result_stems, "csv"), output_file("PROT_analysis_objects", "rds"),
                   output_file("PROT_图序_方法与结果说明", "md"), output_file("PROT_sessionInfo", "txt"))
if (any(file.exists(planned_files))) stop("发现已有同版本文件；请递增VERSION再运行，勿覆盖旧图或旧结果。")
if (!file.exists(pathway_snapshot_file)) write.csv(pathway_definitions, pathway_snapshot_file, row.names = FALSE, fileEncoding = "UTF-8")
write_result(dep, "PROT_DE_results")
write_result(source_comparison, "PROT_source_DEP_comparison")
write_result(sample_qc, "PROT_sample_QC")
write_result(data.frame(gene = gene, log_matrix, check.names = FALSE), "PROT_log2_abundance")
write_result(normalization, "PROT_input_normalization")
write_result(gsea, "PROT_GSEA_results")
write_result(gsea_display, "PROT_GSEA_primary_display")
write_result(data.frame(gene = names(ranks), moderated_t = unname(ranks)), "PROT_GSEA_ranked_proteins")
write_result(burden_definition, "PROT_burden_gene_coverage")
write_result(module_coverage, "PROT_burden_module_coverage")
write_result(burden_data, "PROT_burden_sample_scores")
write_result(component_effects, "PROT_module_effects")
write_result(lomo_effects, "PROT_leave_one_module_out")
write_result(dpp3_data, "PROT_Dpp3_sample_values")
write_result(dpp3_correlations, "PROT_Dpp3_burden_correlations")
write_result(dep[tolower(dep$gene) == "dpp3", ], "PROT_Dpp3_DE_result")

#### Fig S3a #####
# PCA：选完整实测蛋白中方差最高的500个，不缩放每个蛋白的方差。
pca_genes <- gene[complete_variable][order(matrixStats::rowVars(log_matrix[complete_variable, ]), decreasing = TRUE)][1:500]
pca <- prcomp(t(log_matrix[pca_genes, ]), center = TRUE, scale. = FALSE)
pca_variance <- 100 * pca$sdev^2/sum(pca$sdev^2)
pca_data <- data.frame(sample = sample_name, group = group, pca$x[, 1:3])
write_result(pca_data, "PROT_PCA_scores")
write_result(data.frame(gene = pca_genes, variance = matrixStats::rowVars(log_matrix[pca_genes, ])), "PROT_PCA_proteins")
sample_correlation <- cor(log_matrix[complete_variable, ], method = "pearson")
write_result(data.frame(sample = sample_name, sample_correlation, check.names = FALSE), "PROT_sample_correlations")
p_S3a <- ggplot(pca_data, aes(PC1, PC2, colour = group, shape = group)) +
  stat_ellipse(aes(fill = group), type = "norm", level = 0.95, geom = "polygon", alpha = 0.10, colour = NA, show.legend = FALSE) +
  stat_ellipse(type = "norm", level = 0.95, linewidth = 0.8, show.legend = FALSE) +
  geom_point(size = 4, stroke = 0.6) + scale_colour_manual(values = group_colors) + scale_fill_manual(values = group_colors) +
  scale_shape_manual(values = c(Young = 16, Old = 17)) +
  labs(x = sprintf("PC1 (%.1f%%)", pca_variance[1]), y = sprintf("PC2 (%.1f%%)", pca_variance[2]), colour = NULL, shape = NULL) +
  theme_pub() + theme(legend.position = "top", legend.direction = "horizontal")
width_in <- 7.5; height_in <- 6.3
cairo_pdf(filename = output_file("Fig_S3a_蛋白组_PCA", "pdf"), width = width_in, height = height_in, family = "Arial", pointsize = 18, bg = "transparent")
print(p_S3a)
dev.off()
png(filename = output_file("Fig_S3a_蛋白组_PCA", "png"), width = width_in * DPI, height = height_in * DPI, units = "px", res = DPI, type = "cairo-png", family = "Arial", pointsize = 18, bg = "transparent")
print(p_S3a)
dev.off()

#### Fig S3b #####
# 火山图：主阈值FDR<0.05、|log2FC|>=1；全体受检蛋白和Dpp3均显示。
volcano <- dep[is.finite(dep$P_value), ]
volcano$plot_class <- volcano$Regulation
volcano$plot_class[tolower(volcano$gene) == "dpp3"] <- "Dpp3"
volcano$plot_class <- factor(volcano$plot_class, levels = names(change_colors))
volcano$minus_log10_fdr <- -log10(pmax(volcano$FDR_BH, .Machine$double.xmin))
dpp3_volcano <- volcano[tolower(volcano$gene) == "dpp3", ]
volcano_limit <- max(5, ceiling(max(abs(volcano$log2FC))))
p_S3b <- ggplot(volcano, aes(log2FC, minus_log10_fdr, colour = plot_class)) +
  geom_point(size = 1.6, alpha = 0.65) +
  geom_vline(xintercept = c(-1, 1) * LOG2FC_CUTOFF, linetype = 2, linewidth = 0.45, colour = "#777777") +
  geom_hline(yintercept = -log10(FDR_CUTOFF), linetype = 2, linewidth = 0.45, colour = "#777777") +
  geom_point(data = dpp3_volcano, shape = 21, size = 3.3, fill = "#E69F00", colour = "black", stroke = 0.6) +
  ggrepel::geom_text_repel(data = dpp3_volcano, aes(label = "Dpp3"), colour = "black", family = "Arial", size = text_size,
    fontface = "italic", nudge_y = 0.12 * max(volcano$minus_log10_fdr), min.segment.length = 0, seed = SEED, show.legend = FALSE) +
  scale_colour_manual(values = change_colors, breaks = c("Down", "Up", "Dpp3"), drop = FALSE) +
  scale_x_continuous(limits = c(-volcano_limit, volcano_limit), expand = expansion(mult = 0.02)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 3, alpha = 1))) +
  labs(x = "log2 fold change (Old / Young)", y = "-log10 adjusted P", colour = NULL) + theme_pub() +
  theme(legend.position = "top", legend.direction = "horizontal")
width_in <- 8; height_in <- 6.5
cairo_pdf(filename = output_file("Fig_S3b_差异蛋白_火山图", "pdf"), width = width_in, height = height_in, family = "Arial", pointsize = 18, bg = "transparent")
print(p_S3b)
dev.off()
png(filename = output_file("Fig_S3b_差异蛋白_火山图", "png"), width = width_in * DPI, height = height_in * DPI, units = "px", res = DPI, type = "cairo-png", family = "Arial", pointsize = 18, bg = "transparent")
print(p_S3b)
dev.off()

# 热图共用显示方式：已在外部完成聚类，隐藏树；ggplot保证所有文字使用Arial 18pt。
# 上方色条由实际Young/Old分组生成；色标在±2.5饱和，原始Z值另表保存。
heatmap_plot <- function(mat) {
  nr <- nrow(mat)
  d <- data.frame(x = rep(seq_len(ncol(mat)), each = nr), y = rep(nr:1, times = ncol(mat)), z = as.vector(mat))
  ggplot(d, aes(x, y, fill = z)) + geom_tile(width = 0.98, height = 0.98) +
    annotate("rect", xmin = 0.5, xmax = 6.45, ymin = nr + 0.7, ymax = nr + 1.25, fill = group_colors[["Young"]], colour = NA) +
    annotate("rect", xmin = 6.55, xmax = 12.5, ymin = nr + 0.7, ymax = nr + 1.25, fill = group_colors[["Old"]], colour = NA) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-2.5, 2.5), oob = scales::squish,
                         breaks = c(-2, -1, 0, 1, 2), name = "z") +
    scale_x_continuous(breaks = 1:12, labels = colnames(mat), limits = c(0.5, 12.5), expand = c(0, 0)) +
    scale_y_continuous(breaks = 1:nr, labels = rev(rownames(mat)), limits = c(0.5, nr + 1.45), expand = c(0, 0), position = "right") +
    guides(fill = guide_colourbar(barwidth = grid::unit(4, "mm"), barheight = grid::unit(35, "mm"), title.position = "top")) +
    labs(x = NULL, y = NULL) + theme_pub() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18),
          axis.text.y = element_text(size = 18, margin = margin(l = 6)), legend.position = "right")
}
within_group_order <- function(mat) unlist(lapply(levels(group), function(gr) {
  i <- which(group == gr)
  i[hclust(dist(t(mat[, i, drop = FALSE])), method = "complete")$order]
}), use.names = FALSE)

#### Fig S3c #####
# 差异蛋白热图：完整实测的显著蛋白中按FDR选最多15上调+15下调；不补缺失。
heat_candidates <- dep[dep$Regulation %in% c("Up", "Down") & dep$complete_12_samples, ]
heat_candidates <- heat_candidates[order(heat_candidates$FDR_BH, -abs(heat_candidates$log2FC), heat_candidates$gene), ]
top_dep <- rbind(head(heat_candidates[heat_candidates$Regulation == "Up", ], 15), head(heat_candidates[heat_candidates$Regulation == "Down", ], 15))
stopifnot(nrow(top_dep) >= 2L)
top_z <- t(scale(t(log_matrix[top_dep$gene, , drop = FALSE])))
top_z <- top_z[hclust(dist(top_z), method = "complete")$order, within_group_order(top_z), drop = FALSE]
write_result(data.frame(gene = rownames(top_z), top_z, check.names = FALSE), "PROT_top_DEP_heatmap_values")
p_S3c <- heatmap_plot(top_z)
width_in <- 10; height_in <- 10
cairo_pdf(filename = output_file("Fig_S3c_Top差异蛋白_热图", "pdf"), width = width_in, height = height_in, family = "Arial", pointsize = 18, bg = "transparent")
print(p_S3c)
dev.off()
png(filename = output_file("Fig_S3c_Top差异蛋白_热图", "png"), width = width_in * DPI, height = height_in * DPI, units = "px", res = DPI, type = "cairo-png", family = "Arial", pointsize = 18, bg = "transparent")
print(p_S3c)
dev.off()

#### Fig S3d #####
# GSEA：与RNA相同的7个展示词条和顺序；BH在原10通路家族内校正，反向或不显著也显示。
gsea_display$plot_label <- factor(gsea_display$pathway, levels = rev(display_order))
term_labels <- setNames(display_order, display_order)
term_labels["Extracellular matrix organization"] <- "Extracellular matrix\norganization"
p_S3d <- ggplot(gsea_display, aes(NES, plot_label)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.45, colour = "#777777") +
  geom_point(aes(size = -log10(pmax(padj, .Machine$double.xmin)), colour = NES)) +
  scale_y_discrete(labels = term_labels) +
  scale_colour_gradient2(low = "#3C78A8", mid = "#E7E7E7", high = "#C44E52", midpoint = 0, name = "NES") +
  scale_size_continuous(name = "-log10 FDR", range = c(3.5, 9)) +
  guides(colour = guide_colourbar(barwidth = grid::unit(4, "mm"), barheight = grid::unit(35, "mm"), title.position = "top"),
         size = guide_legend(title.position = "top")) +
  labs(x = "Normalized enrichment score\n(Old / Young)", y = NULL) + theme_pub() +
  theme(legend.position = "right", legend.box = "vertical")
width_in <- 10; height_in <- 10
cairo_pdf(filename = output_file("Fig_S3d_衰老通路_GSEA", "pdf"), width = width_in, height = height_in, family = "Arial", pointsize = 18, bg = "transparent")
print(p_S3d)
dev.off()
png(filename = output_file("Fig_S3d_衰老通路_GSEA", "png"), width = width_in * DPI, height = height_in * DPI, units = "px", res = DPI, type = "cairo-png", family = "Arial", pointsize = 18, bg = "transparent")
print(p_S3d)
dev.off()

#### Fig S3e #####
# 方向统一模块热图：列在各年龄组内聚类，行按固定公式排列，不显示树。
module_order <- within_group_order(t(oriented_modules))
module_heat <- rbind(t(oriented_modules), Composite = z_score(aging_burden))[, module_order]
rownames(module_heat) <- c(module_labels, "Composite burden\n(display z)")
write_result(data.frame(component = rownames(module_heat), module_heat, check.names = FALSE), "PROT_module_heatmap_values")
p_S3e <- heatmap_plot(module_heat)
width_in <- 10.8; height_in <- 6.6
cairo_pdf(filename = output_file("Fig_S3e_衰老方向统一_模块热图", "pdf"), width = width_in, height = height_in, family = "Arial", pointsize = 18, bg = "transparent")
print(p_S3e)
dev.off()
png(filename = output_file("Fig_S3e_衰老方向统一_模块热图", "png"), width = width_in * DPI, height = height_in * DPI, units = "px", res = DPI, type = "cairo-png", family = "Arial", pointsize = 18, bg = "transparent")
print(p_S3e)
dev.off()

#### Fig S3f #####
# 4个模块与综合评分的Hedges g及5000次组内bootstrap百分位95%CI；右侧为原始exact P。
forest_data <- component_effects
forest_labels <- c(module_labels, "Composite\naging burden")
forest_data$label <- factor(forest_labels, levels = rev(forest_labels))
forest_data$kind <- ifelse(forest_data$component == "Composite", "Composite", "Module")
forest_left <- min(0, forest_data$g_CI_low); forest_right <- max(0, forest_data$g_CI_high)
forest_span <- max(1, forest_right - forest_left)
p_S3f <- ggplot(forest_data, aes(Hedges_g, label)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5, colour = "#777777") +
  geom_segment(aes(x = g_CI_low, xend = g_CI_high, yend = label), linewidth = 0.9, colour = "#777777") +
  geom_point(aes(shape = kind), size = 4.4, colour = "#C44E52") +
  geom_text(aes(x = forest_right + 0.10 * forest_span, label = paste0("P = ", p_label(exact_permutation_P))),
            hjust = 0, family = "Arial", size = text_size) +
  scale_shape_manual(values = c(Composite = 18, Module = 16), guide = "none") +
  scale_x_continuous(limits = c(forest_left - 0.05 * forest_span, forest_right + 0.55 * forest_span)) +
  labs(x = "Hedges' g (Old - Young)\n95% bootstrap CI; labels: exact permutation P", y = NULL) +
  theme_pub() + theme(axis.ticks.y = element_blank())
width_in <- 11; height_in <- 6.6
cairo_pdf(filename = output_file("Fig_S3f_模块与综合评分_效应量", "pdf"), width = width_in, height = height_in, family = "Arial", pointsize = 18, bg = "transparent")
print(p_S3f)
dev.off()
png(filename = output_file("Fig_S3f_模块与综合评分_效应量", "png"), width = width_in * DPI, height = height_in * DPI, units = "px", res = DPI, type = "cairo-png", family = "Arial", pointsize = 18, bg = "transparent")
print(p_S3f)
dev.off()

#### Fig S3g #####
# 原定义的蛋白组aging burden，不对最终分数再Z化；所有样本保留。
burden_range <- diff(range(aging_burden))
p_S3g <- ggplot(burden_data, aes(group, burden, colour = group)) +
  geom_boxplot(aes(fill = group), width = 0.5, alpha = 0.12, outlier.shape = NA, linewidth = 0.7, show.legend = FALSE) +
  geom_point(aes(shape = group), position = position_jitter(width = 0.10, height = 0, seed = SEED), size = 4, stroke = 0.6) +
  annotate("text", x = 1.5, y = max(aging_burden) + 0.30 * burden_range,
    label = sprintf("Welch P = %s\nExact permutation P = %s", p_label(full_effect$Welch_P), p_label(full_effect$exact_permutation_P)),
    family = "Arial", size = text_size, lineheight = 1.1) +
  scale_colour_manual(values = group_colors) + scale_fill_manual(values = group_colors) + scale_shape_manual(values = c(Young = 16, Old = 17)) +
  scale_x_discrete(labels = c(Young = "Young\n(n = 6)", Old = "Old\n(n = 6)")) +
  scale_y_continuous(expand = expansion(mult = c(0.07, 0.16))) +
  labs(x = NULL, y = "Relative protein cartilage-aging\nburden (a.u.)") + theme_pub() + theme(legend.position = "none")
width_in <- 7; height_in <- 6.5
cairo_pdf(filename = output_file("Fig_S3g_蛋白组_aging_burden", "pdf"), width = width_in, height = height_in, family = "Arial", pointsize = 18, bg = "transparent")
print(p_S3g)
dev.off()
png(filename = output_file("Fig_S3g_蛋白组_aging_burden", "png"), width = width_in * DPI, height = height_in * DPI, units = "px", res = DPI, type = "cairo-png", family = "Arial", pointsize = 18, bg = "transparent")
print(p_S3g)
dev.off()

#### Fig S3h #####
# Dpp3蛋白与burden相关；灰线/带仅为描述性线性拟合/均值95%CI，标注为秩相关检验。
overall <- dpp3_correlations[dpp3_correlations$analysis == "All", ]
fit_grid <- data.frame(log2_abundance = seq(min(dpp3_data$log2_abundance), max(dpp3_data$log2_abundance), length.out = 100))
fit_ci <- predict(lm(burden ~ log2_abundance, data = dpp3_data), newdata = fit_grid, interval = "confidence")
cor_label_y <- max(aging_burden, fit_ci[, "upr"]) + 0.28 * burden_range
p_S3h <- ggplot(dpp3_data, aes(log2_abundance, burden, colour = group)) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, se = TRUE, colour = "#666666", fill = "#D9D9D9", linewidth = 0.8, show.legend = FALSE) +
  geom_point(aes(shape = group), size = 4, stroke = 0.6) +
  annotate("text", x = mean(range(dpp3_data$log2_abundance)), y = cor_label_y, family = "Arial", size = text_size, lineheight = 1.1,
    label = sprintf("Overall rho = %.2f; P = %s\nAge-group adjusted rho = %.2f; P = %s\nn = 12", overall$rho, p_label(overall$P), adjusted_rho, p_label(adjusted_P))) +
  scale_colour_manual(values = group_colors) + scale_shape_manual(values = c(Young = 16, Old = 17)) +
  scale_x_continuous(expand = expansion(mult = c(0.10, 0.10))) + scale_y_continuous(expand = expansion(mult = c(0.07, 0.15))) +
  labs(x = "Dpp3 log2 normalized protein abundance", y = "Relative protein cartilage-aging\nburden (a.u.)", colour = NULL, shape = NULL) +
  theme_pub() + theme(legend.position = "top", legend.direction = "horizontal")
width_in <- 8.5; height_in <- 7.3
cairo_pdf(filename = output_file("Fig_S3h_Dpp3与aging_burden_相关", "pdf"), width = width_in, height = height_in, family = "Arial", pointsize = 18, bg = "transparent")
print(p_S3h)
dev.off()
png(filename = output_file("Fig_S3h_Dpp3与aging_burden_相关", "png"), width = width_in * DPI, height = height_in * DPI, units = "px", res = DPI, type = "cairo-png", family = "Arial", pointsize = 18, bg = "transparent")
print(p_S3h)
dev.off()

# R内输出审计：Cairo PDF用1pt的Tf结合18倍文本矩阵；按矩阵尺度核验实际字号，不把Tf误读为1pt。
audit_cairo_pdf <- function(path) {
  bytes <- readBin(path, "raw", n = file.info(path)$size)
  ascii <- rawToChar(replace(bytes, bytes == as.raw(0), as.raw(32)))
  starts <- gregexpr("(?<!end)stream\\r?\\n", ascii, perl = TRUE, useBytes = TRUE)[[1]]
  ends <- gregexpr("endstream", ascii, fixed = TRUE, useBytes = TRUE)[[1]]
  lengths_start <- attr(starts, "match.length")
  content <- character()
  for (i in seq_along(starts)) {
    a <- starts[i] + lengths_start[i]; b <- ends[ends > a][1] - 1L
    raw_stream <- tryCatch(memDecompress(bytes[a:b], type = "gzip"), error = function(e) raw())
    if (!length(raw_stream)) next
    s <- rawToChar(replace(raw_stream, raw_stream == as.raw(0), as.raw(32)))
    if (grepl(" Tf", s, fixed = TRUE, useBytes = TRUE) && grepl(" Tm", s, fixed = TRUE, useBytes = TRUE)) content <- c(content, s)
  }
  text <- paste(content, collapse = "\n")
  number <- "[-+]?[0-9]*\\.?[0-9]+"
  tf <- regmatches(text, gregexpr(paste0(number, "[[:space:]]+Tf"), text, perl = TRUE))[[1]]
  tm <- regmatches(text, gregexpr(paste0("(?:", number, "[[:space:]]+){6}Tm"), text, perl = TRUE))[[1]]
  stopifnot(length(tf) > 0, length(tm) > 0, all(as.numeric(sub("[[:space:]]+Tf", "", tf)) == 1))
  sizes <- vapply(strsplit(tm, "[[:space:]]+"), function(x) sqrt(sum(as.numeric(x[1:2])^2)), numeric(1))
  fonts <- unique(regmatches(ascii, gregexpr("/BaseFont /[A-Za-z0-9+_-]+", ascii, perl = TRUE, useBytes = TRUE))[[1]])
  stopifnot(all(abs(sizes - font_pt) < 0.1), length(fonts) > 0, all(grepl("Arial", fonts)), grepl("/FontFile", ascii, fixed = TRUE, useBytes = TRUE))
  data.frame(pdf_min_pt = min(sizes), pdf_max_pt = max(sizes), pdf_fonts = paste(fonts, collapse = ";"))
}
export_qa <- do.call(rbind, lapply(seq_len(8), function(i) {
  z <- png::readPNG(output_file(figure_stems[i], "png"), info = TRUE)
  info <- attr(z, "info")
  stopifnot(length(dim(z)) == 3L, dim(z)[3] == 4L, max(z[c(1, dim(z)[1]), c(1, dim(z)[2]), 4]) == 0,
            abs(dim(z)[2] - figure_sizes$width_in[i] * DPI) <= 1, abs(dim(z)[1] - figure_sizes$height_in[i] * DPI) <= 1,
            all(abs(info$dpi - DPI) < 0.1))
  cbind(data.frame(panel = letters[i], width_px = dim(z)[2], height_px = dim(z)[1], dpi = info$dpi[1],
                   transparent_fraction = mean(z[, , 4] == 0), png_corners_transparent = TRUE), audit_cairo_pdf(output_file(figure_stems[i], "pdf")))
}))
write_result(export_qa, "PROT_export_QA")
write_result(figure_sizes, "PROT_figure_manifest")
stopifnot(identical(input_md5, tools::md5sum(names(input_md5))), all(file.exists(figure_files)))
write_result(data.frame(input = names(input_md5), md5_before = unname(input_md5), md5_after = unname(tools::md5sum(names(input_md5))), unchanged = TRUE), "PROT_input_integrity")
source_expected_label <- rep("Nodiff", nrow(source_dep))
source_expected_label[which(source_dep[["Young_vs_Old.p_value"]] < 0.05 & source_log2fc >= 1)] <- "Up"
source_expected_label[which(source_dep[["Young_vs_Old.p_value"]] < 0.05 & source_log2fc <= -1)] <- "Down"
source_audit <- data.frame(metric = c("source_rows", "source_duplicate_names", "FC_is_Old_over_Young_max_log2_error", "source_labels_match_rawP_FC_rule_fraction", "normalized_log2_crosscheck_error"),
                           value = c(nrow(source_dep), sum(duplicated(source_dep$Gene_Name)), source_fc_error,
                                     mean(source_expected_label == source_dep[["Young_vs_Old.Regulation"]]), log_error))
write_result(source_audit, "PROT_source_DEP_audit")
analysis_summary <- data.frame(metric = c("samples", "Young", "Old", "input_proteins", "complete_proteins", "proteins_tested", "not_tested", "up_FDR_abslog2FC1", "down_FDR_abslog2FC1",
  "PCA_PC1_percent", "PCA_PC2_percent", "burden_mean_difference", "burden_Hedges_g", "burden_g_CI_low", "burden_g_CI_high", "burden_Welch_P", "burden_exact_P", "Dpp3_overall_rho", "Dpp3_adjusted_rho", "Dpp3_adjusted_P"),
  value = c(12, 6, 6, nrow(log_matrix), sum(complete), nrow(tested), sum(is.na(dep$P_value)), sum(dep$Regulation == "Up"), sum(dep$Regulation == "Down"),
            pca_variance[1:2], full_effect$Old_minus_Young, full_effect$Hedges_g, full_effect$g_CI_low, full_effect$g_CI_high,
            full_effect$Welch_P, full_effect$exact_permutation_P, overall$rho, adjusted_rho, adjusted_P))
write_result(analysis_summary, "PROT_analysis_summary")
saveRDS(list(dep = dep, gsea = gsea, gsea_display = gsea_display, module_scores = module_scores, module_z = module_z,
             burden = burden_data, effects = component_effects, lomo = lomo_effects, Dpp3 = dpp3_data, Dpp3_correlations = dpp3_correlations,
             plots = list(a = p_S3a, b = p_S3b, c = p_S3c, d = p_S3d, e = p_S3e, f = p_S3f, g = p_S3g, h = p_S3h)),
        output_file("PROT_analysis_objects", "rds"))
writeLines(capture.output(sessionInfo()), output_file("PROT_sessionInfo", "txt"), useBytes = TRUE)

# 图注与方法直接随结果导出，区分事实、假设和未完成验证。
readme <- c("# 蛋白组补充分析：图序、方法与结果", "",
  "本轮仅蛋白组分析；未做RNA–蛋白联合分析、样本配对相关或绝对评分比较。所有图为独立PDF和透明底300dpi PNG，Arial 18pt。", "",
  "## 图序", "",
  "- S3a：500个最高方差、12样本均实测的蛋白做PCA，输入log2标准化丰度，中心化、不缩放；95%正态理论数据椭圆用于描述样本分布，不是组均值置信区间或显著性检验。",
  "- S3b：limma差异蛋白火山图，正方向=Old更高；所有受检蛋白均绘制。主阈值BH FDR<0.05且|log2FC|≥1，Dpp3单独标示，不因其为靶点改变判定阈值。",
  "- S3c：显著且12样本均实测的蛋白中，按FDR选择最多15个上调和15个下调；蛋白行Z。行及各年龄组内列按欧氏距离/complete linkage聚类，隐藏树。此图为差异结果的描述性展示，不是独立验证。",
  "- S3d：原RNA最终7个GSEA词条全部显示，顺序固定；正NES=Old高。点大小=-log10 FDR、颜色=NES。FDR在原10通路家族内BH校正；没有按预期方向或显著性筛选展示。",
  "- S3e：四个模块均转换为越高越aging-oriented的方向；列在组内欧氏距离/complete linkage聚类，隐藏树；综合行仅显示时做Z。热图色标在±2.5饱和，原值另存。",
  "- S3f：模块及综合评分的Hedges g（Old−Young），精确gamma小样本校正；线为5000次组内bootstrap百分位95%CI；标签是双侧精确置换原始P。",
  "- S3g：队列内相对protein cartilage-aging burden，箱体IQR、中位数线、1.5×IQR须；全部个体绘制。显示双侧Welch及精确置换P；综合分数不再Z化。",
  "- S3h：Dpp3 log2蛋白丰度与burden的总体Spearman相关，以及年龄组校正偏秩相关；灰线和灰带是描述性线性拟合/均值95%CI，不是Spearman置信带。", "",
  "## 输入和差异模型", "",
  "读取标准化Excel的Abundance、Log2_abundance和Normalization。两张丰度表逐值核对；已有代表蛋白选择和标准化不重复进行。保留NA，不补零或随机低值。",
  "按Young和Old各6个独立样本处理。每组至少3个有效观测且非零方差的蛋白用limma::lmFit、Old−Young对比及eBayes(trend=TRUE,robust=TRUE)分析；BH覆盖全部受检蛋白。未检验蛋白保留在总表并明确标记。",
  "原始DEP实际FC为Old_Mean/Young_Mean，虽然列名写Young_vs_Old；其样本丰度已补全、来源填补方法未知。来源Up/Down标签基于原始P，不当作本轮FDR结论。新主分析的FC/P/FDR均来自同一个标准化矩阵和模型，没有拼接来源统计值。重名来源仅作审计，不强行匹配代表蛋白。", "",
  "## 固定GSEA与aging burden", "",
  "GSEA按全部受检蛋白的有符号moderated t预排序，fgseaMultilevel(minSize=10,maxSize=500,eps=0)，串行执行、固定随机种子；背景不是差异蛋白子集。9个MSigDB通路沿用msigdbr 7.5.1的Mus musculus映射，另加原R脚本122基因的固定SenMayo名单；已导出冻结定义，后续复现优先读取快照。",
  "aging burden严格保持原公式：[Z(Inflammation)−Z(Cartilage development)−Z(ECM organization)−Z(Mitochondrial programs)]/4。每蛋白先在12样本中Z标准化，再在模块内求均值，各模块再次在12样本中Z标准化；四模块等权，Dpp3排除；不加入SenMayo、不按组差异选基因或调整符号。",
  "模块成员取原RNA冻结定义与本次蛋白组完整实测、非零方差蛋白的交集。固定成员可避免不同样本用不同蛋白计算模块；代价是评分仅覆盖可测子集，缺失非随机及覆盖不足仍会影响结果。蛋白组评分不是RNA评分的数值替代，不能当作绝对生物学年龄或衰老细胞比例。完整成员、纳入原因及覆盖率见PROT_burden_gene_coverage和PROT_burden_module_coverage。",
  paste(capture.output(print(module_coverage, row.names = FALSE)), collapse = "\n"), "",
  "## 评分推断与限制", "",
  "bootstrap在各组6个固定分数内有放回抽样5000次；不重估上游标准化或两层Z。精确置换枚举924种6:6分配，统计量为绝对组均值差；条件于固定分数及可交换性假设，不能消除批次、细胞组成或缺失机制的影响。",
  "Welch与置换各在4模块+综合的5项家族内BH校正；留一模块4项另成家族。图中为原始P，完整FDR另表提供。bootstrap区间与置换P是不同推断方法，小样本下可能不一致。",
  "相关P为双侧近似检验。偏秩相关先分别将两个秩变量对年龄组回归，再相关残差；t近似df=9，已与秩回归P交叉核验。总体、Young、Old和年龄组校正共4项另做BH。Dpp3虽不参与评分，关联仍不能证明因果或独立外部有效性。",
  "AUTHOR_INPUT_NEEDED：确认每列表征独立小鼠还是混样/技术重复，以及是否存在配对、批次或窝别。鉴定FDR、肽段数和缺失机制未提供；目前不能据此使用肽段精度加权或评估全部实验层面质量。", "",
  "## 本次结果", "",
  sprintf("共%d个输入蛋白，%d个受检，%d个未检验；FDR<0.05且|log2FC|≥1：Old上调%d，下调%d。", nrow(dep), nrow(tested), sum(is.na(dep$P_value)), sum(dep$Regulation == "Up"), sum(dep$Regulation == "Down")),
  sprintf("综合burden Old−Young=%.4f，Hedges g=%.3f（95%%CI %.3f–%.3f），Welch P=%s，exact P=%s。", full_effect$Old_minus_Young, full_effect$Hedges_g,
          full_effect$g_CI_low, full_effect$g_CI_high, p_label(full_effect$Welch_P), p_label(full_effect$exact_permutation_P)),
  paste(capture.output(print(gsea_display[, c("pathway", "NES", "padj", "size")], row.names = FALSE)), collapse = "\n"),
  sprintf("Dpp3总体rho=%.3f，P=%s；年龄组校正rho=%.3f，P=%s。", overall$rho, p_label(overall$P), adjusted_rho, p_label(adjusted_P)),
  "SenMayo/SASP正向但FDR未达0.05；线粒体翻译GSEA未达显著；ECM模块单独组间差异也未达显著。均如实保留，不能把所有衰老条目写成显著升高。", "",
  "## 复现和输出", "",
  "直接Source本目录完整R脚本；输入路径、依赖、阈值、种子在开头。每图下方都有可单独执行的cairo_pdf和png保存代码。重跑时递增VERSION（例如-v2），不新建文件夹、不覆盖旧文件。",
  "所有关键计算均在交付的R文件内。PROT_export_QA由R检查PNG透明通道、尺寸、300dpi和Cairo PDF字体/18pt尺度；未以白底图片替代透明输出。R和包版本见sessionInfo。",
  "旧转录组脚本含有既存语法问题，本轮仅查阅定义、读取冻结成员/显示词条，不执行或修改旧脚本；新蛋白组脚本不依赖旧R会话对象。", "",
  "## 方法来源", "",
  "- limma方法与缺失值说明：https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf ; https://support.bioconductor.org/p/9142613/",
  "- fgsea：https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html",
  "- SenMayo原始研究：https://www.nature.com/articles/s41467-022-32552-1"
)
writeLines(readme, output_file("PROT_图序_方法与结果说明", "md"), useBytes = TRUE)
message("蛋白组分析完成；全部输出直接保存至：", OUTPUT_DIR)
