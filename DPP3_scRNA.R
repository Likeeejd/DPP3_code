#!/usr/bin/env Rscript

suppressWarnings(invisible(Sys.setlocale("LC_ALL", "zh_CN.UTF-8")))

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(ggrepel)
  library(edgeR)
  library(limma)
  library(UCell)
  library(fgsea)
  library(babelgene)
  library(ragg)
  library(svglite)
})

#### 运行设置：RStudio 中可直接点击 Source 运行 -------------------------------

input_file <- Sys.getenv(
  "DPP3_SCRNA_INPUT",
  unset = "/Volumes/zzData/scRNA_实验_JAR/单细胞学习/liuchuanju_scRNA_data_202474.rdata"
)
mouse_signature_file <- Sys.getenv(
  "DPP3_MOUSE_AGING_SIGNATURE",
  unset = "/Volumes/zzData/DPP3/Figure-原始数据整理/Figure 1/DPP3-Figure1-补充图/results/RNA_cartilage_aging_burden_gene_sets.csv"
)
output_root <- Sys.getenv(
  "DPP3_SCRNA_OUTPUT",
  unset = "/Volumes/zzData/DPP3/Figure-原始数据整理/Figure 1/DPP3-单细胞完整分析"
)

figure_main_dir <- file.path(output_root, "Figures", "Main")
figure_supp_dir <- file.path(output_root, "Figures", "Supplementary")
result_dir <- file.path(output_root, "results")
source_dir <- file.path(result_dir, "Source_Data")
de_dir <- file.path(result_dir, "Pseudobulk_DE")
for (path in c(figure_main_dir, figure_supp_dir, result_dir, source_dir, de_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

base_family <- "Arial"
condition_colors <- c(Normal = "#3C78A8", OA = "#C44E52")
celltype_levels <- c("EC", "FC", "HomC", "ProC", "RegC", "preHTC", "HTC")
celltype_colors <- c(
  EC = "#009E73", FC = "#7F7F7F", HomC = "#0072B2", ProC = "#E69F00",
  RegC = "#8C6BB1", preHTC = "#D55E00", HTC = "#CC79A7"
)
minimum_cells_per_donor_subtype <- 20L

theme_pub <- function(base_size = 14) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.title = element_text(size = base_size, colour = "black"),
      axis.text = element_text(size = base_size - 1, colour = "black"),
      axis.line = element_line(linewidth = 0.45, colour = "black"),
      axis.ticks = element_line(linewidth = 0.45, colour = "black"),
      legend.title = element_text(size = base_size - 1),
      legend.text = element_text(size = base_size - 2),
      strip.text = element_text(size = base_size - 1, face = "bold", colour = "black"),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

save_plot <- function(plot, stem, directory, width_mm = 180, height_mm = 135) {
  pdf_file <- file.path(directory, paste0(stem, ".pdf"))
  svg_file <- file.path(directory, paste0(stem, ".svg"))
  png_file <- file.path(directory, paste0(stem, ".png"))
  png_temp <- tempfile(pattern = "ragg_", tmpdir = directory, fileext = ".png")

  grDevices::cairo_pdf(
    pdf_file, width = width_mm / 25.4, height = height_mm / 25.4,
    family = base_family, bg = "white"
  )
  print(plot)
  grDevices::dev.off()

  svglite::svglite(
    svg_file, width = width_mm / 25.4, height = height_mm / 25.4,
    bg = "white"
  )
  print(plot)
  grDevices::dev.off()

  ragg::agg_png(
    png_temp, width = width_mm, height = height_mm, units = "mm",
    res = 600, background = "white"
  )
  print(plot)
  grDevices::dev.off()
  stopifnot(file.rename(png_temp, png_file))
}

write_csv_gz <- function(x, filename, row.names = FALSE) {
  connection <- gzfile(filename, open = "wt")
  on.exit(close(connection), add = TRUE)
  write.csv(x, connection, row.names = row.names)
}

z_safe <- function(x) {
  if (sum(is.finite(x)) < 2L || stats::sd(x, na.rm = TRUE) == 0) {
    return(rep(0, length(x)))
  }
  as.numeric(scale(x))
}

group_test <- function(data, value, subtype_name) {
  x <- data[[value]][data$condition == "Normal"]
  y <- data[[value]][data$condition == "OA"]
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) < 2L || length(y) < 2L) {
    return(data.frame(
      Celltype = subtype_name, n_Normal = length(x), n_OA = length(y),
      mean_Normal = mean(x), mean_OA = mean(y), OA_minus_Normal = mean(y) - mean(x),
      ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_
    ))
  }
  test <- t.test(y, x)
  data.frame(
    Celltype = subtype_name, n_Normal = length(x), n_OA = length(y),
    mean_Normal = mean(x), mean_OA = mean(y), OA_minus_Normal = mean(y) - mean(x),
    ci_low = unname(test$conf.int[1]), ci_high = unname(test$conf.int[2]),
    p_value = test$p.value
  )
}

#### 1. 读入已清洗对象并锁定统计单位 -----------------------------------------

stopifnot(file.exists(input_file), file.exists(mouse_signature_file))
loaded_objects <- load(input_file)
stopifnot("srsc" %in% loaded_objects, inherits(srsc, "Seurat"))
stopifnot(
  identical(DefaultAssay(srsc), "RNA"),
  all(c("counts", "data") %in% slotNames(srsc[["RNA"]])),
  all(c("orig.ident", "group", "Celltype", "nCount_RNA", "nFeature_RNA", "percent.mt") %in%
        colnames(srsc@meta.data)),
  all(c("pca", "harmony", "umap", "tsne") %in% Reductions(srsc)),
  "DPP3" %in% rownames(srsc)
)

metadata <- srsc@meta.data
metadata$cell_id <- rownames(metadata)
metadata$donor <- sub("_raw$", "", sub("^GSM[0-9]+_", "", metadata$orig.ident))
metadata$condition <- factor(metadata$group, levels = c("Normal", "OA"))
metadata$Celltype <- factor(metadata$Celltype, levels = celltype_levels)
stopifnot(
  !anyNA(metadata$condition), !anyNA(metadata$Celltype),
  all(vapply(split(as.character(metadata$condition), metadata$donor),
             function(x) length(unique(x)) == 1L, logical(1)))
)

srsc$donor <- metadata$donor
srsc$condition <- metadata$condition
srsc$Celltype <- metadata$Celltype
Idents(srsc) <- "Celltype"

availability <- data.frame(
  variable = c("donor", "disease_condition", "cell_subtype", "age", "sex", "batch", "doublet_call"),
  available = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),
  analysis_use = c(
    "Independent biological unit", "Normal versus OA covariate", "Within-subtype strata",
    "Not available; not modeled", "Not available; not modeled", "Not available; not modeled",
    "Not available; object was supplied as already cleaned"
  )
)
write.csv(availability, file.path(result_dir, "metadata_availability.csv"), row.names = FALSE)

donor_summary <- metadata %>%
  group_by(donor, condition) %>%
  summarise(
    n_cells = n(),
    median_nCount_RNA = median(nCount_RNA),
    median_nFeature_RNA = median(nFeature_RNA),
    median_percent_mt = median(percent.mt),
    median_percent_HB = if ("percent.HB" %in% names(metadata)) median(percent.HB) else NA_real_,
    .groups = "drop"
  )
donor_summary$condition <- factor(donor_summary$condition, levels = c("Normal", "OA"))
donor_summary$donor <- factor(donor_summary$donor, levels = donor_summary$donor[order(donor_summary$condition)])
write.csv(donor_summary, file.path(result_dir, "S5A_供者QC汇总.csv"), row.names = FALSE)

cell_count_table <- metadata %>%
  count(donor, condition, Celltype, name = "n_cells") %>%
  complete(
    donor = unique(metadata$donor), Celltype = factor(celltype_levels, levels = celltype_levels),
    fill = list(n_cells = 0L)
  ) %>%
  mutate(
    condition = factor(ifelse(grepl("^normal", donor), "Normal", "OA"), levels = c("Normal", "OA")),
    eligible_for_donor_subtype_analysis = n_cells >= minimum_cells_per_donor_subtype
  )
write.csv(cell_count_table, file.path(result_dir, "S5B_供者亚型细胞数.csv"), row.names = FALSE)

qc_long <- metadata %>%
  select(cell_id, donor, condition, nCount_RNA, nFeature_RNA, percent.mt) %>%
  pivot_longer(
    cols = c(nCount_RNA, nFeature_RNA, percent.mt),
    names_to = "metric", values_to = "value"
  )
qc_long$metric <- factor(
  qc_long$metric,
  levels = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  labels = c("UMI counts", "Detected genes", "Mitochondrial reads (%)")
)

p_qc <- ggplot(qc_long, aes(donor, value, fill = condition, colour = condition)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.25, alpha = 0.28) +
  geom_boxplot(width = 0.13, outlier.shape = NA, fill = "white", linewidth = 0.35) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = condition_colors, drop = FALSE) +
  scale_colour_manual(values = condition_colors, drop = FALSE) +
  labs(x = NULL, y = NULL, fill = NULL, colour = NULL) +
  theme_pub(13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
save_plot(p_qc, "补充图S5A_样本与细胞QC", figure_supp_dir, 260, 105)

p_donor_counts <- ggplot(donor_summary, aes(donor, n_cells, fill = condition)) +
  geom_col(width = 0.72, colour = "black", linewidth = 0.35) +
  geom_text(aes(label = scales::comma(n_cells)), vjust = -0.35, size = 4.1, family = base_family) +
  scale_fill_manual(values = condition_colors, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.13)), labels = scales::comma) +
  labs(x = NULL, y = "Cells", fill = NULL) +
  theme_pub(14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
save_plot(p_donor_counts, "补充图S5B_各供者细胞数", figure_supp_dir, 165, 125)

#### 2. 现有注释、批次校正降维与标志基因核查 -----------------------------------

p_umap <- DimPlot(
  srsc, reduction = "umap", group.by = "Celltype", cols = celltype_colors,
  label = TRUE, repel = TRUE, pt.size = 0.08, raster = TRUE, raster.dpi = c(1200, 1200)
) +
  labs(x = "UMAP 1", y = "UMAP 2", colour = "Cell subtype") +
  theme_pub(14) +
  theme(legend.position = "right", aspect.ratio = 1)
save_plot(p_umap, "主图1j_细胞类型UMAP", figure_main_dir, 180, 155)
save_plot(p_umap, "补充图S5C_Harmony校正UMAP", figure_supp_dir, 180, 155)

marker_panel <- list(
  EC = c("C7", "CFH", "CYTL1"),
  FC = c("COL1A1", "COL3A1", "DCN"),
  HomC = c("COL2A1", "ACAN", "SOX9"),
  ProC = c("STMN1", "TUBA1B", "HMGB2"),
  RegC = c("CHI3L1", "CLU", "PRG4"),
  preHTC = c("SPP1", "COL10A1", "ALPL"),
  HTC = c("COL10A1", "MMP13", "IBSP")
)
marker_panel <- lapply(marker_panel, intersect, y = rownames(srsc))
stopifnot(all(lengths(marker_panel) >= 2L))
marker_features <- unique(unlist(marker_panel, use.names = FALSE))
write.csv(
  data.frame(expected_subtype = rep(names(marker_panel), lengths(marker_panel)),
             marker = unlist(marker_panel, use.names = FALSE)),
  file.path(result_dir, "S5D_细胞亚型标志基因列表.csv"), row.names = FALSE
)

p_markers <- DotPlot(
  srsc, features = marker_features, group.by = "Celltype", assay = "RNA",
  cols = c("#E7E7E7", "#B2182B"), dot.scale = 8, scale.by = "radius"
) +
  scale_size(range = c(0, 7.5)) +
  labs(x = NULL, y = NULL, colour = "Scaled expression", size = "Expressing cells (%)") +
  theme_pub(13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.35)
  )
save_plot(p_markers, "补充图S5D_细胞亚型标志基因DotPlot", figure_supp_dir, 260, 145)

#### 3. 供者级细胞组成 ------------------------------------------------------

cell_proportions <- cell_count_table %>%
  group_by(donor) %>%
  mutate(total_cells = sum(n_cells), proportion = n_cells / total_cells) %>%
  ungroup()
write.csv(cell_proportions, file.path(result_dir, "S5E_供者级细胞亚型比例.csv"), row.names = FALSE)

proportion_tests <- bind_rows(lapply(celltype_levels, function(subtype) {
  current <- filter(cell_proportions, Celltype == subtype)
  group_test(current, "proportion", subtype)
}))
proportion_tests$adjusted_p_value <- p.adjust(proportion_tests$p_value, method = "BH")
write.csv(proportion_tests, file.path(result_dir, "S5E_供者级细胞比例组间检验.csv"), row.names = FALSE)

p_proportions <- ggplot(cell_proportions, aes(condition, 100 * proportion, colour = condition)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, fill = "white", linewidth = 0.45) +
  geom_point(size = 2.4, position = position_jitter(width = 0.08, height = 0)) +
  facet_wrap(~Celltype, ncol = 4, scales = "free_y") +
  scale_colour_manual(values = condition_colors, drop = FALSE) +
  labs(x = NULL, y = "Cells per donor (%)", colour = NULL) +
  theme_pub(13) +
  theme(legend.position = "top")
save_plot(p_proportions, "补充图S5E_供者级细胞亚型比例", figure_supp_dir, 245, 175)

#### 4. DPP3 定位：所有细胞可视化，供者为定量单位 -----------------------------

rna_data <- GetAssayData(srsc, assay = "RNA", slot = "data")
metadata$DPP3_log1p <- as.numeric(rna_data["DPP3", ])
metadata$DPP3_detected <- metadata$DPP3_log1p > 0
srsc$DPP3_detected <- metadata$DPP3_detected

p_dpp3_feature <- FeaturePlot(
  srsc, features = "DPP3", reduction = "umap", order = TRUE,
  min.cutoff = "q05", max.cutoff = "q99.5", cols = c("#E7E7E7", "#F4A261", "#B2182B"),
  pt.size = 0.08, raster = TRUE, raster.dpi = c(1200, 1200)
) +
  labs(x = "UMAP 1", y = "UMAP 2", colour = "DPP3") +
  theme_pub(14) +
  theme(legend.position = "right", aspect.ratio = 1)
save_plot(p_dpp3_feature, "主图1j_DPP3表达定位", figure_main_dir, 180, 155)

p_dpp3_split <- FeaturePlot(
  srsc, features = "DPP3", reduction = "umap", split.by = "condition", keep.scale = "all",
  order = TRUE, min.cutoff = "q05", max.cutoff = "q99.5",
  cols = c("#E7E7E7", "#F4A261", "#B2182B"), pt.size = 0.08,
  raster = TRUE, raster.dpi = c(1200, 1200)
) &
  theme_pub(13) &
  theme(aspect.ratio = 1)
save_plot(p_dpp3_split, "补充图S5F_DPP3按Normal与OA拆分定位", figure_supp_dir, 285, 145)

dpp3_cell_summary <- metadata %>%
  group_by(donor, condition, Celltype) %>%
  summarise(
    n_cells = n(),
    mean_log_normalized_expression = mean(DPP3_log1p),
    fraction_positive = mean(DPP3_detected),
    .groups = "drop"
  ) %>%
  mutate(eligible_for_donor_subtype_analysis = n_cells >= minimum_cells_per_donor_subtype)
write.csv(dpp3_cell_summary, file.path(result_dir, "S5F_DPP3供者亚型表达与阳性比例.csv"), row.names = FALSE)

p_dpp3_mean <- ggplot(dpp3_cell_summary, aes(condition, mean_log_normalized_expression, colour = condition)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, fill = "white", linewidth = 0.4) +
  geom_point(aes(shape = eligible_for_donor_subtype_analysis), size = 2.2,
             position = position_jitter(width = 0.08, height = 0)) +
  facet_wrap(~Celltype, ncol = 4, scales = "free_y") +
  scale_colour_manual(values = condition_colors, drop = FALSE) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = paste0("At least ", minimum_cells_per_donor_subtype, " cells")) +
  labs(x = NULL, y = "Mean log-normalized DPP3", colour = NULL) +
  theme_pub(12) +
  theme(legend.position = "top")

p_dpp3_fraction <- ggplot(dpp3_cell_summary, aes(condition, 100 * fraction_positive, colour = condition)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, fill = "white", linewidth = 0.4) +
  geom_point(aes(shape = eligible_for_donor_subtype_analysis), size = 2.2,
             position = position_jitter(width = 0.08, height = 0)) +
  facet_wrap(~Celltype, ncol = 4, scales = "free_y") +
  scale_colour_manual(values = condition_colors, drop = FALSE) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  labs(x = NULL, y = "DPP3-positive cells (%)", colour = NULL) +
  theme_pub(12) +
  theme(legend.position = "top")

save_plot(
  p_dpp3_mean / p_dpp3_fraction + plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = 16, face = "bold")),
  "补充图S5G_DPP3供者平均表达与阳性比例", figure_supp_dir, 245, 310
)

#### 5. 供者×细胞亚型伪批量与 OA/Normal 差异 -------------------------------

rna_counts <- GetAssayData(srsc, assay = "RNA", slot = "counts")
pb_key <- paste(metadata$donor, metadata$Celltype, sep = "__")
pb_levels <- unique(pb_key)
pb_factor <- factor(pb_key, levels = pb_levels)
aggregation_matrix <- Matrix::sparse.model.matrix(~0 + pb_factor)
colnames(aggregation_matrix) <- sub("^pb_factor", "", colnames(aggregation_matrix))
pseudobulk_counts <- rna_counts %*% aggregation_matrix
colnames(pseudobulk_counts) <- pb_levels

pb_meta <- data.frame(key = pb_levels, stringsAsFactors = FALSE)
pb_meta$donor <- sub("__.*$", "", pb_meta$key)
pb_meta$Celltype <- factor(sub("^.*__", "", pb_meta$key), levels = celltype_levels)
donor_condition <- setNames(as.character(donor_summary$condition), as.character(donor_summary$donor))
pb_meta$condition <- factor(unname(donor_condition[pb_meta$donor]), levels = c("Normal", "OA"))
pb_meta$n_cells <- as.integer(table(pb_factor)[pb_levels])
pb_meta$eligible <- pb_meta$n_cells >= minimum_cells_per_donor_subtype
stopifnot(!anyNA(pb_meta$condition), identical(colnames(pseudobulk_counts), pb_meta$key))
write.csv(pb_meta, file.path(result_dir, "供者亚型伪批量样本信息.csv"), row.names = FALSE)

write_csv_gz(
  data.frame(gene = rownames(pseudobulk_counts), as.matrix(pseudobulk_counts), check.names = FALSE),
  file.path(result_dir, "供者亚型伪批量原始counts.csv.gz")
)

pseudobulk_logcpm <- matrix(
  NA_real_, nrow = nrow(pseudobulk_counts), ncol = ncol(pseudobulk_counts),
  dimnames = dimnames(pseudobulk_counts)
)
dpp3_effects <- list()
de_results <- list()

for (subtype in celltype_levels) {
  selected <- which(pb_meta$Celltype == subtype & pb_meta$eligible)
  current_meta <- pb_meta[selected, , drop = FALSE]
  group_n <- table(current_meta$condition)
  if (length(selected) >= 2L) {
    y_all <- DGEList(counts = pseudobulk_counts[, selected, drop = FALSE])
    y_all <- calcNormFactors(y_all)
    pseudobulk_logcpm[, selected] <- cpm(y_all, log = TRUE, prior.count = 2)
  }

  if (length(group_n) < 2L || any(group_n < 2L)) {
    dpp3_effects[[subtype]] <- data.frame(
      Celltype = subtype, n_Normal = unname(group_n["Normal"]), n_OA = unname(group_n["OA"]),
      log2FC_OA_vs_Normal = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      p_value = NA_real_, adjusted_p_value = NA_real_, tested = FALSE,
      reason = "Fewer than two eligible donors in at least one condition"
    )
    next
  }

  y <- DGEList(counts = pseudobulk_counts[, selected, drop = FALSE])
  design <- model.matrix(~condition, data = current_meta)
  keep_genes <- filterByExpr(y, design = design)
  y <- y[keep_genes, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  voom_data <- voom(y, design, plot = FALSE)
  fit <- eBayes(lmFit(voom_data, design), trend = FALSE, robust = TRUE)
  coefficient <- "conditionOA"
  table_de <- topTable(fit, coef = coefficient, number = Inf, sort.by = "none")
  table_de$gene <- rownames(table_de)
  moderated_se <- abs(table_de$logFC / table_de$t)
  moderated_se[!is.finite(moderated_se)] <- NA_real_
  degrees_freedom <- fit$df.total[match(table_de$gene, rownames(fit$coefficients))]
  critical_value <- qt(0.975, df = degrees_freedom)
  table_de$ci_low <- table_de$logFC - critical_value * moderated_se
  table_de$ci_high <- table_de$logFC + critical_value * moderated_se
  table_de$Celltype <- subtype
  table_de$n_Normal <- unname(group_n["Normal"])
  table_de$n_OA <- unname(group_n["OA"])
  table_de <- table_de[, c(
    "gene", "Celltype", "n_Normal", "n_OA", "logFC", "ci_low", "ci_high",
    "AveExpr", "t", "P.Value", "adj.P.Val", "B"
  )]
  names(table_de)[names(table_de) == "P.Value"] <- "p_value"
  names(table_de)[names(table_de) == "adj.P.Val"] <- "adjusted_p_value"
  write_csv_gz(table_de, file.path(de_dir, paste0("伪批量差异_", subtype, ".csv.gz")))
  de_results[[subtype]] <- table_de

  if ("DPP3" %in% table_de$gene) {
    dpp3_row <- table_de[table_de$gene == "DPP3", ]
    dpp3_effects[[subtype]] <- data.frame(
      Celltype = subtype, n_Normal = unname(group_n["Normal"]), n_OA = unname(group_n["OA"]),
      log2FC_OA_vs_Normal = dpp3_row$logFC,
      ci_low = dpp3_row$ci_low, ci_high = dpp3_row$ci_high,
      p_value = dpp3_row$p_value, adjusted_p_value = dpp3_row$adjusted_p_value,
      tested = TRUE, reason = "Tested by donor-level limma-voom"
    )
  } else {
    dpp3_effects[[subtype]] <- data.frame(
      Celltype = subtype, n_Normal = unname(group_n["Normal"]), n_OA = unname(group_n["OA"]),
      log2FC_OA_vs_Normal = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      p_value = NA_real_, adjusted_p_value = NA_real_, tested = FALSE,
      reason = "DPP3 did not pass pseudobulk expression filtering"
    )
  }
}

dpp3_effects <- bind_rows(dpp3_effects)
dpp3_effects$Celltype <- factor(dpp3_effects$Celltype, levels = rev(celltype_levels))
write.csv(dpp3_effects, file.path(result_dir, "S5H_DPP3伪批量OA_vs_Normal效应.csv"), row.names = FALSE)
write_csv_gz(
  data.frame(gene = rownames(pseudobulk_logcpm), pseudobulk_logcpm, check.names = FALSE),
  file.path(result_dir, "供者亚型伪批量logCPM.csv.gz")
)

dpp3_pb <- pb_meta
dpp3_pb$DPP3_logCPM <- as.numeric(pseudobulk_logcpm["DPP3", ])
dpp3_pb <- dpp3_pb[dpp3_pb$eligible & is.finite(dpp3_pb$DPP3_logCPM), ]
write.csv(dpp3_pb, file.path(result_dir, "S5G_DPP3供者伪批量logCPM.csv"), row.names = FALSE)

p_dpp3_pb <- ggplot(dpp3_pb, aes(condition, DPP3_logCPM, colour = condition)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, fill = "white", linewidth = 0.42) +
  geom_point(size = 2.35, position = position_jitter(width = 0.08, height = 0)) +
  facet_wrap(~Celltype, ncol = 3, scales = "free_y") +
  scale_colour_manual(values = condition_colors, drop = FALSE) +
  labs(x = NULL, y = expression("DPP3 pseudobulk log"[2]*"CPM"), colour = NULL) +
  theme_pub(12) +
  theme(legend.position = "top")
save_plot(p_dpp3_pb, "补充图S5G_DPP3供者伪批量表达", figure_supp_dir, 215, 180)

p_dpp3_forest <- ggplot(filter(dpp3_effects, tested), aes(log2FC_OA_vs_Normal, Celltype)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "#666666", linewidth = 0.45) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.16, linewidth = 0.55, colour = "#444444") +
  geom_point(aes(fill = adjusted_p_value < 0.05), shape = 21, size = 3.4, colour = "black") +
  scale_fill_manual(values = c(`TRUE` = "#C44E52", `FALSE` = "white"), name = "FDR < 0.05") +
  labs(x = expression(log[2]*" fold change (OA / Normal)"), y = NULL) +
  theme_pub(14) +
  theme(legend.position = "top")
save_plot(p_dpp3_forest, "补充图S5H_DPP3亚型内OA_vs_Normal效应", figure_supp_dir, 165, 125)

#### 6. 小鼠 RNA 衰老签名的一对一同源映射与 UCell 投射 ----------------------

mouse_sets <- read.csv(mouse_signature_file, stringsAsFactors = FALSE)
stopifnot(all(c("module", "gene") %in% names(mouse_sets)))
mouse_sets <- mouse_sets[nzchar(mouse_sets$module) & nzchar(mouse_sets$gene), ]

ortholog_table <- babelgene::orthologs(
  unique(mouse_sets$gene), species = "mouse", human = FALSE,
  min_support = 3, top = FALSE
)
one_mouse <- names(which(table(ortholog_table$symbol) == 1L))
one_human <- names(which(table(ortholog_table$human_symbol) == 1L))
ortholog_one_to_one <- ortholog_table %>%
  filter(symbol %in% one_mouse, human_symbol %in% one_human) %>%
  distinct(symbol, human_symbol, .keep_all = TRUE)

signature_mapping <- mouse_sets %>%
  left_join(
    select(ortholog_one_to_one, mouse_gene = symbol, human_gene = human_symbol, support_n),
    by = c("gene" = "mouse_gene")
  ) %>%
  mutate(
    mapped_one_to_one = !is.na(human_gene),
    present_in_scRNA = human_gene %in% rownames(srsc),
    excluded_as_DPP3 = human_gene == "DPP3"
  )
write.csv(signature_mapping, file.path(result_dir, "小鼠RNA衰老签名_人一对一同源映射.csv"), row.names = FALSE)

signature_coverage <- signature_mapping %>%
  group_by(module) %>%
  summarise(
    mouse_input_genes = n(),
    one_to_one_mapped = sum(mapped_one_to_one),
    present_in_scRNA = sum(present_in_scRNA & !excluded_as_DPP3),
    mapping_coverage = one_to_one_mapped / mouse_input_genes,
    scRNA_coverage = present_in_scRNA / mouse_input_genes,
    .groups = "drop"
  )
write.csv(signature_coverage, file.path(result_dir, "S5I_跨物种衰老签名覆盖率.csv"), row.names = FALSE)

score_sets <- split(
  signature_mapping$human_gene[
    signature_mapping$present_in_scRNA & !signature_mapping$excluded_as_DPP3
  ],
  signature_mapping$module[
    signature_mapping$present_in_scRNA & !signature_mapping$excluded_as_DPP3
  ]
)
required_modules <- c(
  "Inflammation", "Cartilage_homeostasis", "ECM_homeostasis",
  "Mitochondrial_homeostasis"
)
score_sets <- score_sets[required_modules]
stopifnot(all(lengths(score_sets) >= 50L), !"DPP3" %in% unlist(score_sets, use.names = FALSE))

senmayo_mouse <- strsplit(
  "Acvr1b Ang Angpt1 Angptl4 Areg Axl Bex3 Bmp2 Bmp6 C3 Ccl1 Ccl13 Ccl16 Ccl2 Ccl20 Ccl24 Ccl26 Ccl3 Ccl3l1 Ccl4 Ccl5 Ccl7 Ccl8 Cd55 Cd9 Csf1 Csf2 Csf2rb Cst4 Ctnnb1 Ctsb Cxcl1 Cxcl10 Cxcl12 Cxcl16 Cxcl2 Cxcl3 Cxcl8 Cxcr2 Dkk1 Edn1 Egf Egfr Ereg Esm1 Ets2 Fas Fgf1 Fgf2 Fgf7 Gdf15 Gem Gmfg Hgf Hmgb1 Icam1 Icam5 Igf1 Igfbp1 Igfbp2 Igfbp3 Igfbp4 Igfbp5 Igfbp6 Igfbp7 Il10 Il13 Il15 Il18 Il1a Il1b Il2 Il6 Il6st Il7 Inha Iqgap2 Itga2 Itpka Jun Kitl Lcp1 Mif Mmp10 Mmp12 Mmp13 Mmp14 Mmp2 Mmp3 Mmp9 Nap1l4 Nrg1 Pappa Pecam1 Pgf Plat Plau Plaur Ptbp1 Ptger2 Ptges Rps6ka5 Scamp4 Selplg Sema3f Serpinb3c Serpine1 Serpine2 Spp1 Spx Timp2 Tnf Tnfrsf10c Tnfrsf11b Tnfrsf1a Tnfrsf1b Tubgcp2 Vegfa Vegfc Vgf Wnt16 Wnt2",
  " ", fixed = TRUE
)[[1]]
senmayo_orthologs <- babelgene::orthologs(
  senmayo_mouse, species = "mouse", human = FALSE, min_support = 3, top = FALSE
)
sen_one_mouse <- names(which(table(senmayo_orthologs$symbol) == 1L))
sen_one_human <- names(which(table(senmayo_orthologs$human_symbol) == 1L))
senmayo_map <- senmayo_orthologs %>%
  filter(symbol %in% sen_one_mouse, human_symbol %in% sen_one_human) %>%
  distinct(symbol, human_symbol, .keep_all = TRUE) %>%
  mutate(present_in_scRNA = human_symbol %in% rownames(srsc))
write.csv(senmayo_map, file.path(result_dir, "SenMayo_人一对一同源映射.csv"), row.names = FALSE)
score_sets$SenMayo <- setdiff(intersect(senmayo_map$human_symbol, rownames(srsc)), "DPP3")
stopifnot(length(score_sets$SenMayo) >= 50L)

srsc <- AddModuleScore_UCell(
  srsc, features = score_sets, assay = "RNA", slot = "data", name = "_UCell",
  maxRank = 1500, chunk.size = 2000, ncores = 1, force.gc = TRUE
)
score_columns <- paste0(names(score_sets), "_UCell")
stopifnot(all(score_columns %in% colnames(srsc@meta.data)))
metadata <- srsc@meta.data
metadata$cell_id <- rownames(metadata)
metadata$condition <- factor(metadata$condition, levels = c("Normal", "OA"))
metadata$Celltype <- factor(metadata$Celltype, levels = celltype_levels)
metadata$DPP3_log1p <- as.numeric(rna_data["DPP3", ])
metadata$DPP3_detected <- metadata$DPP3_log1p > 0

score_aggregates <- metadata %>%
  group_by(donor, condition, Celltype) %>%
  summarise(
    n_cells = n(),
    across(all_of(score_columns), list(median = median, mean = mean), .names = "{.col}_{.fn}"),
    .groups = "drop"
  ) %>%
  mutate(eligible = n_cells >= minimum_cells_per_donor_subtype)

score_aggregates <- score_aggregates %>%
  group_by(Celltype) %>%
  mutate(
    Inflammation_z_median = z_safe(Inflammation_UCell_median),
    Cartilage_z_median = z_safe(Cartilage_homeostasis_UCell_median),
    ECM_z_median = z_safe(ECM_homeostasis_UCell_median),
    Mitochondrial_z_median = z_safe(Mitochondrial_homeostasis_UCell_median),
    SenMayo_z_median = z_safe(SenMayo_UCell_median),
    Inflammation_z_mean = z_safe(Inflammation_UCell_mean),
    Cartilage_z_mean = z_safe(Cartilage_homeostasis_UCell_mean),
    ECM_z_mean = z_safe(ECM_homeostasis_UCell_mean),
    Mitochondrial_z_mean = z_safe(Mitochondrial_homeostasis_UCell_mean)
  ) %>%
  ungroup() %>%
  mutate(
    RNA_aging_burden_median = rowMeans(cbind(
      Inflammation_z_median, -Cartilage_z_median, -ECM_z_median, -Mitochondrial_z_median
    )),
    RNA_aging_burden_plus_SenMayo = rowMeans(cbind(
      SenMayo_z_median, Inflammation_z_median,
      -Cartilage_z_median, -ECM_z_median, -Mitochondrial_z_median
    )),
    RNA_aging_burden_mean = rowMeans(cbind(
      Inflammation_z_mean, -Cartilage_z_mean, -ECM_z_mean, -Mitochondrial_z_mean
    ))
  )
write.csv(score_aggregates, file.path(result_dir, "S5I_供者亚型UCell衰老评分.csv"), row.names = FALSE)

score_eligible <- filter(score_aggregates, eligible)
aging_group_tests <- bind_rows(lapply(celltype_levels, function(subtype) {
  current <- filter(score_eligible, Celltype == subtype)
  group_test(current, "RNA_aging_burden_median", subtype)
}))
aging_group_tests$adjusted_p_value <- p.adjust(aging_group_tests$p_value, method = "BH")
write.csv(aging_group_tests, file.path(result_dir, "S5I_衰老负荷OA_vs_Normal供者级检验.csv"), row.names = FALSE)

p_aging <- ggplot(score_eligible, aes(condition, RNA_aging_burden_median, colour = condition)) +
  geom_hline(yintercept = 0, colour = "#BDBDBD", linewidth = 0.35) +
  geom_boxplot(width = 0.5, outlier.shape = NA, fill = "white", linewidth = 0.42) +
  geom_point(size = 2.35, position = position_jitter(width = 0.08, height = 0)) +
  facet_wrap(~Celltype, ncol = 3, scales = "free_y") +
  scale_colour_manual(values = condition_colors, drop = FALSE) +
  labs(x = NULL, y = "RNA-derived cartilage aging burden", colour = NULL) +
  theme_pub(12) +
  theme(legend.position = "top")
save_plot(p_aging, "补充图S5I_小鼠RNA衰老签名人单细胞投射", figure_supp_dir, 215, 180)

main_1k <- p_dpp3_pb + p_aging +
  plot_layout(ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "top", plot.tag = element_text(size = 16, face = "bold"))
save_plot(main_1k, "主图1k_DPP3供者伪批量与RNA衰老签名", figure_main_dir, 360, 185)

#### 7. 亚型内 DPP3—衰老负荷供者级关联 --------------------------------------

association_data <- score_eligible %>%
  left_join(select(dpp3_pb, donor, Celltype, DPP3_logCPM), by = c("donor", "Celltype")) %>%
  filter(is.finite(DPP3_logCPM), is.finite(RNA_aging_burden_median))
write.csv(association_data, file.path(result_dir, "S5J_DPP3与衰老负荷供者级数据.csv"), row.names = FALSE)

association_results <- bind_rows(lapply(celltype_levels, function(subtype) {
  current <- filter(association_data, Celltype == subtype)
  group_n <- table(current$condition)
  if (nrow(current) < 5L || length(group_n) < 2L || any(group_n < 2L)) {
    return(data.frame(
      Celltype = subtype, n_donors = nrow(current), n_Normal = unname(group_n["Normal"]),
      n_OA = unname(group_n["OA"]), spearman_rho = NA_real_, spearman_p = NA_real_,
      disease_adjusted_partial_rho = NA_real_, disease_adjusted_p = NA_real_,
      tested = FALSE, reason = "Insufficient donor coverage"
    ))
  }
  simple_test <- suppressWarnings(cor.test(
    current$DPP3_logCPM, current$RNA_aging_burden_median,
    method = "spearman", exact = FALSE
  ))
  dpp3_residual <- residuals(lm(rank(DPP3_logCPM) ~ condition, data = current))
  burden_residual <- residuals(lm(rank(RNA_aging_burden_median) ~ condition, data = current))
  adjusted_test <- suppressWarnings(cor.test(dpp3_residual, burden_residual, method = "pearson"))
  data.frame(
    Celltype = subtype, n_donors = nrow(current), n_Normal = unname(group_n["Normal"]),
    n_OA = unname(group_n["OA"]), spearman_rho = unname(simple_test$estimate),
    spearman_p = simple_test$p.value,
    disease_adjusted_partial_rho = unname(adjusted_test$estimate),
    disease_adjusted_p = adjusted_test$p.value,
    tested = TRUE, reason = "Donor-level within-subtype association"
  )
}))
association_results$adjusted_p_value <- p.adjust(association_results$disease_adjusted_p, method = "BH")
write.csv(association_results, file.path(result_dir, "S5J_DPP3与衰老负荷亚型内关联.csv"), row.names = FALSE)

association_plot_data <- association_data %>%
  filter(Celltype %in% association_results$Celltype[association_results$tested])
association_labels <- association_results %>%
  filter(tested) %>%
  mutate(label = sprintf("partial rho = %.2f\nn = %d donors", disease_adjusted_partial_rho, n_donors))

p_association <- ggplot(
  association_plot_data,
  aes(DPP3_logCPM, RNA_aging_burden_median, colour = condition)
) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "#555555", fill = "#D9D9D9", linewidth = 0.55) +
  geom_point(size = 2.6) +
  geom_text(
    data = association_labels, aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE, hjust = -0.08, vjust = 1.15, size = 3.4, family = base_family
  ) +
  facet_wrap(~Celltype, ncol = 3, scales = "free") +
  scale_colour_manual(values = condition_colors, drop = FALSE) +
  labs(
    x = expression("DPP3 pseudobulk log"[2]*"CPM"),
    y = "RNA-derived cartilage aging burden", colour = NULL
  ) +
  theme_pub(12) +
  theme(legend.position = "top")
save_plot(p_association, "补充图S5J_DPP3与衰老负荷供者级相关", figure_supp_dir, 220, 190)

#### 8. 签名定义与聚合方式敏感性 ---------------------------------------------

sensitivity_results <- bind_rows(lapply(celltype_levels, function(subtype) {
  current <- filter(score_eligible, Celltype == subtype)
  if (nrow(current) < 4L) return(NULL)
  data.frame(
    Celltype = subtype,
    comparison = c("Primary vs +SenMayo", "Median vs mean aggregation"),
    spearman_rho = c(
      suppressWarnings(cor(current$RNA_aging_burden_median,
                           current$RNA_aging_burden_plus_SenMayo,
                           method = "spearman", use = "pairwise.complete.obs")),
      suppressWarnings(cor(current$RNA_aging_burden_median,
                           current$RNA_aging_burden_mean,
                           method = "spearman", use = "pairwise.complete.obs"))
    ),
    n_donors = nrow(current)
  )
}))
write.csv(sensitivity_results, file.path(result_dir, "S5K_衰老评分敏感性分析.csv"), row.names = FALSE)

p_sensitivity <- ggplot(sensitivity_results, aes(comparison, Celltype, fill = spearman_rho)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("%.2f", spearman_rho)), size = 4, family = base_family) +
  scale_fill_gradient2(
    low = "#3C78A8", mid = "white", high = "#C44E52", midpoint = 0,
    limits = c(-1, 1), name = "Spearman rho"
  ) +
  labs(x = NULL, y = NULL) +
  theme_pub(13) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), axis.line = element_blank(), axis.ticks = element_blank())
save_plot(p_sensitivity, "补充图S5K_衰老签名定义敏感性", figure_supp_dir, 175, 135)

#### 9. 探索性亚型内 DPP3 相关通路 -------------------------------------------

pathway_association <- list()
fgsea_sets <- score_sets

for (subtype in association_results$Celltype[association_results$tested]) {
  subtype <- as.character(subtype)
  current_keys <- association_data$key[association_data$Celltype == subtype]
  current_keys <- intersect(current_keys, colnames(pseudobulk_logcpm))
  current_meta <- pb_meta[match(current_keys, pb_meta$key), , drop = FALSE]
  current_matrix <- pseudobulk_logcpm[, current_keys, drop = FALSE]
  complete_columns <- colSums(is.finite(current_matrix)) == nrow(current_matrix)
  current_matrix <- current_matrix[, complete_columns, drop = FALSE]
  current_meta <- current_meta[complete_columns, , drop = FALSE]
  if (ncol(current_matrix) < 5L || length(unique(current_meta$condition)) < 2L) next

  ranked_matrix <- t(apply(current_matrix, 1, rank, ties.method = "average"))
  design <- model.matrix(~condition, data = current_meta)
  gene_residuals <- t(lm.fit(design, t(ranked_matrix))$residuals)
  dpp3_residuals <- residuals(lm(rank(current_matrix["DPP3", ]) ~ condition, data = current_meta))
  denominator <- sqrt(rowSums(gene_residuals^2) * sum(dpp3_residuals^2))
  gene_rho <- rowSums(sweep(gene_residuals, 2, dpp3_residuals, `*`)) / denominator
  gene_rho[!is.finite(gene_rho)] <- 0
  gene_rho["DPP3"] <- NA_real_
  ranks <- sort(gene_rho[is.finite(gene_rho)], decreasing = TRUE)

  set.seed(1)
  pathway_result <- as.data.frame(fgseaMultilevel(
    pathways = fgsea_sets, stats = ranks, minSize = 10, maxSize = 500, eps = 0
  ))
  pathway_result$Celltype <- subtype
  pathway_result$n_donors <- ncol(current_matrix)
  pathway_result$leadingEdge <- vapply(pathway_result$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  pathway_association[[subtype]] <- pathway_result
}

pathway_association <- bind_rows(pathway_association)
if (nrow(pathway_association) > 0L) {
  pathway_association$pathway <- factor(
    pathway_association$pathway,
    levels = c("SenMayo", "Inflammation", "Cartilage_homeostasis", "ECM_homeostasis", "Mitochondrial_homeostasis")
  )
  pathway_association$Celltype <- factor(pathway_association$Celltype, levels = celltype_levels)
  write.csv(pathway_association, file.path(result_dir, "S5L_亚型内DPP3相关通路_fgsea.csv"), row.names = FALSE)

  p_pathways <- ggplot(pathway_association, aes(Celltype, pathway, fill = NES)) +
    geom_tile(colour = "white", linewidth = 0.9) +
    geom_point(
      aes(size = -log10(pmax(padj, 1e-12))), shape = 21,
      fill = "white", colour = "black", stroke = 0.35
    ) +
    scale_fill_gradient2(low = "#3C78A8", mid = "white", high = "#C44E52", midpoint = 0, name = "NES") +
    scale_size_continuous(range = c(1.2, 5), name = expression(-log[10]~FDR)) +
    labs(x = NULL, y = NULL) +
    theme_pub(12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.line = element_blank(), axis.ticks = element_blank())
  save_plot(p_pathways, "补充图S5L_探索性亚型内DPP3相关通路", figure_supp_dir, 205, 140)
}

#### 10. 完整源数据、图件清单与可复现记录 ------------------------------------

cell_metadata_export <- metadata %>%
  select(
    cell_id, orig.ident, donor, condition, Celltype, nCount_RNA, nFeature_RNA,
    percent.mt, any_of("percent.HB"), DPP3_log1p, DPP3_detected, all_of(score_columns)
  )
write_csv_gz(cell_metadata_export, file.path(result_dir, "细胞元数据与UCell评分.csv.gz"))

figure_manifest <- data.frame(
  figure = c(
    "Main 1j cell-type UMAP", "Main 1j DPP3 FeaturePlot", "Main 1k donor pseudobulk/signature",
    paste0("S5", c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"))
  ),
  evidence_role = c(
    "Cell-state localization", "DPP3 localization", "Donor-level quantitative projection",
    "QC", "Donor cell counts", "Corrected embedding", "Annotation validation",
    "Donor-level composition", "Condition-split localization", "Donor pseudobulk DPP3",
    "Within-subtype OA versus Normal", "Cross-species aging score", "DPP3-score association",
    "Signature sensitivity", "Exploratory pathway association"
  ),
  independent_unit = c(
    "Cells shown; no inferential test", "Cells shown; no inferential test", "Donor",
    "Cell descriptive", "Donor", "Cells shown", "Cells shown", "Donor", "Cells shown",
    "Donor", "Donor", "Donor", "Donor", "Donor", "Donor"
  )
)
write.csv(figure_manifest, file.path(result_dir, "figure_manifest.csv"), row.names = FALSE)

analysis_summary <- data.frame(
  item = c(
    "Cells", "Features", "Donors", "Normal donors", "OA donors", "Cell subtypes",
    "Minimum cells per donor-subtype", "Age metadata", "Sex metadata", "Batch metadata",
    "Doublet metadata", "Primary aging signature", "Primary inference unit"
  ),
  value = c(
    ncol(srsc), nrow(srsc), length(unique(metadata$donor)),
    length(unique(metadata$donor[metadata$condition == "Normal"])),
    length(unique(metadata$donor[metadata$condition == "OA"])),
    length(unique(metadata$Celltype)), minimum_cells_per_donor_subtype,
    "Unavailable", "Unavailable", "Unavailable", "Unavailable",
    "RNA-derived mouse cartilage aging modules mapped one-to-one to human",
    "Donor within cell subtype"
  )
)
write.csv(analysis_summary, file.path(result_dir, "analysis_summary.csv"), row.names = FALSE)

qa_lines <- c(
  "DPP3 scRNA-seq analysis QA notes",
  "",
  "Figure contract",
  "Core conclusion: DPP3 is localized across human cartilage cell states and is evaluated against OA status and an RNA-derived cartilage aging signature without treating cells as independent replicates.",
  "Archetype: quantitative grid.",
  "Backend: R only.",
  "Primary unit: donor within annotated cell subtype.",
  "",
  "Data integrity",
  sprintf("All %s cells and %s genes were retained for visualization and descriptive summaries.", scales::comma(ncol(srsc)), scales::comma(nrow(srsc))),
  sprintf("Inferential donor-subtype analyses required at least %d cells in that donor-subtype combination.", minimum_cells_per_donor_subtype),
  "No additional cell filtering, re-clustering, or re-annotation was performed because the supplied object was already cleaned and annotated.",
  "",
  "Reviewer-risk boundaries",
  "The cohort contains 3 Normal and 4 OA donors; OA is a disease contrast and must not be described as physiological aging.",
  "Age, sex, batch, and doublet-call metadata were unavailable and therefore were not invented or modeled.",
  "DPP3 is sparsely detected at cell level; donor-pseudobulk logCPM and fraction-positive summaries are both supplied.",
  "The projected signature is RNA-derived because matched proteomics were not available at this stage.",
  "All DPP3-associated pathway results are exploratory because within-subtype donor n is small.",
  "No DPP3-high versus DPP3-low cell-cluster comparison was performed."
)
writeLines(qa_lines, file.path(output_root, "QA_notes.md"))
writeLines(capture.output(sessionInfo()), file.path(result_dir, "sessionInfo.txt"))

stopifnot(
  nrow(donor_summary) == 7L,
  nrow(cell_proportions) == 7L * length(celltype_levels),
  all(signature_coverage$present_in_scRNA >= 50L),
  all(file.exists(file.path(figure_main_dir, c(
    "主图1j_细胞类型UMAP.png", "主图1j_DPP3表达定位.png",
    "主图1k_DPP3供者伪批量与RNA衰老签名.png"
  ))))
)

message("DPP3 single-cell analysis completed: ", output_root)
