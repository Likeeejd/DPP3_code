#!/usr/bin/env Rscript

#### Fig 1 #####
# 本轮只优化截图中的个体表达图和两张效应量统计图，不重新定义统计方法。
# 在RStudio中Source即可；读取同目录已保存的完整分析对象，不依赖当前工作区。
suppressWarnings(invisible(Sys.setlocale("LC_ALL", "zh_CN.UTF-8")))
required_packages <- c("ggplot2", "patchwork", "ashr", "systemfonts", "png")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("请先安装R包：", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages(library(ggplot2))

OUTPUT_DIR <- "/Volumes/zzData/DPP3/Figure-原始数据整理/Figure 1/Figure1-转录组-蛋白组联合分析"
VERSION <- "-v6"
font_pt <- 12
text_size <- font_pt / ggplot2::.pt
SEED <- 20260828L
stopifnot(dir.exists(OUTPUT_DIR), capabilities("cairo"), any(systemfonts::system_fonts()$family == "Arial"))

input_file <- file.path(OUTPUT_DIR, "MAIN_analysis_objects-v3.rds")
stopifnot(file.exists(input_file))
analysis <- readRDS(input_file)
dpp3_values <- analysis$DPP3_values
dpp3_effects <- analysis$DPP3_effects
module_effects <- analysis$module_effects
stopifnot(nrow(dpp3_values) == 24L, all(table(dpp3_values$assay, dpp3_values$group) == 6L),
  nrow(dpp3_effects) == 2L, nrow(module_effects) == 10L)

group_colors <- c(Young = "#2F7FB7", Old = "#D95A59")
theme_pub <- function() {
  theme_classic(base_size = font_pt, base_family = "Arial") + theme(
    text = element_text(family = "Arial", size = font_pt, colour = "black"),
    axis.title = element_text(size = font_pt), axis.text = element_text(size = font_pt, colour = "black"),
    axis.line = element_line(linewidth = 0.45, colour = "black"), axis.ticks = element_line(linewidth = 0.45, colour = "black"),
    legend.title = element_text(size = font_pt), legend.text = element_text(size = font_pt),
    strip.text = element_text(size = font_pt), plot.tag = element_text(size = font_pt),
    plot.title = element_text(size = font_pt, hjust = 0.5), plot.subtitle = element_blank(),
    plot.caption = element_text(size = font_pt),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.box.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA),
    strip.background = element_rect(fill = "transparent", colour = NA),
    panel.grid = element_blank(), plot.margin = margin(5, 7, 5, 7))
}
theme_set(theme_pub())

format_fdr <- function(x) ifelse(x < 1e-4, format(x, scientific = TRUE, digits = 3), sprintf("%.4f", x))
expression_panel <- function(d, assay_name, y_label) {
  e <- dpp3_effects[dpp3_effects$assay == assay_name, ]
  r <- range(d$expression); span <- diff(r)
  ggplot(d, aes(factor(group, levels = c("Young", "Old")), expression, colour = group)) +
    geom_boxplot(aes(fill = group), width = 0.48, outlier.shape = NA, alpha = 0.13, linewidth = 0.55) +
    geom_point(aes(shape = group), size = 2.5,
      position = position_jitter(width = 0.075, height = 0, seed = SEED)) +
    annotate("text", x = 1.5, y = r[2] + 0.25 * span,
      label = sprintf("log2FC = %.2f\nFDR = %s", e$log2FC_display, format_fdr(e$FDR)),
      family = "Arial", size = text_size, lineheight = 1.12) +
    scale_colour_manual(values = group_colors) + scale_fill_manual(values = group_colors) +
    scale_shape_manual(values = c(Young = 16, Old = 17)) +
    scale_y_continuous(limits = c(r[1] - 0.10 * span, r[2] + 0.43 * span),
      breaks = pretty(r, n = 4), expand = expansion(mult = 0)) +
    labs(title = assay_name, x = NULL, y = y_label) +
    theme_pub() + theme(legend.position = "none")
}

#### Fig 1g #####
# 12只小鼠的Dpp3/DPP3个体表达；箱体为IQR、中线为中位数、须为1.5×IQR。
rna_plot <- dpp3_values[dpp3_values$assay == "Transcriptome", ]
protein_plot <- dpp3_values[dpp3_values$assay == "Proteome", ]
p_rna <- expression_panel(rna_plot, "Transcriptome",
  expression(paste("Dpp3 log"[displaystyle(2)], "(normalized count + 1)")))
p_protein <- expression_panel(protein_plot, "Proteome",
  expression(paste("DPP3 log"[displaystyle(2)], " normalized protein abundance")))
p_expression <- patchwork::wrap_plots(p_rna, p_protein, nrow = 1, widths = c(1, 1))

width_cm <- 18; height_cm <- 7
png_file <- file.path(OUTPUT_DIR, paste0("Figure1g_DPP3_RNA蛋白个体表达", VERSION, ".png"))
pdf_file <- file.path(OUTPUT_DIR, paste0("Figure1g_DPP3_RNA蛋白个体表达", VERSION, ".pdf"))
if (file.exists(png_file) || file.exists(pdf_file)) stop("v6表达图已存在，请递增VERSION。")
ggplot2::ggsave(filename = png_file, plot = p_expression, width = width_cm, height = height_cm,
  units = "cm", dpi = 600, bg = "transparent", scale = 1)
cairo_pdf(filename = pdf_file, width = width_cm / 2.54, height = height_cm / 2.54,
  family = "Arial", pointsize = font_pt, bg = "transparent")
print(p_expression)
dev.off()

#### Fig 1 aging-burden effect #####
# Composite cartilage-aging burden：Hedges g及5000次bootstrap 95%CI；P为924种分组的双侧exact permutation P。
burden_effect <- module_effects[module_effects$component == "Composite", ]
burden_effect$y <- ifelse(burden_effect$assay == "Transcriptome", 2, 1)
burden_effect$P_label <- sprintf("Exact P = %.3f", burden_effect$exact_permutation_P)
p_burden_effect <- ggplot(burden_effect, aes(Hedges_g, y)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.45, colour = "#808080") +
  geom_segment(aes(x = g_CI_low, xend = g_CI_high, yend = y), linewidth = 0.8, colour = "#8E8E8E") +
  geom_point(size = 3.2, colour = "#D95A59") +
  geom_text(aes(x = 11.3, y = y + 0.20, label = P_label), hjust = 1,
    family = "Arial", size = text_size) +
  scale_x_continuous(limits = c(-0.4, 11.7), breaks = c(0, 4, 8), expand = expansion(mult = 0)) +
  scale_y_continuous(breaks = c(2, 1), labels = c("Transcriptome", "Proteome"),
    expand = expansion(mult = c(0.20, 0.23))) +
  labs(x = "Cartilage-aging burden\nHedges' g (Old - Young)\n95% bootstrap CI", y = NULL) +
  theme_pub() + theme(axis.ticks.y = element_blank())

width_cm <- 8.5; height_cm <- 7.5
png_file <- file.path(OUTPUT_DIR, paste0("Figure1_aging_burden_跨组学效应量", VERSION, ".png"))
pdf_file <- file.path(OUTPUT_DIR, paste0("Figure1_aging_burden_跨组学效应量", VERSION, ".pdf"))
if (file.exists(png_file) || file.exists(pdf_file)) stop("v6 aging-burden统计图已存在，请递增VERSION。")
ggplot2::ggsave(filename = png_file, plot = p_burden_effect, width = width_cm, height = height_cm,
  units = "cm", dpi = 600, bg = "transparent", scale = 1)
cairo_pdf(filename = pdf_file, width = width_cm / 2.54, height = height_cm / 2.54,
  family = "Arial", pointsize = font_pt, bg = "transparent")
print(p_burden_effect)
dev.off()

#### Fig 1 DPP3 effect #####
# RNA为DESeq2+ashr后验均值及95% CrI；蛋白为limma log2FC及95%CI。
rna_de <- analysis$rna_de
valid <- is.finite(rna_de$log2FC) & is.finite(rna_de$lfcSE) & rna_de$lfcSE > 0
ash_fit <- ashr::ash(rna_de$log2FC[valid], rna_de$lfcSE[valid], mixcompdist = "normal", method = "shrink", outputlevel = 2)
rna_index <- which(tolower(rna_de$gene[valid]) == "dpp3")
rna_ci_all <- ashr::ashci(ash_fit, betaindex = rna_index)
protein_de <- analysis$protein_de[tolower(analysis$protein_de$gene) == "dpp3", ]
dpp3_effect_plot <- data.frame(
  assay = c("Transcriptome", "Proteome"), y = c(2, 1),
  estimate = c(ashr::get_pm(ash_fit)[rna_index], protein_de$logFC),
  lower = c(rna_ci_all[rna_index, 1], protein_de$CI.L),
  upper = c(rna_ci_all[rna_index, 2], protein_de$CI.R),
  interval = c("95% CrI", "95% CI"), FDR = dpp3_effects$FDR)
stopifnot(abs(dpp3_effect_plot$estimate[1] - dpp3_effects$log2FC_display[1]) < 1e-8)
dpp3_effect_plot$label <- paste0(dpp3_effect_plot$interval, "\nFDR = ", format_fdr(dpp3_effect_plot$FDR))
p_dpp3_effect <- ggplot(dpp3_effect_plot, aes(estimate, y)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.45, colour = "#808080") +
  geom_segment(aes(x = lower, xend = upper, yend = y), linewidth = 0.8, colour = "#8E8E8E") +
  geom_point(shape = 21, size = 3.2, stroke = 0.55, fill = "#E69F00", colour = "black") +
  geom_text(aes(x = 0.77, y = y + 0.20, label = label), hjust = 1,
    family = "Arial", size = text_size, lineheight = 1.12) +
  scale_x_continuous(limits = c(-0.04, 0.80), breaks = c(0, 0.2, 0.4), expand = expansion(mult = 0)) +
  scale_y_continuous(breaks = c(2, 1), labels = c("Transcriptome", "Proteome"),
    expand = expansion(mult = c(0.20, 0.28))) +
  labs(x = expression(atop(paste("DPP3 log"[displaystyle(2)], " fold change"), "(Old / Young)")), y = NULL) +
  theme_pub() + theme(axis.ticks.y = element_blank())

width_cm <- 8.5; height_cm <- 7.5
png_file <- file.path(OUTPUT_DIR, paste0("Figure1_Dpp3_跨组学效应量", VERSION, ".png"))
pdf_file <- file.path(OUTPUT_DIR, paste0("Figure1_Dpp3_跨组学效应量", VERSION, ".pdf"))
if (file.exists(png_file) || file.exists(pdf_file)) stop("v6 DPP3统计图已存在，请递增VERSION。")
ggplot2::ggsave(filename = png_file, plot = p_dpp3_effect, width = width_cm, height = height_cm,
  units = "cm", dpi = 600, bg = "transparent", scale = 1)
cairo_pdf(filename = pdf_file, width = width_cm / 2.54, height = height_cm / 2.54,
  family = "Arial", pointsize = font_pt, bg = "transparent")
print(p_dpp3_effect)
dev.off()

# 最小导出检查：PNG必须具有透明通道、600dpi，并符合指定厘米尺寸。
png_files <- file.path(OUTPUT_DIR, paste0(c("Figure1g_DPP3_RNA蛋白个体表达",
  "Figure1_aging_burden_跨组学效应量", "Figure1_Dpp3_跨组学效应量"), VERSION, ".png"))
expected_cm <- rbind(c(18, 7), c(8.5, 7.5), c(8.5, 7.5))
for (i in seq_along(png_files)) {
  z <- png::readPNG(png_files[i], info = TRUE); info <- attr(z, "info")
  stopifnot(length(dim(z)) == 3L, dim(z)[3] == 4L, all(abs(info$dpi - 600) < 0.1),
    abs(dim(z)[2] / info$dpi[1] * 2.54 - expected_cm[i, 1]) < 0.01,
    abs(dim(z)[1] / info$dpi[2] * 2.54 - expected_cm[i, 2]) < 0.01,
    max(z[c(1, dim(z)[1]), c(1, dim(z)[2]), 4]) == 0)
}
message("完成：表达图18×7 cm；两张效应量图各8.5×7.5 cm；PDF和透明PNG已输出。")
