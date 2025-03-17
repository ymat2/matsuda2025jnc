library(conflicted)
library(tidyverse)

source("./rstats/common_settings.R")


##### PCA before removing batch effect #####

eigenvec = readr::read_delim("out/batch/popstr.pca.before.eigenvec") |>
  dplyr::left_join(sample2group, by = dplyr::join_by(IID == Run)) |>
  dplyr::mutate(
    Group = forcats::fct_relevel(group, names(color_scale)), 
    Sequencer = system) |>
  tidyr::drop_na()

eigenval = readr::read_csv("out/batch/popstr.pca.before.eigenval", col_names = "V1")
df = data.frame(pc = 1:nrow(eigenval), eigenval/sum(eigenval)*100)

eigenvec |>
  dplyr::mutate(SEQ = dplyr::if_else(Sequencer == "NovaSeq", 1, 2)) |>
  dplyr::select(`#FID`, IID, SEQ) |>
  readr::write_tsv("out/batch/sequencer.pheno")

fa = ggplot(eigenvec) +
  aes(PC1, PC2, color = Group) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_color_manual(values = color_scale) +
  #coord_fixed() +
  xlab(paste0("PC1 (", round(df[1,2], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(df[2,2], digits = 2), "%)")) +
  labs(tag = "a") +
  theme_test(base_size = 16) +
  theme(plot.tag.location = "plot")
fa

fb = ggplot(eigenvec) +
  aes(PC1, PC2, color = Sequencer) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_color_discrete(type = palette.colors(palette = "Set1")) +
  #coord_fixed() +
  xlab(paste0("PC1 (", round(df[1,2], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(df[2,2], digits = 2), "%)")) +
  theme_test(base_size = 16) +
  theme(plot.tag.location = "plot")
fb


##### Check distribution of GWAS p-value and Fst #####

gwas_result = readr::read_tsv("out/batch/gwas-result.txt") |>
  dplyr::mutate(extract = dplyr::if_else(FDR_BH > 0.05, "Y", "N"))

nall = nrow(gwas_result)
nexc = nrow(gwas_result |> dplyr::filter(extract == "N"))
pexc = round((nexc/nall)*100, digits = 2)

txt = stringr::str_c(
  "Number of removed SNPs:\n",
  nexc |> format(big.mark=",", scientific = FALSE),
  "/",
  nall |> format(big.mark=",", scientific = FALSE),
  " (",
  pexc,
  "%)"
)

fc = ggplot(gwas_result) +
  geom_histogram(aes(x = FDR_BH, fill = extract), binwidth = 0.01, boundary = 0) +
  geom_vline(xintercept = 0.05, linetype = "longdash", color = "#BBBBBB") +
  ggplot2::annotate("text", label = txt, x = 0, y = 250000, hjust = 0, size = 5) +
  scale_x_continuous(breaks = c(0.05), labels = c("p-adjusted = 0.05")) +
  scale_y_continuous(limit = c(0, 300000), labels = scales::label_number(scale_cut = scales::cut_long_scale())) +
  labs(x = "adjusted p-value of GWAS", y = "", tag = "b") +
  scale_fill_manual(values = c("#4d4d4d", "#BBBBBB")) +
  theme_test(base_size = 16) +
  theme(legend.position = "none", plot.tag.location = "plot")
fc


##### PCA after removing batch effect #####

eigenvec2 = readr::read_delim("out/batch/popstr.pca.after.eigenvec") |>
  dplyr::left_join(sample2group, by = dplyr::join_by(IID == Run)) |>
  dplyr::mutate(
    Group = forcats::fct_relevel(group, names(color_scale)), 
    Sequencer = system) |>
  tidyr::drop_na()

eigenval2 = readr::read_csv("out/batch/popstr.pca.after.eigenval", col_names = "V1")
df2 = data.frame(pc = 1:nrow(eigenval2), eigenval2/sum(eigenval2)*100)

fd = ggplot(eigenvec2) +
  aes(PC1, PC2, color = Group) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_color_manual(values = color_scale) +
  #coord_fixed() +
  xlab(paste0("PC1 (", round(df2[1,2], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(df2[2,2], digits = 2), "%)")) +
  labs(tags = "c") +
  theme_test(base_size = 16) +
  theme(plot.tag.location = "plot")
fd

fe = ggplot(eigenvec2) +
  aes(PC1, PC2, color = Sequencer) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_color_discrete(type = palette.colors(palette = "Set1")) +
  #coord_fixed() +
  xlab(paste0("PC1 (", round(df2[1,2], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(df2[2,2], digits = 2), "%)")) +
  theme_test(base_size = 16)
fe


## cowplot

library(patchwork)

row1 = fa + fb + patchwork::plot_layout(guides = 'collect')
row2 = fd + fe + patchwork::plot_layout(guides = 'collect')
p = row1 / fc / row2 +
  patchwork::plot_layout(nrow = 3, heights = c(1, .7, 1)) &
  theme(plot.tag = element_text(face = "bold", size = 22))
p
ggsave("images/figS1_revise.png", p, w = 10, h = 10, bg = "#FFFFFF")
