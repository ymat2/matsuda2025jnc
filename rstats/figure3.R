library(conflicted)
library(tidyverse)
library(phytools)
library(ggtree)
library(patchwork)

source("./rstats/common_settings.R")


### PCA

eigenvec2 = readr::read_delim("out/batch/popstr.pca.after.eigenvec") |>
  dplyr::left_join(sample2group, by = dplyr::join_by(IID == Run)) |>
  dplyr::mutate(
    Group = forcats::fct_relevel(group, names(color_scale)), 
    Sequencer = system) |>
  tidyr::drop_na()

eigenval2 = readr::read_csv("out/batch/popstr.pca.after.eigenval", col_names = "V1")
df2 = data.frame(pc = 1:nrow(eigenval2), eigenval2/sum(eigenval2)*100)


fa = ggplot(eigenvec2) +
  aes(PC1, PC2, color = Group) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = color_scale) +
  coord_fixed() +
  xlab(paste0("PC1 (", round(df2[1,2], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(df2[2,2], digits = 2), "%)")) +
  theme_test(base_size = BASESIZE) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  guides(color = guide_legend(nrow = 2))
fa

fb = ggplot(eigenvec2) +
  aes(PC3, PC4, color = Group) +
  geom_point(size = 2, alpha = 0.9) +
  scale_color_manual(values = color_scale) +
  coord_fixed() +
  xlab(paste0("PC3 (", round(df2[3,2], digits = 2), "%)")) +
  ylab(paste0("PC4 (", round(df2[4,2], digits = 2), "%)")) +
  theme_test(base_size = BASESIZE) +
  theme(legend.position = "none")
fb


### NJtree

rapidnj = ape::read.tree("out/phylo/popstr.sub.nj.tree")
tip2lab = dplyr::tibble(label = rapidnj$tip.label) |>
  dplyr::mutate(label = stringr::str_remove_all(label, "'")) |>
  dplyr::left_join(sample2group, by = dplyr::join_by(label == Run)) |>
  dplyr::mutate(label = sample_id)
rapidnj$tip.label = tip2lab$label

### Admixture

fam = readr::read_tsv("out/admixture/popstr.snp.filt.fam", col_names = FALSE) |>
  dplyr::rename(Run = X1) |>
  dplyr::left_join(sample2group, by = "Run") |>
  dplyr::select(Run, sample_id, group)

q2 = readr::read_delim("out/admixture/popstr.snp.filt.2.Q", col_names = FALSE) |>
  dplyr::bind_cols(fam) |>
  dplyr::mutate(K = "K = 2") |>
  tidyr::pivot_longer(dplyr::starts_with("X"), names_to = "anc", values_to = "prop")

q3 = readr::read_delim("out/admixture/popstr.snp.filt.3.Q", col_names = FALSE) |>
  dplyr::bind_cols(fam) |>
  dplyr::mutate(K = "K = 3") |>
  tidyr::pivot_longer(dplyr::starts_with("X"), names_to = "anc", values_to = "prop")

q4 = readr::read_delim("out/admixture/popstr.snp.filt.4.Q", col_names = FALSE) |>
  dplyr::bind_cols(fam) |>
  dplyr::mutate(K = "K = 4") |>
  tidyr::pivot_longer(dplyr::starts_with("X"), names_to = "anc", values_to = "prop")

df = dplyr::bind_rows(q2, q3, q4) 

tip2grp = df |>
  dplyr::select(sample_id, group) |>
  dplyr::distinct(sample_id, .keep_all = TRUE) |>
  dplyr::rename(label = sample_id)

pnj2 = rapidnj |>
  dplyr::full_join(tip2grp, by = "label") |>
  ggtree(aes(color = group), linewidth = .33, branch.length = "none") +
  geom_tiplab(size = 2) +
  #geom_nodelab(aes(label = label), size = 1.5, hjust = -.1) +
  scale_x_continuous(expand = expansion(mult = c(.0, .4))) +
  scale_color_manual(values = color_scale) +
  theme(legend.position = "none")
pnj2

adm = ggplot(df) +
  aes(x = sample_id, y = prop, fill = anc) +
  geom_bar(stat = "identity", position = "fill", alpha = .66) +
  scale_fill_brewer(palette = "Set1") +
  facet_grid(cols = vars(K), scales = "free_x", space = "free_x") +
  coord_flip() +
  theme_minimal(base_size = BASESIZE) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 14),
    axis.title = element_blank(),
    axis.text = element_blank()
  )
adm

padm = adm |> aplot::insert_left(pnj2, width = .67) |> aplot::as.patchwork()
padm

cve = readr::read_delim("out/admixture/CV-error.txt",
                        col_names = c("c", "e", "K", "error")) |>
  dplyr::mutate(K = as.integer(stringr::str_extract(K, "[0-9]+"))) |>
  #dplyr::filter(K > 1) |>
  ggplot() +
  aes(x = as.integer(K), y = error) +
  geom_point(size = 2) +
  geom_line(linewidth = 1) +
  scale_x_continuous(breaks = c(5, 10, 15), labels = c(5, 10, 15)) +
  labs(x = "K", y = "CV Error") +
  theme_classic(base_size = BASESIZE)
cve


### Align plots

plb = cowplot::plot_grid(
  fb, cve,
  nrow = 1,
  labels = c("b", "c"), label_size = LABELSIZE, rel_widths = c(1, 1)
)
pl = cowplot::plot_grid(
  fa, plb, 
  nrow = 2, 
  labels = c("a", ""), label_size = LABELSIZE, rel_heights = c(2, 1)
)
p = cowplot::plot_grid(
  pl, padm, 
  ncol = 2,
  labels = c("", "d"), label_size = LABELSIZE, rel_widths = c(1, .5)
)
p
ggsave("images/structure_analysis_revise.jpg", p, w = 10, h = 10, bg = "#FFFFFF")

