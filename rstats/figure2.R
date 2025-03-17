library(conflicted)
library(tidyverse)
library(phytools)
library(ggtree)
library(patchwork)

source("./rstats/common_settings.R")


### NJtree of mitchondrial whole genome -----

nj = ape::read.tree("out/phylo/popstr.mt.nj.tree")
tip2sym = dplyr::tibble(label = stringr::str_split(nj$tip.label, "/", simplify = TRUE)[,2]) |>
  dplyr::left_join(sample2group, by = dplyr::join_by(label == Run))
nj$tip.label = tip2sym$sample_id

njtr = ggtree(nj, layout = "circular", branch.length = "none", linewidth = .5) + 
  geom_tiplab(size = 1.5, offset = .2)

### Haplogroup classification by MitoToolPy -----

mitotool = readr::read_tsv("out/phylo/popstr.mt.mitotool", skip = 4, comment = "---") |>
  dplyr::mutate(label = stringr::str_split(Sample_Name, "/", simplify = TRUE)[,2]) |>
  dplyr::left_join(sample2group, by = dplyr::join_by(label == Run)) |>
  dplyr::rename("Haplogroup (MitoToolPy)" = Haplogroup) |>
  dplyr::select(sample_id, `Haplogroup (MitoToolPy)`) |>
  tibble::column_to_rownames(var = "sample_id")

df_group = tip2sym |>
  dplyr::select(sample_id, group) |>
  tibble::column_to_rownames(var = "sample_id") |>
  dplyr::mutate(group = forcats::fct_relevel(group, names(color_scale)))

df_check = mitotool |> 
  tibble::rownames_to_column(var = "sample_id") |> 
  dplyr::left_join(sample2group, by = "sample_id")
  

pnj = gheatmap(njtr, data = mitotool, offset = 8, width = .1, colnames = FALSE) +
  scale_fill_viridis_d(option = "H") +
  labs(fill = "Haplogroup (MitoToolPy)") +
  guides(fill = guide_legend(ncol = 2))

pnj = pnj + ggnewscale::new_scale_fill()
pnj2 = gheatmap(pnj, data = df_group, offset = 6, width = .1, colnames = FALSE) +
  scale_fill_manual(values = color_scale) +
  labs(fill = "Group") +
  guides(fill = guide_legend(ncol = 2))
pnj2


### NJtree for Japanese breeds -----

jpool_tips = tip2sym |> dplyr::filter(group == "Japan")
njjp = nj |> ape::keep.tip(jpool_tips$sample_id)

### Use aplot

mitotool_jp = mitotool |>
  tibble::rownames_to_column(var = "sample_id") |>
  dplyr::filter(sample_id %in% jpool_tips$sample_id) |>
  tidyr::separate_rows(`Haplogroup (MitoToolPy)`, sep = ",\\s*") |>
  dplyr::mutate(obs = "1") |>
  dplyr::arrange(`Haplogroup (MitoToolPy)`) |>
  tidyr::pivot_wider(names_from = `Haplogroup (MitoToolPy)`, values_from = obs, values_fill = "0") |>
  tidyr::pivot_longer(-1, names_to = "Haplogroup", values_to = "obs") |>
  dplyr::mutate(Haplogroup = Haplogroup |> forcats::as_factor() |> forcats::fct_rev())

pnjjp = ggtree(njjp, branch.length = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0))) + 
  layout_dendrogram() +
  theme()

haplo_dot = ggplot(mitotool_jp|> dplyr::filter(obs == "1")) +
  aes(x = sample_id, y = Haplogroup) +
  geom_point(aes(color = obs), size = 2.5) +
  geom_line(data = mitotool_jp |> dplyr::filter(obs == "1"), 
            aes(group = sample_id), color = "#444444", linewidth = 1) +
  scale_color_manual(values = c("0" = "#DDDDDD", "1" = "#444444")) +
  theme_minimal() +
  theme(
    #panel.background = element_blank(),
    legend.position = "none",
    #panel.grid.major.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
)

pp = haplo_dot |> aplot::insert_top(pnjjp, height = .3) |> aplot::as.patchwork()

### arrange plots -----

p = cowplot::plot_grid(
  pnj2, pp, 
  nrow = 2, labels = c("a", "b"), rel_heights = c(2, 1), 
  label_size = LABELSIZE)
p
ggsave(file = "images/mitochondria_haplogroup_revise.jpg", p, w = 10, h = 10, bg = "#FFFFFF")

