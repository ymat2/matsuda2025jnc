library(conflicted)
library(tidyverse)
library(phytools)
library(ggtree)
library(patchwork)


accession2name = readr::read_tsv("data/sra_accession.tsv") |>
  dplyr::select(Run, sample_id, country) |>
  dplyr::rename(label = Run)

jpool_names = readxl::read_excel(path = "data/Japanese_Chicken_DNAseq.xlsx", sheet = "Sheet1") |>
  dplyr::select(sample, sample_name)

### NJ -----

nj = ape::read.tree("out/phylo/popstr.mt.nj.tree")
tip2sym = data.frame(label = stringr::str_split(nj$tip.label, "/", simplify = TRUE)[,2]) |>
  dplyr::left_join(accession2name, by = "label") |>
  dplyr::left_join(jpool_names, by = dplyr::join_by(label == sample)) |>
  dplyr::mutate(
    sample_id = dplyr::if_else(is.na(sample_name), sample_id, sample_name),
    country = dplyr::if_else(is.na(country), "Japan", country)
)
nj$tip.label = tip2sym$sample_id

njtr = ggtree(nj, layout = "circular", branch.length = "none", linewidth = .5) + 
  geom_tiplab(size = 1.5, offset = .2)

### ML -----

ml = ape::read.tree("out/phylo/popstr.mt.treefile")
tip2sym = data.frame(label = stringr::str_split(ml$tip.label, "/", simplify = TRUE)[,2]) |>
  dplyr::left_join(accession2name, by = "label") |>
  dplyr::left_join(jpool_names, by = dplyr::join_by(label == sample)) |>
  dplyr::mutate(
    sample_id = dplyr::if_else(is.na(sample_name), sample_id, sample_name),
    country = dplyr::if_else(is.na(country), "Japan", country)
  )
ml$tip.label = tip2sym$sample_id

mltr = ggtree(ml, layout = "circular", branch.length = "none", linewidth = .5) + 
  geom_tiplab(size = 1.5, offset = .2)

### Haplogroup classification by MitoToolPy -----

mitotool = readr::read_tsv("out/phylo/popstr.mt.mitotool", skip = 4, comment = "---") |>
  dplyr::mutate(label = stringr::str_split(Sample_Name, "/", simplify = TRUE)[,2]) |>
  dplyr::left_join(accession2name, by = "label") |>
  dplyr::left_join(jpool_names, by = dplyr::join_by(label == sample)) |>
  dplyr::mutate(sample_id = dplyr::if_else(is.na(sample_id), sample_name, sample_id)) |>
  dplyr::rename("Haplogroup (MitoToolPy)" = Haplogroup) |>
  dplyr::select(sample_id, `Haplogroup (MitoToolPy)`) |>
  tibble::column_to_rownames(var = "sample_id")

df_country = tip2sym |>
  dplyr::select(sample_id, country) |>
  tibble::column_to_rownames(var = "sample_id")

df_check = mitotool |> 
  tibble::rownames_to_column(var = "sample_id") |> 
  dplyr::left_join(accession2name, by = "sample_id")
  

pnj = gheatmap(njtr, data = mitotool, offset = 8, width = .1, colnames = FALSE) +
  scale_fill_viridis_d(option = "H") +
  labs(fill = "Haplogroup (MitoToolPy)") +
  guides(fill = guide_legend(ncol = 2))

pnj = pnj + ggnewscale::new_scale_fill()
pnj2 = gheatmap(pnj, data = df_country, offset = 6, width = .1, colnames = FALSE) +
  scale_fill_discrete(type = palette.colors(palette = "Okabe-Ito")[-1]) +
  labs(fill = "Country") +
  guides(fill = guide_legend(ncol = 2))

pml = gheatmap(mltr, data = mitotool, offset = 8, width = .1, colnames = FALSE) +
  scale_fill_viridis_d(option = "H") +
  labs(fill = "Haplogroup (MitoToolPy)")

p = pnj + theme(legend.position = "none") + pml +
  patchwork::plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 24))

ggsave(file = "images/haplogroups_by_mitotool.png", p, w = 13.5, h = 6)


### MLtree for Japanese breeds -----

jpool_tips = tip2sym |> dplyr::filter(country == "Japan")
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

p = cowplot::plot_grid(pnj2, pp, nrow = 2, labels = c("a", "b"), rel_heights = c(3, 2), label_size = 20)
ggsave(file = "images/mitochondria_haplogroup.tif", p, w = 9, h = 9, bg = "#FFFFFF")

