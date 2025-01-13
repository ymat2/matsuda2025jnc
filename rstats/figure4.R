library(conflicted)
library(tidyverse)
library(phytools)
library(ggtree)
library(OptM)
library(myrrr)


## NJtree

accession2name = readr::read_tsv("data/sra_accession.tsv") |>
  dplyr::select(Run, sample_id) |>
  dplyr::rename(label = Run)

jpool_names = readxl::read_excel(path = "data/Japanese_Chicken_DNAseq.xlsx", sheet = "Sheet1") |>
  dplyr::select(sample, sample_name)

rapidnj = ape::read.tree("out/phylo/popstr.sub.nj.tree")
tip2lab = dplyr::tibble(label = rapidnj$tip.label) |>
  dplyr::mutate(label = stringr::str_remove_all(label, "'")) |>
  dplyr::left_join(jpool_names, by = dplyr::join_by(label == sample)) |>
  dplyr::left_join(accession2name, by = "label") |>
  dplyr::mutate(sample_id = dplyr::if_else(is.na(sample_id), sample_name, sample_id))
rapidnj$tip.label = tip2lab$sample_id

jpool_tips = tip2lab |> dplyr::filter(stringr::str_detect(label, "P$"))
rapidnj = ape::keep.tip(rapidnj, jpool_tips$sample_name)

grp = dplyr::tibble(
  nodes = c(43, 53, 66),
  names = c("Japan1\n(Shokoku group)", "Japan2\n(Gamecocks)", "Japan3\n(Others)"),
  colors = palette.colors(palette = "Paired")[1:3]
)
  
pnj = ggtree(rapidnj) + 
  geom_tiplab(size = 3) +
  geom_cladelab(
    node = grp$nodes, 
    label = grp$names, 
    barcolor = grp$colors, 
    textcolor = grp$colors, 
    fontsize = 4,
    barsize = 1, 
    offset = .05, 
    hjust = -.2, 
    align = TRUE
    ) + 
  scale_x_continuous(expand = expansion(mult = c(.1, .75)))
pnj

### 

folder = "out/treemix"
prefix = stringr::str_c(folder, "treemix", sep  ="/")


## Decide optimal M using OptM -------------------------------------------------

test.optM = OptM::optM(folder, method = "Evanno")

deltaM = ggplot(test.optM) +
  aes(x = m, y = Deltam) +
  geom_point(color = "#444444", size = 3) +
  geom_line(color = "#444444", linewidth = 1) +
  scale_x_continuous(labels = seq(0, 8), breaks = seq(0, 8)) +
  labs(
    x = expression("m (migration edges)"),
    y = expression(paste({Delta}, "m", sep=""))
  ) +
  theme_classic(base_size = 14) 

deltaM

## Draw Tree of optimal M ------------------------------------------------------

opt_M = 5
stem = paste0(prefix, ".1.", opt_M)
tmobj = myrrr::read_treemix(stem = stem)

new_pop = c("Japan2 (Gamecocks)", "Shandong", "Shanxi", "Gamecock (China)", 
            "Yunnan", "Gamecock (Thailand, Vietnam)", "Vietnam", "Guangxi", "Jiangsu", "RJF",
            "Japan1 (Shokoku group)", "Jiangxi", "South Korea", "Henan",
            "Beijing", "Thailand", "Japan3 (Others)")
tmobj$layout$tips$pop = new_pop

ptm = myrrr::plot_treemix(tmobj) +
  scale_x_continuous(expand = expansion(mult = c(.05, .2))) +
  theme_treemix(base_size = 14) +
  theme(
    legend.position = "inside",
    legend.justification = c(1, 1),
    legend.background = element_blank()
  )
ptm

p1 = cowplot::plot_grid(pnj, deltaM, nrow = 2, rel_heights = c(2, 1),
                        labels = c("a", "b"), label_size = 20)
p = cowplot::plot_grid(p1, ptm, ncol = 2, rel_widths = c(2, 3), 
                       labels = c("", "c"), label_size = 20)
ggsave("images/treemix_result.tif", p, h = 9, w = 9, bg = "#FFFFFF")
