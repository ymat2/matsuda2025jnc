library(conflicted)
library(tidyverse)
library(phytools)
library(ggtree)
library(OptM)
library(myrrr)

source("./rstats/common_settings.R")


## NJtree

rapidnj = ape::read.tree("out/phylo/popstr.sub.nj.tree")
tip2lab = dplyr::tibble(label = rapidnj$tip.label) |>
  dplyr::mutate(label = stringr::str_remove_all(label, "'")) |>
  dplyr::left_join(sample2group, by = dplyr::join_by(label == Run)) |>
  dplyr::mutate(label = sample_id)
rapidnj$tip.label = tip2lab$label

tip2grp = tip2lab |> dplyr::select(label, group)

jnode = data.frame(
  nodes = c(234, 212, 218, 222, 191),
  nodename = c("1", "2", "3", "4", "5")
)

pnj = rapidnj |>
  phytools::reroot(node.number = 151, position = .0005) |>
  dplyr::full_join(tip2grp, by = "label") |>
  ggtree(aes(color = group), linewidth = .33) +
  geom_hline(yintercept = c(40.5, 52.5, 58.5, 63.5, 76.5, 96.5), linetype = "dashed", color = "#333333", linewidth = .25) +
  geom_tiplab(size = 2, align = TRUE) +
  geom_nodelab(aes(label = label), size = 1.5, hjust = -.1) +
  geom_cladelab(node = 234, label = "1", geom = "label", textcolour = "#FFFFFF", fill = "#333333", size = 5, label.r = unit(.5, "lines"), label.size = 0, offset = .02, barsize = 0, align = TRUE) +
  geom_cladelab(node = 212, label = "2", geom = "label", textcolour = "#FFFFFF", fill = "#333333", size = 5, label.r = unit(.5, "lines"), label.size = 0, offset = .02, barsize = 0, align = TRUE) +
  geom_cladelab(node = 218, label = "3", geom = "label", textcolour = "#FFFFFF", fill = "#333333", size = 5, label.r = unit(.5, "lines"), label.size = 0, offset = .02, barsize = 0, align = TRUE) +
  geom_cladelab(node = 222, label = "4", geom = "label", textcolour = "#FFFFFF", fill = "#333333", size = 5, label.r = unit(.5, "lines"), label.size = 0, offset = .02, barsize = 0, align = TRUE) +
  geom_cladelab(node = 191, label = "5", geom = "label", textcolour = "#FFFFFF", fill = "#333333", size = 5, label.r = unit(.5, "lines"), label.size = 0, offset = .02, barsize = 0, align = TRUE) +
  geom_treescale(x = .005, y = 131, fontsize = 3, linesize = 1, offset = 1) +
  scale_x_continuous(expand = expansion(mult = c(.05, .1))) +
  scale_y_continuous(expand = expansion(mult = c(.05, .05))) +
  scale_color_manual(values = color_scale) +
  theme(legend.position = "none")
pnj

###

ja_grp = dplyr::tibble(
  x = c(2, 2, 2, 2, 2),
  y = seq(1, 5),
  label = c("5", "4", "3", "2", "1"),
  txt = c(
    "Commercial breeds",
    "Other breeds",
    "Shokoku group",
    "Nongame breeds",
    "Gamecocks (Shamo)"
  )) |>
  ggplot() +
  aes(y = y) +
  geom_label(x = 1, aes(label = label), color = "#FFFFFF", fill = "#333333", size = 5, label.r = unit(.5, "lines"), label.size = 0) +
  geom_text(aes(x = x, label = txt), hjust = 0, color ="#333333", size = 5) +
  #ggplot2::annotate("segment", x = 1.5, y = 2, xend = 2, yend = 3, color = "#333333") +
  #ggplot2::annotate("segment", x = 1.5, y = 3, xend = 2.2, yend = 3, color = "#333333") +
  #ggplot2::annotate("segment", x = 1.5, y = 4, xend = 2, yend = 3, color = "#333333") +
  ggplot2::annotate("point", x = 3, y = 0.25, shape = 25, size = 3, color = "#FFFFFF", fill = "darkred") +
  ggplot2::annotate("text", x = 1, y = -.5, label = "Group for TreeMix", size = 6, fontface = "bold", color ="darkred", hjust = 0) +
  scale_x_continuous(limits = c(0, 10)) +
  scale_y_continuous(limits = c(-1, 6)) +
  theme_void()

ja_grp

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

opt_M = 1
stem = paste0(prefix, ".1.", opt_M)
tmobj = myrrr::read_treemix(stem = stem)

new_pop = c(
  "Japan (Gamecocks)", "South korea", "RJF", "Vietnam", "Thailand", "Yunnan", 
  "Gamecock (China)", "Commercial, Other regions (Node 5)", "Japan (Shokoku)", "Jiangxi", 
  "Henan", "Japan (Node 4)", "Jiangsu", "Guangxi", "Beijing", "Shandong", 
  "Japan (Node 2)", "Gamecock (Southeast Asia)", "Shanxi"  
)
tmobj$layout$tips$pop = new_pop

ptm = myrrr::plot_treemix(tmobj) +
  scale_x_continuous(expand = expansion(mult = c(.05, .8))) +
  theme_treemix(base_size = 14) +
  theme(
    legend.position = "inside",
    legend.justification = c(1, 0),
    legend.background = element_blank()
  )
ptm

prt = cowplot::plot_grid(ja_grp, deltaM, ncol = 2, rel_widths = c(1, 1),
                         labels = c("b", "c"), label_size = LABELSIZE)
pr = cowplot::plot_grid(prt, ptm, nrow = 2, rel_heights = c(1, 2),
                        labels = c("", "d"), label_size = LABELSIZE)
p = cowplot::plot_grid(pnj, pr, ncol = 2, rel_widths = c(2, 3), 
                       labels = c("a", ""), label_size = LABELSIZE)
ggsave("images/treemix_result_revise.jpg", p, h = 10, w = 10, bg = "#FFFFFF")

