library(conflicted)
library(tidyverse)
library(cowplot)

source("./rstats/common_settings.R")


## Photos -----

ong = cowplot::ggdraw() +
  cowplot::draw_image("photos/ONG.jpg", scale = .95) +
  ggtitle("Onagadori\n(ONG)") +
  theme_test(base_size = BASESIZE) +
  theme(
    plot.title = element_text(hjust = .5, margin = margin(0,0,-5,0)), 
    panel.border = element_blank()
  )

koe = cowplot::ggdraw() +
  cowplot::draw_image("photos/KOE.jpg", scale = .95) +
  ggtitle("Koeyoshi\n(KOE)") +
  theme_test(base_size = BASESIZE) +
  theme(
    plot.title = element_text(hjust = .5, margin = margin(0,0,-5,0)), 
    panel.border = element_blank()
  )

tmr = cowplot::ggdraw() +
  cowplot::draw_image("photos/TMR.jpg", scale = .95) +
  ggtitle("Tomaru\n(TMR)") +
  theme_test(base_size = BASESIZE) +
  theme(
    plot.title = element_text(hjust = .5, margin = margin(0,0,-5,0)), 
    panel.border = element_blank()
  )

tot = cowplot::ggdraw() +
  cowplot::draw_image("photos/TOT1.jpg", scale = .95) +
  ggtitle("Totenko\n(TOT)") + 
  theme_test(base_size = BASESIZE) +
  theme(
    plot.title = element_text(hjust = .5, margin = margin(0,0,-5,0)), 
    panel.border = element_blank()
  )

osm = cowplot::ggdraw() +
  cowplot::draw_image("photos/OSM.jpg", scale = .95) +
  ggtitle("Oshamo\n(OSM)") +
  theme_test(base_size = BASESIZE) +
  theme(
    plot.title = element_text(hjust = .5, margin = margin(0,0,-5,0)), 
    panel.border = element_blank()
  )

shk = cowplot::ggdraw() +
  cowplot::draw_image("photos/SHK.jpg", scale = .95) +
  ggtitle("Shokoku\n(SHK)") +
  theme_test(base_size = BASESIZE) +
  theme(
    plot.title = element_text(hjust = .5, margin = margin(0,0,-5,0)), 
    panel.border = element_blank()
  )

chb = cowplot::ggdraw() +
  cowplot::draw_image("photos/CHB.jpg", scale = .95) +
  ggtitle("Chabo\n(CHB)") +
  theme_test(base_size = BASESIZE) +
  theme(
    plot.title = element_text(hjust = .5, margin = margin(0,0,-5,0)), 
    panel.border = element_blank()
  )

p1 = cowplot::plot_grid(ong, shk, osm, chb, nrow = 1, 
                        rel_widths = c(2.46, 1.65, 1.1, 1.65))

## Sample location map -----

world_map = ggplot2::map_data("world") |> dplyr::as_tibble()

file_path = "data/Japanese_Chicken_DNAseq.xlsx"
sheet = "Sheet1"
sample_info = readxl::read_excel(path = file_path, sheet = sheet) |>
  dplyr::filter(poolseq > 1) |>
  dplyr::select(sample, sample_name, common_name, poolseq, sampling_location, long, lat)

jmap = sample_info |>
  dplyr::filter(long != "NA") |>
  dplyr::mutate(long = as.numeric(long), lat = as.numeric(lat)) |>
  ggplot() + 
  aes(long, lat) +
  geom_map(data = world_map, aes(map_id = region), map = world_map, fill = "#bbbbbb", linewidth = 0) +
  geom_point(size = 3, position = "jitter", shape = 21, color = "#333333", fill = "#FFFFFF") +
  scale_x_continuous(limits = c(125, 150), labels = scales::label_number(suffix = "°E")) +
  scale_y_continuous(limits = c(25, 50), labels = scales::label_number(suffix = "°N")) +
  scale_color_brewer(palette = "Set1") +
  coord_fixed() +
  theme_bw(base_size = 20) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    legend.title = element_blank(),
    legend.position = "inside", 
    legend.justification = c(1, 0),
    legend.background = element_blank()
  )

jpool_snp_count = readr::read_tsv("../jpool/out/snp_count.txt") 
abrc_snp_count = readr::read_tsv("../ABRC/out/snp_count.txt") 
snp_count = dplyr::bind_rows(jpool_snp_count, abrc_snp_count) |>
  dplyr::left_join(sample_info, by = "sample") |>
  tidyr::pivot_longer(dplyr::starts_with("n_"), names_to = "type", values_to = "count") |>
  tidyr::drop_na() |>
  dplyr::mutate(sample_name = forcats::fct_reorder(sample_name, count, .fun = sum)) |>
  ggplot() +
  aes(x = sample_name, y = count) +
  geom_bar(aes(fill = type), stat = "identity") +
  scale_y_continuous(
    labels = scales::label_number(scale_cut = scales::cut_long_scale()),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_fill_brewer(
    palette = "Paired",
    labels = c("#heteroSNPs", "#homoSNPs", "#indels")
  ) +
  labs(y = "Variant count") +
  coord_flip() +
  theme_bw(base_size = 16) +
  theme(
    axis.text.y = element_text(size = 8),
    axis.title.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "inside",
    legend.justification = c(1, 0),
    legend.title = element_blank(),
    legend.background = element_blank()
  )

p2 = cowplot::plot_grid(jmap, NULL, snp_count, ncol = 3, rel_widths = c(1.4, .05, 1), scale = .95,
                        labels = c("b", "c", ""), label_size = LABELSIZE, align = "h", axis = "tb")

p = cowplot::plot_grid(p1, p2, nrow = 2, rel_heights = c(1, 2), labels = c("a", ""), label_size = LABELSIZE)
ggsave(file = "images/sample_information_revise.jpg", p, h = 8, w = 12, bg = "#FFFFFF")

