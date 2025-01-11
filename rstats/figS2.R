library(conflicted)
library(tidyverse)
library(readxl)


## Maps 

world_map = ggplot2::map_data("world") |> dplyr::as_tibble()

map_asia = ggplot(world_map) + 
  aes(long, lat) +
  geom_map(aes(map_id = region), map = world_map, fill = "#a6bddb", linewidth = 0) +
  scale_x_continuous(limits = c(90, 150), labels = scales::label_number(suffix = "°N")) +
  scale_y_continuous(limits = c(0, 50), labels = scales::label_number(suffix = "°E")) + 
  coord_fixed(ratio = 1) +
  theme_test(base_size = 14) +
  theme(axis.title = element_blank())
map_asia

map_east_asia = world_map |>
  dplyr::filter(region %in% c("China", "South Korea", "North Korea")) |>
  ggplot() + 
  aes(long, lat) +
  geom_map(aes(map_id = region), map = world_map, fill = "#a6bddb", linewidth = 0) +
  coord_fixed(ratio = 1) +
  theme_test(base_size = 14) +
  theme(axis.title = element_blank())
map_east_asia

map_thai = world_map |>
  dplyr::filter(region == c("Thailand")) |>
  ggplot() + 
  aes(long, lat) +
  geom_map(aes(map_id = region), map = world_map, fill = "#a6bddb", linewidth = 0) +
  coord_fixed(ratio = 1) +
  theme_test(base_size = 14) +
  theme(axis.title = element_blank())
map_thai


## Samples

file_path = "data/chicken_dna_sample.xlsx"
sheet = "Sheet1"
samples = readxl::read_excel(file_path, sheet = sheet) |>
  dplyr::select(sample, individuals, long, lat) |>
  tidyr::drop_na() |>
  dplyr::mutate(group = dplyr::case_when(
    stringr::str_detect(sample, "SM|YKD") ~ "Gamecocks",
    stringr::str_detect(sample, "RIR|WKN|WLH|WPR|AKN") ~ "Commercial lines",
    .default = "Native breeds"
  ))

df_sample = dplyr::tibble(
  sample = stringr::str_split(samples$sample, "_", simplify = TRUE)[,1] |>
    stringr::str_remove("P$") |> unique()
  ) |> dplyr::mutate(
    group = dplyr::case_when(
      sample %in% c("HNI", "GJD", "TSA", "MIE", "JTK") ~ "Jidori",
      sample %in% c("ONG", "TOT", "SHK", "MHK", "KKW", "TMR", "KOE", "STM") ~ "Shokoku",
      stringr::str_detect(sample, "SM|YKD") ~ "Shamo",
      .default = "Commercial"
    ),
    ancestor = dplyr::case_when(
      sample %in% c("ONG", "TOT", "SHK", "MHK", "KKW", "TMR", "KOE", "STM") ~ "SHK",
      stringr::str_detect(sample, "SM|YKD") ~ "OSM",
      sample %in% c("TKK", "NGY", "KMT") ~ "NGY",
      .default = sample
    ),
    ya = dplyr::case_when(
      sample == "SHK" ~ 1000,
      sample == "OSM" ~ 400,
      sample == "NGY" ~ 150,
      ancestor == "NGY" ~ 100,
      group == "Jidori" ~ 2000,
      sample %in% c("UKK", "CHB") ~ 250,
      group == "Commercial" ~ 150,
      .default = 100
  )) |>
  dplyr::group_by(group) |>
  dplyr::mutate(sample = forcats::fct_infreq(sample))
  
 
ggplot(df_sample) +
  aes(x = ya, y = sample, group = group) +
  geom_rect(aes(xmin = 1700, xmax = 2300, ymin = -Inf, ymax = Inf), fill = "#DDDDDD") +
  geom_rect(aes(xmin = 832, xmax = 1230, ymin = -Inf, ymax = Inf), fill = "#DDDDDD") +
  geom_rect(aes(xmin = 420, xmax = 186, ymin = -Inf, ymax = Inf), fill = "#DDDDDD") +
  geom_point(aes(color = group)) +
  geom_segment(aes(x = ya, xend = 0, y = ancestor, yend = sample, color = group)) +
  scale_x_reverse(
    breaks = c(2300, 1700, 1230, 832, 420, 186),
    labels = c("300 BC", " 300", "794", "1185", "1603", "1867"),
    expand = expansion(mult = c(0.05, 0))
    ) +
  scale_y_discrete(position = 'right') +
  scale_color_manual(values = c("Jidori" = "#F5C710", "Shamo" = "#D55E00", 
                                "Shokoku" = "#0072B2", "Commercial" = "#009E73")) +
  theme_test(base_size = 14) +
  theme(
    axis.title.y = element_blank(),
    legend.position = "none"
  )
