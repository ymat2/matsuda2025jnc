library(conflicted)
library(tidyverse)
library(readxl)


## Maps 

world_map = ggplot2::map_data("world") |> dplyr::as_tibble()

map_asia = ggplot(world_map) + 
  aes(long, lat) +
  geom_map(aes(map_id = region), map = world_map, fill = "#a6bddb", linewidth = 0) +
  ggplot2::annotate("point", x = 100, y = 20, size = 40, color = "#F0E442", alpha = .5, shape = 16) +
  ggplot2::annotate("text", label = "Early\ndomestication", x = 100, y = 20, size = 5, color = "#444444") +
  ggplot2::annotate("point", x = 140, y = 10, size = 40, color = "#999999", alpha = .5, shape = 16) +
  ggplot2::annotate("text", label = "Other regions\n(USA, Europe, etc.)", x = 140, y = 10, size = 5, color = "#444444") +
  scale_x_continuous(limits = c(90, 150), labels = scales::label_number(suffix = "°N")) +
  scale_y_continuous(limits = c(0, 50), labels = scales::label_number(suffix = "°E")) + 
  coord_fixed(ratio = 1) +
  theme_test(base_size = 14) +
  theme(axis.title = element_blank())
map_asia

ggsave("images/map_asia2.png", map_asia, w = 6, h = 5)


#####

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
