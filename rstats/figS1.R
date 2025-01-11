library(conflicted)
library(tidyverse)


acc2chr= readr::read_tsv("data/sequence_report.tsv") |>
  dplyr::rename(CHROM = `RefSeq seq accession`, chr = `Chromosome name`) |>
  dplyr::select(CHROM, chr)

accession2name = readr::read_tsv("data/sra_accession.tsv") |>
  dplyr::select(Run, sample_id, Instrument, system)


##### PCA before removing batch effect #####

eigenvec = readr::read_delim("out/batch/popstr.pca.before.eigenvec") |>
  dplyr::rename(Run = IID) |>
  dplyr::full_join(accession2name, by = "Run") |>
  dplyr::mutate(sample_id = dplyr::if_else(is.na(sample_id), Run, sample_id)) |>
  dplyr::mutate(group = dplyr::case_when(
    stringr::str_detect(sample_id, "KNC") ~ "Korea",
    stringr::str_detect(sample_id, "RJF") ~ "RJF",
    stringr::str_detect(sample_id, "gamecock_") ~ "Other gamecocks",
    stringr::str_detect(sample_id, "SM|YKD") ~ "Japanese Shamo",
    stringr::str_detect(sample_id, "Thailand|Vietnam|SGN|OMK|NKH") ~ "South-East Asia",
    stringr::str_detect(sample_id, "[1-6]$") ~ "China",
    .default = "Japan"
  )) |>
  tidyr::replace_na(list(Instrument = "Illumina NovaSeq 6000", system = "NovaSeq")) |>
  dplyr::mutate(Group = group, Sequencer = system) |>
  tidyr::drop_na()

eigenval = readr::read_csv("out/batch/popstr.pca.before.eigenval", col_names = "V1")
df = data.frame(pc = 1:nrow(eigenval), eigenval/sum(eigenval)*100)


fa = ggplot(eigenvec) +
  aes(PC1, PC2, color = Group) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_color_discrete(type = palette.colors(palette = "Okabe-Ito")[-1]) +
  coord_fixed() +
  xlab(paste0("PC1 (", round(df[1,2], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(df[2,2], digits = 2), "%)")) +
  labs(tag = "a") +
  theme_bw(base_size = 16) +
  theme(plot.tag.location = "plot")
fa

fb = ggplot(eigenvec) +
  aes(PC1, PC2, color = Sequencer) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_color_discrete(type = palette.colors(palette = "Set1")) +
  coord_fixed() +
  xlab(paste0("PC1 (", round(df[1,2], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(df[2,2], digits = 2), "%)")) +
  theme_bw(base_size = 16) +
  theme(plot.tag.location = "plot")
fb


##### Check distribution of GWAS p-value and Fst #####

gwas_result = readr::read_tsv("out/batch/gwas-result.txt") |>
  dplyr::mutate(extract = dplyr::if_else(FDR_BH > 0.05, "Y", "N"))

nall = nrow(gwas_result)
nexc = nrow(gwas_result |> dplyr::filter(extract == "N"))
pexc = round((nexc/nall)*100, digits = 2)

txt = "Number of removed SNPs:\n356,946/6,180,729 (5.78%)"

fc = ggplot(gwas_result) +
  geom_histogram(aes(x = FDR_BH, fill = extract), binwidth = 0.01, boundary = 0) +
  geom_vline(xintercept = 0.05, linetype = "longdash", color = "#BBBBBB") +
  annotate("text", label = txt, x = 0, y = 300000, hjust = 0, size = 5) +
  scale_x_continuous(breaks = c(0, 0.05, 1), labels = c(0, 0.05, 1)) +
  scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_long_scale())) +
  labs(x = "adjusted p-value of GWAS", y = "", tag = "b") +
  scale_fill_manual(values = c("#4d4d4d", "#BBBBBB")) +
  theme_test(base_size = 16) +
  theme(legend.position = "none", plot.tag.location = "plot")
fc


##### PCA after removing batch effect #####

eigenvec2 = readr::read_delim("out/batch/popstr.pca.after.eigenvec") |>
  dplyr::rename(Run = IID) |>
  dplyr::full_join(accession2name, by = "Run") |>
  dplyr::mutate(sample_id = dplyr::if_else(is.na(sample_id), Run, sample_id)) |>
  dplyr::mutate(group = dplyr::case_when(
    stringr::str_detect(sample_id, "KNC") ~ "Korea",
    stringr::str_detect(sample_id, "RJF") ~ "RJF",
    stringr::str_detect(sample_id, "gamecock_") ~ "Other gamecocks",
    stringr::str_detect(sample_id, "SM|YKD") ~ "Japanese Shamo",
    stringr::str_detect(sample_id, "Thailand|Vietnam|SGN|OMK|NKH") ~ "South-East Asia",
    stringr::str_detect(sample_id, "[1-6]$") ~ "China",
    .default = "Japan"
  )) |>
  tidyr::replace_na(list(Instrument = "Illumina NovaSeq 6000", system = "NovaSeq")) |>
  dplyr::mutate(Group = group, Sequencer = system) |>
  tidyr::drop_na()

eigenval2 = readr::read_csv("out/batch/popstr.pca.after.eigenval", col_names = "V1")
df2 = data.frame(pc = 1:nrow(eigenval2), eigenval2/sum(eigenval2)*100)


fd = ggplot(eigenvec2) +
  aes(PC1, PC2, color = Group) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_color_discrete(type = palette.colors(palette = "Okabe-Ito")[-1]) +
  coord_fixed() +
  xlab(paste0("PC1 (", round(df2[1,2], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(df2[2,2], digits = 2), "%)")) +
  labs(tags = "c") +
  theme_bw(base_size = 16) +
  theme(plot.tag.location = "plot")
fd

fe = ggplot(eigenvec2) +
  aes(PC1, PC2, color = Sequencer) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_color_discrete(type = palette.colors(palette = "Set1")) +
  coord_fixed() +
  xlab(paste0("PC1 (", round(df2[1,2], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(df2[2,2], digits = 2), "%)")) +
  theme_bw(base_size = 16)
fe


## cowplot

library(patchwork)

row1 = fa + fb + patchwork::plot_layout(guides = 'collect')
row2 = fd + fe + patchwork::plot_layout(guides = 'collect')
p = row1 / fc / row2 +
  patchwork::plot_layout(nrow = 3, heights = c(1, .7, 1)) &
  theme(plot.tag = element_text(face = "bold", size = 22))
ggsave("images/figS1.tif", p, w = 9, h = 10.5, bg = "#FFFFFF")

