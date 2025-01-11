library(conflicted)
library(tidyverse)
library(phytools)
library(ggtree)
library(patchwork)


### Const and Files

BASESIZE = 16
LABELSIZE = 18

acc2chr= readr::read_tsv("data/sequence_report.tsv") |>
  dplyr::rename(CHROM = `RefSeq seq accession`, chr = `Chromosome name`) |>
  dplyr::select(CHROM, chr)

accession2name = readr::read_tsv("data/sra_accession.tsv") |>
  dplyr::select(Run, sample_id, Instrument, system) |>
  dplyr::rename(label = Run)

jpool_names = readxl::read_excel(path = "data/Japanese_Chicken_DNAseq.xlsx", sheet = "Sheet1") |>
  dplyr::select(sample, sample_name)

### PCA

eigenvec2 = readr::read_delim("out/batch/popstr.pca.after.eigenvec") |>
  dplyr::rename(label = IID) |>
  dplyr::full_join(accession2name, by = "label") |>
  dplyr::mutate(sample_id = dplyr::if_else(is.na(sample_id), label, sample_id)) |>
  dplyr::mutate(group = dplyr::case_when(
    stringr::str_detect(sample_id, "KNC") ~ "South Korea",
    stringr::str_detect(sample_id, "RJF") ~ "RJF",
    stringr::str_detect(sample_id, "Thailand|Vietnam|SGN|OMK|NKH") ~ "South-East Asia",
    sample_id %in% jpool_names$sample ~ "Japan",
    .default = "China"
  )) |>
  tidyr::replace_na(list(Instrument = "Illumina NovaSeq 6000", system = "NovaSeq")) |>
  dplyr::mutate(Group = group,Sequencer = system) |>
  tidyr::drop_na()

eigenval2 = readr::read_csv("out/batch/popstr.pca.after.eigenval", col_names = "V1")
df2 = data.frame(pc = 1:nrow(eigenval2), eigenval2/sum(eigenval2)*100)


fa = ggplot(eigenvec2) +
  aes(PC1, PC2, color = Group) +
  geom_point(size = 2, alpha = 0.9) +
  scale_color_discrete(type = palette.colors(palette = "Okabe-Ito")[-1]) +
  coord_fixed() +
  xlab(paste0("PC1 (", round(df2[1,2], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(df2[2,2], digits = 2), "%)")) +
  theme_test(base_size = BASESIZE) +
  theme(legend.position = "bottom", legend.title = element_blank()) + 
  guides(color = guide_legend(nrow = 3))  
fa


### NJtree

rapidnj = ape::read.tree("out/phylo/popstr.sub.nj.tree") |>
  dplyr::as_tibble() |>
  dplyr::mutate(label = stringr::str_remove_all(label, "'")) |>
  dplyr::left_join(accession2name, by = "label") |>
  dplyr::left_join(jpool_names, dplyr::join_by(label == sample)) |>
  dplyr::mutate(label = dplyr::if_else(is.na(sample_id), sample_name, sample_id)) |>
  dplyr::select(parent, node, branch.length, label) |>
  ape::as.phylo()


### Admixture

fam = readr::read_tsv("out/admixture/popstr.snp.filt.fam", col_names = FALSE) |>
  dplyr::rename(Run = X1) |>
  dplyr::left_join(accession2name, by = dplyr::join_by(Run == label)) |>
  dplyr::left_join(jpool_names, dplyr::join_by(Run == sample)) |>
  dplyr::mutate(sample_id = dplyr::if_else(is.na(sample_id), sample_name, sample_id)) |>
  dplyr::select(Run, sample_id)

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

df = dplyr::bind_rows(q2, q3, q4) |>
  dplyr::mutate(group = dplyr::case_when(
    stringr::str_detect(sample_id, "KNC") ~ "South Korea",
    stringr::str_detect(sample_id, "RJF") ~ "RJF",
    stringr::str_detect(sample_id, "Thailand|Vietnam|SGN|OMK|NKH") ~ "South-East Asia",
    sample_id %in% jpool_names$sample_name ~ "Japan",
    .default = "China"
  ))

tip2grp = df |>
  dplyr::select(sample_id, group) |>
  dplyr::distinct(sample_id, .keep_all = TRUE) |>
  dplyr::rename(label = sample_id)

pnj2 = rapidnj |>
  dplyr::full_join(tip2grp, by = "label") |>
  ggtree(aes(color = group)) +
  geom_tiplab(size = 2, align = TRUE) +
  scale_x_continuous(expand = expansion(mult = c(.0, .5))) +
  scale_color_discrete(type = palette.colors(palette = "Okabe-Ito")[-1]) +
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

padm = adm |> aplot::insert_left(pnj2, width = .5) |> aplot::as.patchwork()
padm

cve = readr::read_delim("out/admixture/CV-error.txt",
                        col_names = c("c", "e", "K", "error")) |>
  dplyr::mutate(K = as.integer(stringr::str_extract(K, "[0-9]+"))) |>
  #dplyr::filter(K > 1) |>
  ggplot() +
  aes(x = as.integer(K), y = error) +
  geom_point(size = 2) +
  geom_line(linewidth = 1) +
  # scale_x_continuous(breaks = seq(1, 8), labels = seq(1, 8)) +
  labs(x = "K", y = "CV Error") +
  theme_classic(base_size = BASESIZE)
cve


### Align plots

p1 = cowplot::plot_grid(fa, cve, nrow = 2, labels = c("a", "b"), label_size = LABELSIZE, rel_heights = c(2, 1))
p2 = cowplot::plot_grid(p1, NULL, padm, ncol = 3, labels = c("", "", "c"), label_size = LABELSIZE, rel_widths = c(1, .1, 1))
ggsave("images/structure_analysis.tif", p2, w = 8, h = 8, bg = "#FFFFFF")

