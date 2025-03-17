library(conflicted)
library(tidyverse)
library(readxl)

source("./rstats/common_settings.R")


#### Table 1 : Japanese breeds used in this study ####

sample_info_1 = readxl::read_excel(path = "data/Japanese_Chicken_DNAseq.xlsx", sheet = "Sheet1") |>
  #dplyr::filter(poolseq > 1) |>
  dplyr::filter(country == "Japan") |>
  dplyr::mutate(
    breed = dplyr::if_else(
      stringr::str_detect(common_name, "shamo|ABRC"),
      stringr::str_split(common_name, "-", simplify = TRUE)[,1],
      common_name
    ),
    Source = dplyr::if_else(poolseq > 1, "This study", "Bendesky et al. 2024")
  ) |>
  dplyr::select(sample_name, common_name, breed, poolseq, Source) |>
  dplyr::arrange(desc(Source), sample_name) |>
  dplyr::rename(
    "Sample name" = sample_name,
    "Population" = common_name,
    "Breed" = breed,
    "Number of individuals" = poolseq,
  )

knitr::kable(sample_info_1)


#### Table 2 : Other regions breeds sequenced in this study ####

sample_info_2 = readxl::read_excel(path = "data/Japanese_Chicken_DNAseq.xlsx", sheet = "Sheet1") |>
  dplyr::filter(poolseq > 1) |>
  dplyr::filter(country != "Japan") |>
  dplyr::mutate(Source = "This study") |>
  dplyr::select(sample_name, common_name, poolseq, country, Source) |>
  dplyr::arrange(desc(Source), sample_name) |>
  dplyr::rename(
    "Sample name" = sample_name,
    "Population" = common_name,
    "Origin" = country,
    "Number of individuals" = poolseq,
  )

knitr::kable(sample_info_2)


#### Table S1 ####

sample_info_S = readxl::read_excel(path = "data/Japanese_Chicken_DNAseq.xlsx", sheet = "Sheet1") |>
  dplyr::filter(poolseq > 1) |>
  dplyr::mutate(
    provider = stringr::str_replace_all(provider, "Breeder", "1"),
    provider = stringr::str_replace_all(provider, "Hiroshima University", "2"),
    provider = stringr::str_replace_all(provider, "Livestock Experiment Station", "3"),
    provider = stringr::str_replace_all(provider, "National Livestock Breeding Center", "4"),
    provider = stringr::str_replace_all(provider, "Avian Bioscience Research Center", "5"),
    breed = dplyr::if_else(
      stringr::str_detect(common_name, "shamo"),
      stringr::str_split(common_name, "-", simplify = TRUE)[,1],
      common_name
    )
    ) |>
  dplyr::mutate(Description = dplyr::case_when(
    stringr::str_detect(sample_name, "SM") ~ "Gamecocks",
    country != "Japan" ~ "Other regions",
    .default = "Native breeds"
  )) |>
  dplyr::mutate(Sequencer = dplyr::if_else(provider == "5", "Illumina HiSeq Ten X", "Illumina NovaSeq 6000")) |>
  dplyr::select(sample_name, common_name, breed, sampling_location, provider, Sequencer, Description) |>
  dplyr::arrange(Description, sample_name) |>
  dplyr::rename(
    "Sample name" = sample_name,
    "Population" = common_name,
    "Breed" = breed,
    "Sampling location" = sampling_location,
    "Provider" = provider
  )

knitr::kable(sample_info_S)


#### Table S2 ####

sample_info = readxl::read_excel(path = "data/Japanese_Chicken_DNAseq.xlsx", sheet = "Sheet1") |>
  tidyr::drop_na() |>
  dplyr::select(sample, sample_name, poolseq)

jpool_flagstat = readr::read_tsv("../jpool/out/flag_summary.tsv")
jpool_mapping_stat = readr::read_tsv("../jpool/out/mapping_summary.tsv")
jpool_stat = dplyr::full_join(jpool_flagstat, jpool_mapping_stat, by = "sample")

abrc_flagstat = readr::read_tsv("../ABRC/out/abrc_flag_summary.tsv")
abrc_mapping_stat = readr::read_tsv("../ABRC/out/abrc_mapping_summary.tsv")
abrc_stat = dplyr::full_join(abrc_flagstat, abrc_mapping_stat, by = "sample")

sample_stats = dplyr::bind_rows(jpool_stat, abrc_stat) |>
  dplyr::inner_join(sample_info, by = "sample") |>
  dplyr::filter(poolseq > 1) |>
  dplyr::arrange(sample_name) |>
  dplyr::select(sample_name, poolseq, num_reads, prop_mapped, meandp, meancov) |>
  dplyr::mutate(num_reads = format(num_reads, big.mark = ",", big.interval = 3L)) |>
  dplyr::rename(
    "Sample" = sample_name,
    "Number of individuals" = poolseq,
    "Number of reads" = num_reads,
    "Properly mapped reads (%)" = prop_mapped,
    "Mean depth" = meandp,
    "Mean coverage (%)" = meancov
  )

knitr::kable(sample_stats, digits = 2L)


#### Table S3 ####

ext_mapping_summary = readr::read_tsv("out/mapping_summary.tsv") |>
  dplyr::rename(Run = sample) |>
  dplyr::select(Run, meandp, meancov)

jpool_mapping_summary = readr::read_tsv("../jpool/out/mapping_summary.tsv") |>
  dplyr::rename(sample_id = sample) |>
  dplyr::select(sample_id, meandp, meancov)

extern_sample = readr::read_tsv("data/sra_accession.tsv") |>
  dplyr::select(!system) |>
  dplyr::left_join(ext_mapping_summary, by = "Run") |>
  dplyr::rows_patch(jpool_mapping_summary, by = c("sample_id"), unmatched = "ignore") |>
  dplyr::left_join(jpool_names, by = join_by(sample_id == sample)) |>
  dplyr::mutate(
    sample_name = dplyr::if_else(is.na(sample_name), sample_id, sample_name),
    Instrument = stringr::str_remove(Instrument, "^Illumina ")
    ) |>
  dplyr::select(Run, sample_name, source, country, Instrument, meandp, meancov) |>
  dplyr::rename(
    "Sample name" = sample_name, 
    "Source" = source, 
    "Country" = country,
    "Mean depth" = meandp,
    "Mean coverage (%)" = meancov
  )

knitr::kable(extern_sample, digits = 2L)
