# common settings


### Const and Files

BASESIZE = 16
LABELSIZE = 18

### sample names

accession2name = readr::read_tsv("data/sra_accession.tsv") |>
  dplyr::filter(country != "Japan")

abrc = c("BLE", "BMC", "CAL", "DD", "GB", "GSP", "JB", "RIR", "SIL", "WLG")

jpool_names = readxl::read_excel(path = "data/Japanese_Chicken_DNAseq.xlsx", sheet = "Sheet1") |>
  # dplyr::filter(poolseq > 1) |>
  dplyr::select(sample, sample_name, country, poolseq) |>
  dplyr::rename(Run = sample, sample_id = sample_name) |>
  dplyr::mutate(source = dplyr::if_else(poolseq > 1, "This study", "Bendesky et al. 2024")) |>
  dplyr::mutate(
    Instrument = dplyr::if_else(Run %in% abrc, "Illumina HiSeq Ten X", "Illumina NovaSeq 6000"), 
    system = dplyr::if_else(Run %in% abrc, "HiSeq", "NovaSeq")
  ) |>
  dplyr::select(!poolseq)

sample2group = dplyr::bind_rows(accession2name, jpool_names) |>
  dplyr::mutate(
    group = dplyr::case_when(
      stringr::str_detect(sample_id, "RJF") ~ "RJF",
      stringr::str_detect(country, "Thailand|Vietnam") ~ "Southeast Asia",
      stringr::str_detect(country, "UK|USA|Egypt|Italy|Chile") ~ "Other regions",
      .default = country
    )
  )
rm(accession2name, jpool_names)

color_scale = c(
  "RJF" = "#009E73",
  "Southeast Asia" = "#F0E442",
  "China" = "#E69F00",
  "South Korea" = "#CC79A7",
  "Japan" = "#56B4E9",
  "Other regions" = "#999999"
)


