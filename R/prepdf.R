#convert seconds key into an integer
get_seconds <- function(x){
  idx <- x %in% c("01", "05")
  x[idx] <- 0
  idx <- x %in% c("02", "06")
  x[idx] <- 150
  idx <- x %in% c("03", "03r", "07")
  x[idx] <- 300
  idx <- x %in% c("04", "04r", "08")
  x[idx] <- 600
  return(x)
}

#convert tech rep key into a number
get_tech_rep <- function(x){
  idx <- x == "01"
  x[idx] <- 1
  idx <- x == "01a"
  x[idx] <- 2
  idx <- x == "02"
  x[idx] <- 3
  idx <- x == "02a"
  x[idx] <- 4
  return(x)
}
library(readxl)
library(here)
library(magrittr)
library(dplyr)

df <- read_excel(here("test_data","Data_for_Dan_200303.xlsx"),sheet = "Microsomes") %>%
  dplyr::rename_all(list(~ tolower(stringr::str_replace_all(., "\\s", "_")))) %>%
  tidyr::separate(replicate_name, into = c("bio_rep_key", "tech_rep_key", "bio_sample_key"), sep = "_") %>%
  dplyr::mutate(from_redo = dplyr::if_else(stringr::str_detect(bio_sample_key, "r$"), TRUE, FALSE) ) %>%
  dplyr::mutate(bio_sample_key = stringr::str_replace_all(bio_sample_key, "r", "")) %>%
  dplyr::mutate(genotype = dplyr::if_else(bio_sample_key %in% c("01", "02", "03", "03r", "04", "04r"), "Col-0", "bik1pbl")) %>%
  dplyr::mutate(seconds = get_seconds(bio_sample_key)) %>%
  dplyr::mutate(bio_rep = dplyr::if_else(bio_rep_key == "pde180416", 1, dplyr::if_else(bio_rep_key == "pde180620", 2,3))) %>%
  dplyr::mutate(tech_rep = get_tech_rep(tech_rep_key) ) %>%
  dplyr::mutate(pep_id = paste0(molecule_list_name, peptide_modified_sequence)) %>%
  dplyr::mutate(tag = paste(molecule_list_name, peptide_modified_sequence,bio_rep,tech_rep,genotype,seconds)) %>%
  dplyr::mutate(total_area = dplyr::if_else(total_area == 0, NA_real_, total_area)) %>%
  dplyr::filter(! molecule_list_name == 'AT3G13530.1', ! molecule_list_note == 'TOM') %>%
  dplyr::distinct()

  readr::write_csv(df, here("test_data", "large.csv"))
