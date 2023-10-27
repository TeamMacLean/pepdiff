generate_sample_data <- function(n=50000, npeps = 30, treatment = 2, tech_reps = 5, bio_reps = 3, seconds = 4, missing_obs_proportion = 0.05, na_proportion= 0.1) {
  set.seed(123)  # for reproducibility

  peptides <- stringi::stri_rand_strings(npeps, 9, pattern = "[A-Z]")
  # Sample data for testing
  d <- tibble::tibble(
    treatment = sample(paste0("treatment", letters[1:treatment]), n, replace = TRUE),
    bio_rep = factor(sample(1:bio_reps, n, replace = TRUE)),
    tech_rep = sample(letters[1:tech_reps], n, replace = TRUE),
    quant = runif(n, min = 0, max = 300),
    seconds = factor(sample((1:seconds) * 60, n, replace = TRUE)),
    gene_id = sample(paste0("G", 1:npeps), n, replace = TRUE),
    peptide = sample(peptides, n, replace=TRUE)
  )

  # Determine the number of rows to remove to achieve the specified missing_proportion
  num_missing_rows <- round(n * missing_obs_proportion)

  # Randomly select and remove rows to create missing quant values
  if (num_missing_rows > 0) {
    remove_indices <- sample(1:n, num_missing_rows)
    d <- d[-remove_indices, ]
  }

  if (na_proportion > 0){
    nas <- sample(c(1, NA), nrow(d),prob = c(1-na_proportion, na_proportion), replace=TRUE)
    d$quant <- d$quant * nas
  }

  return(d)
}

