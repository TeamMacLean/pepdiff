# Read proteomics data into a pepdiff_data object

Imports a CSV file containing PRM proteomics data and creates a
pepdiff_data object suitable for analysis with \[compare()\].

## Usage

``` r
read_pepdiff(file, id, gene, value, factors, replicate, tech_rep = NULL)
```

## Arguments

- file:

  Path to CSV file

- id:

  Column name containing peptide identifiers

- gene:

  Column name containing gene identifiers

- value:

  Column name containing abundance values

- factors:

  Character vector of column names to use as experimental factors

- replicate:

  Column name containing biological replicate identifiers

- tech_rep:

  Optional column name containing technical replicate identifiers. If
  provided, data will NOT be automatically combined - use
  \[combine_tech_reps()\] explicitly after import.

## Value

A pepdiff_data object with components:

- data:

  Tibble with columns: peptide, gene_id, \[factors\], bio_rep, value

- factors:

  Character vector of factor names

- design:

  Tibble of factor combinations with n_reps, n_peptides

- missingness:

  Tibble of peptide missingness statistics

- peptides:

  Character vector of unique peptide IDs

- call:

  The original function call

## Examples

``` r
if (FALSE) { # \dontrun{
# Simple import with one factor
dat <- read_pepdiff(
  "data.csv",
  id = "peptide_sequence",
  gene = "gene_name",
  value = "intensity",
  factors = "treatment",
  replicate = "bio_rep"
)

# Multi-factor import
dat <- read_pepdiff(
  "data.csv",
  id = "peptide",
  gene = "gene_id",
  value = "total_area",
  factors = c("treatment", "timepoint"),
  replicate = "bio_rep"
)
} # }
```
