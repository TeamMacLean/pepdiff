# read data from a file

reads data, renames columns appropriately, discards unused columns,
factors and reorders, discards duplicate rows

## Usage

``` r
import_data(
  file,
  treatment = "genotype",
  bio_rep = "bio_rep",
  tech_rep = "tech_rep",
  quant = "total_area",
  seconds = "seconds",
  gene_id = "gene_id",
  peptide = "peptide_sequence"
)
```

## Arguments

- file:

  Path to the file to load - must be a csv file

- treatment:

  Column containing the treatment of the observation

- bio_rep:

  Column containing the biological replicate of the observation

- tech_rep:

  Column containing the technical replicate of the observation

- quant:

  Column containing the quantitation data

- seconds:

  Column containing timepoint of observation

- gene_id:

  Column containing the id of the gene this hit

- peptide:

  Column containing the sequence of this peptide

## Value

tibble with columns id, gene_id, peptide, treatment, seconds, bio_rep,
tech_rep, quant
