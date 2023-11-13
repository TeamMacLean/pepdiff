anon <- import_data("~/Desktop/pepdiff/inst/extdata/anon.csv",                  gene_id = "gene_name",
treatment = "treatment_name")

comparisons <- data.frame(
control = c('665e6428','665e6428'),
c_seconds = c(0,150),
treatment = c('665e6428','665e6428'),
t_seconds = c(150, 0)
)

many <- compare_many(anon, comparisons, tests = c("bootstrap_t", "rank_product"))
