




test_get_eb <- function(treatment, control) {
  #check this is in correct order wrt peptides... should be.
  m <- log(cbind(treatment, control))
  tc <- c(rep("treatment", 3), rep("control", 3))
  design <- model.matrix(~factor(tc))
  colnames(design) <- c("trt","ctrl")

  fit <- limma::lmFit(m, design)
  fit <- limma::eBayes(fit)
  betahat <- fit$coefficients[,2]
  sehat <- fit$stdev.unscaled[,2]*fit$sigma
  fit.vash <- vashr::vash(sehat=sehat, df=fit$df.residual[1], betahat=betahat)


  data.frame(
    eb_vs_p_val = fit.vash$pvalue,
    eb_vs_fdr = fit.vash$qvalue,
    eb_p_val = fit$p.value[,2],
    eb_fdr = p.adjust(fit$p.value[,2], method="fdr")
  )

}

