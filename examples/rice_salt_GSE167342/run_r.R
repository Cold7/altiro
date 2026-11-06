# Rice salt stress: salt-sensitive PL12 (group A) vs salt-tolerant Saltol introgression line IL22 (group B).
# Data: GEO GSE167342 (leaves, control vs 24 h 80 mM NaCl, 3 replicates), bundled with the package.
# Run from this folder:  Rscript run_r.R

library(altiro)

ex <- load_example("rice_salt")   # counts: 36,844 genes x 12 samples (S40 ... S72)
ex$samples[, c("line", "treatment")]

res <- altiro(ex$counts, ex$samples,
              group_col = "line", condition_col = "treatment",
              group_a = "PL12", group_b = "IL22",   # A = sensitive, B = tolerant
              control = "Control", stress = "Salt")
res
# tables + figures (model_patterns.png/pdf, model_summary.png/pdf)
save_altiro(res, "results_r", condition_labels = c("Ctrl", "NaCl"))

strict <- classify_responses(res$contrasts, preset = "strict", group_names = c("PL12", "IL22"))
cat("\nPrimed genes, lenient vs strict preset:",
    sum(res$classification$model == "primed"), "->", sum(strict$model == "primed"), "\n")
