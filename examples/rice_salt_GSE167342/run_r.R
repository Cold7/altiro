# Rice salt stress: salt-sensitive PL12 (group A) vs salt-tolerant Saltol introgression line IL22 (group B).
# Data: GEO GSE167342 (leaves, control vs 24 h 80 mM NaCl, 3 replicates), bundled with the package.
# Run from this folder:  Rscript run_r.R

library(orthoresponse)

ex <- load_example("rice_salt")   # counts: 36,844 genes x 12 samples (S40 ... S72)
ex$samples[, c("line", "treatment")]

res <- orthoresponse(ex$counts, ex$samples,
                     group_col = "line", condition_col = "treatment",
                     group_a = "PL12", group_b = "IL22",   # A = sensitive, B = tolerant
                     control = "Control", stress = "Salt")
res
save_orthoresponse(res, "results_r")

strict <- classify_responses(res$contrasts, preset = "strict", group_names = c("PL12", "IL22"))
cat("\nStress-ready genes, paper vs strict preset:",
    sum(res$classification$model == "stress_ready"), "->", sum(strict$model == "stress_ready"), "\n")

png(file.path("results_r", "model_patterns.png"), width = 1600, height = 1500, res = 170)
plot_model_patterns(res$classification, condition_labels = c("Ctrl", "NaCl"))
dev.off()
