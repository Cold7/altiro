# Human ovarian cancer: OVCAR5 cisplatin-sensitive cells (group A) vs their cisplatin-resistant derivative
# (group B), treated with vehicle or the DHODH inhibitor brequinar (blocks de novo pyrimidine synthesis).
# Data: GEO GSE283381 (raw HTSeq counts deposited by the authors, 3 replicates), bundled with the package.
# Run from this folder:  Rscript run_r.R

library(orthoresponse)

ex <- load_example("ovarian_brequinar")   # counts: 59,453 genes x 12 samples
ex$samples[, c("line", "treatment")]

res <- orthoresponse(ex$counts, ex$samples,
                     group_col = "line", condition_col = "treatment",
                     group_a = "CisS", group_b = "CisR",        # A = sensitive, B = resistant
                     control = "Vehicle", stress = "Brequinar")
res
# tables + figures (model_patterns.png/pdf, model_summary.png/pdf)
save_orthoresponse(res, "results_r", condition_labels = c("Veh", "BQ"))

strict <- classify_responses(res$contrasts, preset = "strict", group_names = c("CisS", "CisR"))
cat("\nStress-ready genes, paper vs strict preset:",
    sum(res$classification$model == "stress_ready"), "->", sum(strict$model == "stress_ready"), "\n")
print(table(strict$label))
print(primed_asymmetry(strict))

sr3 <- subset(strict, label == "stress_ready_3")
sr3 <- sr3[order(-sr3$a_response_log2FoldChange), c("a_response_log2FoldChange", "basal_log2FoldChange")]
print(round(head(sr3, 20), 2))
