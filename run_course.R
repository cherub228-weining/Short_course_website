# Run from the repository root: Rscript code/R/run_course.R
if (!file.exists("code/R/course_demo.R")) {
  stop("Run this script from the Short_course repository root.")
}
source("code/R/course_demo.R")
dir.create("data/generated", recursive = TRUE, showWarnings = FALSE)
write.csv(d, "data/generated/synthetic_panel.csv", row.names = FALSE)
cat("Synthetic panel saved to data/generated/synthetic_panel.csv\n")

