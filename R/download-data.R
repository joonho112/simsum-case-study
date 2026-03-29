# download-data.R
# Helper script to verify/obtain required data files for the replication course
#
# The simulation results and TSL15 dataset are from the OSF repository
# for Pustejovsky & Tipton (2022): https://osf.io/x8yre/
#
# Required files in data/:
#   - Tanner-Smith-Lipsey-2015-subset.rds  (42 KB)
#   - dv-model-data.Rdata                   (8 KB)
#   - Ctype-model-data.Rdata                (7 KB)
#   - RVE-simulation-dvcat-results.Rdata    (3.5 MB)
#   - RVE-simulation-Ctype-results.Rdata    (2.8 MB)

required_files <- c(
  "data/Tanner-Smith-Lipsey-2015-subset.rds",
  "data/dv-model-data.Rdata",
  "data/Ctype-model-data.Rdata",
  "data/RVE-simulation-dvcat-results.Rdata",
  "data/RVE-simulation-Ctype-results.Rdata"
)

cat("Checking required data files...\n\n")
for (f in required_files) {
  if (file.exists(f)) {
    cat(sprintf("  [OK] %s (%.1f KB)\n", f, file.size(f) / 1024))
  } else {
    cat(sprintf("  [MISSING] %s\n", f))
  }
}

missing <- required_files[!file.exists(required_files)]
if (length(missing) > 0) {
  cat("\n", length(missing), "file(s) missing. Download from:\n")
  cat("  https://osf.io/x8yre/\n\n")
  cat("Place the files in the data/ directory of this project.\n")
} else {
  cat("\nAll data files present. Ready to render.\n")
}
