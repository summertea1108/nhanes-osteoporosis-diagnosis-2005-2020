# Run ONCE on an OSC login node (it needs internet):  Rscript install_packages.R
lib <- Sys.getenv("R_LIBS_USER")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
pkgs <- c("haven", "survey", "dplyr", "ggplot2")
need <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(need)) install.packages(need, lib = lib, repos = "https://cloud.r-project.org")
for (p in pkgs) cat(sprintf("%-8s %s\n", p, if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else "NOT INSTALLED"))
