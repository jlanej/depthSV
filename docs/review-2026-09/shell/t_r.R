suppressPackageStartupMessages(library(data.table))
S <- "/private/tmp/claude-501/-Users-Kitty-git-depthSV/1db72424-c4a1-45f4-b81e-5b466f3ec8a4/scratchpad/review_shell"
REPO <- "/Users/Kitty/git/depthSV/.claude/worktrees/agent-a08479f00845fba52"

cat("=== E3: choose_ndim.R fit on a steep spectrum (all 200 reported values are signal) ===\n")
# Extract fit_one/mp_quantile_fn from the script by sourcing only the function definitions.
src <- readLines(file.path(REPO, "example/1000G_highcov/R/choose_ndim.R"))
start <- grep("^mp_quantile_fn <- function", src); end <- grep("^runs <- strsplit", src) - 1
eval(parse(text = src[start:end]))
n <- 3202; p <- 142070
set.seed(1)
# spectrum decaying like a power law, never reaching the MP bulk within 200 ranks
sv_steep <- 3000 * (1:200)^(-0.9)
r <- tryCatch(fit_one(sv_steep, n, p, 0.01, 20L), error = function(e) paste("ERROR:", conditionMessage(e)))
if (is.character(r)) cat(r, "\n") else cat("k =", r$k, "status =", r$status, "\n")
# Also a spectrum where exactly ~185 values are above the edge
sigma <- 1; edge <- sigma * (sqrt(p) + sqrt(n))
sv_185 <- c(edge * 1.5 * seq(3, 1.05, length.out = 185), edge * 0.98 * seq(1, 0.9, length.out = 15))
r2 <- tryCatch(fit_one(sv_185, n, p, 0.01, 20L), error = function(e) paste("ERROR:", conditionMessage(e)))
if (is.character(r2)) cat(r2, "\n") else cat("k =", r2$k, "status =", r2$status, "\n")
# And the normal case for reference: 40 signal values then a bulk
qf <- mp_quantile_fn(n / p); s_unit <- sqrt(p * qf(1 - (seq_len(200) - 0.5) / n))
sv_ok <- c(edge * seq(4, 1.2, length.out = 40), s_unit[41:200])
r3 <- tryCatch(fit_one(sv_ok, n, p, 0.01, 20L), error = function(e) paste("ERROR:", conditionMessage(e)))
if (is.character(r3)) cat(r3, "\n") else cat("reference: k =", r3$k, "status =", r3$status, "\n")

cat("\n=== E4: example read_shards on a header-only shard and mixed shards ===\n")
d <- file.path(S, "shards"); dir.create(d, showWarnings = FALSE)
hdr <- "#CHROM\tSTART\tEND\tRegion\tN\tNCase\tNControl\tEstimate\tStd. Error\tt value\tPr(>|t|)"
writeLines(hdr, file.path(d, "a.txt")); system(paste("gzip -f", shQuote(file.path(d, "a.txt")))); file.rename(file.path(d, "a.txt.gz"), file.path(d, "mtdna_cn.linear.chr20_1-1000000.txt.gz"))
writeLines(c(hdr, "chr20\t0\t1000\tchr20:0-1000\t64\t64\t64\t0.1\t0.05\t2\t0.04"), file.path(d, "b.txt")); system(paste("gzip -f", shQuote(file.path(d, "b.txt")))); file.rename(file.path(d, "b.txt.gz"), file.path(d, "mtdna_cn.linear.chr20_1000001-2000000.txt.gz"))
# stray files that a glob might catch
file.create(file.path(d, "mtdna_cn.linear.chr20_1-1000000.txt.gz.tmp.gz"))
file.create(file.path(d, "mtdna_cn.linear.chr20_1-1000000.txt.gz.tbi"))
file.create(file.path(d, "mtdna_cn.linear.chr20_1-1000000.log"))
file.create(file.path(d, "mtdna_cn_adj.linear.chr20_1-1000000.txt.gz"))
file.create(file.path(d, "mtdna_cn_null_adj.linear.chr20_1-1000000.txt.gz"))
file.create(file.path(d, "log2_mtdna_cn.linear.chr20_1-1000000.txt.gz"))
cat("glob for mtdna_cn.linear:", basename(Sys.glob(file.path(d, "mtdna_cn.linear.*.txt.gz"))), "\n")
src2 <- readLines(file.path(REPO, "example/1000G_highcov/R/evaluate.R"))
s2 <- grep("^read_shards <- function", src2); e2 <- grep("^chrom_class <- function", src2) - 1
eval(parse(text = src2[s2:e2]))
res <- tryCatch(read_shards(d, "mtdna_cn", "linear"), error = function(e) paste("ERROR:", conditionMessage(e)), warning = function(w) paste("WARNING:", conditionMessage(w)))
if (is.character(res)) cat(res, "\n") else { cat("files:", length(res$files), "rows:", nrow(res$dt), "\n"); print(res$dt) }
# header-only alone
d2 <- file.path(S, "shards2"); dir.create(d2, showWarnings = FALSE)
file.copy(file.path(d, "mtdna_cn.linear.chr20_1-1000000.txt.gz"), file.path(d2, "mtdna_cn.linear.chr20_1-1000000.txt.gz"), overwrite = TRUE)
res2 <- tryCatch(read_shards(d2, "mtdna_cn", "linear"), error = function(e) paste("ERROR:", conditionMessage(e)), warning = function(w) paste("WARNING:", conditionMessage(w)))
if (is.character(res2)) cat("header-only alone:", res2, "\n") else cat("header-only alone: files", length(res2$files), "rows", nrow(res2$dt), "cols", ncol(res2$dt), "\n")

cat("\n=== E5: fread of the analyses.tsv that 01_prepare_inputs.sh writes ===\n")
an <- c("# depthSV 1000G example analyses. Format: name<TAB>method<TAB>model",
        "# The depth term must be named cov_resids (see conf/phenotypes.example.tsv).", "#",
        "# SEX appears twice on purpose. The linear run carries the truth check:",
        "# sex separates depth on chrX/Y so completely that the logistic Wald z",
        "# collapses there (Hauck-Donner), so the logistic run exercises that",
        "# engine and documents the collapse rather than asserting rank.",
        "mtdna_cn\tlinear\tMTDNA_CN~cov_resids", "log2_mtdna_cn\tlinear\tLOG2_MTDNA_CN~cov_resids",
        "mtdna_cn_null\tlinear\tMTDNA_CN_NULL~cov_resids", "sex_linear\tlinear\tSEX~cov_resids",
        "inferred_sex\tlogistic\tSEX~cov_resids")
writeLines(an, file.path(S, "analyses.tsv"))
m <- withCallingHandlers(fread(file.path(S, "analyses.tsv"), header = FALSE, sep = "\t", col.names = c("name", "method", "model")),
                         warning = function(w) { cat("WARNING:", conditionMessage(w), "\n"); invokeRestart("muffleWarning") })
m <- m[!grepl("^#", name) & nzchar(name)]
cat("rows read:", nrow(m), "->", paste(m$name, collapse = ","), "\n")
# with a comment line in the middle (a user edit)
an2 <- append(an, "# a later comment", after = 9)
writeLines(an2, file.path(S, "analyses2.tsv"))
m2 <- tryCatch(withCallingHandlers(fread(file.path(S, "analyses2.tsv"), header = FALSE, sep = "\t", col.names = c("name", "method", "model")),
                         warning = function(w) { cat("WARNING:", conditionMessage(w), "\n"); invokeRestart("muffleWarning") }), error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL })
if (!is.null(m2)) { m2 <- m2[!grepl("^#", name) & nzchar(name)]; cat("with mid-file comment, rows read:", nrow(m2), "->", paste(m2$name, collapse = ","), "\n") }

cat("\n=== E19: sprintf %d with whole doubles (profile_report.R) ===\n")
cat(sprintf("%d", 3), sprintf("%d", sum(c(1, 2))), "\n")
r <- tryCatch(sprintf("%d", 2.5), error = function(e) conditionMessage(e)); cat("2.5 ->", r, "\n")

cat("\n=== calibration ratio edge cases ===\n")
dist <- function(r) 1 - r
ratio <- c(dist(0.999) / dist(1), dist(1) / dist(1), dist(0.99) / dist(0.995))
verdict <- fifelse(!is.finite(ratio), "undetermined", fifelse(ratio <= 1.5, "within seed noise", "exceeds seed noise"))
print(data.table(ratio, verdict))
