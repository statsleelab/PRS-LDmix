#!/usr/bin/env Rscript
# Supplementary tables validating the genome-wide PRS-MIX reference panels, in
# the paper/<study> layout: LaTeX bodies, an HTML preview, every number in
# tables.tsv, frozen inputs and a manifest.
#
#   table_panel_checks.tex     S1a. Per layout and superpopulation: chromosomes
#                              passing the structural and matrix checks, blocks
#                              and SNPs, and agreement of within-block r with
#                              PRS-CS's published 1000 Genomes panels (median
#                              and range over chromosomes, shared SNP pairs).
#   table_panel_samples.tex    S1b. Per superpopulation: reference sample size,
#                              mean r^2 between neighbouring SNPs, median MAF and
#                              share of SNPs with MAF < 0.05 (250 kb panels).
#   table_upstream.tex         S1c. PRS-MIX's PRS-CS script with a single panel
#                              against unmodified PRS-CS: posterior-effect files
#                              compared and identical.
#   tables.html, tables.tsv, inputs/, manifest.json
#
# Inputs are the per-chromosome reports of validate_prscs_panels.R
# (<layout_dir>/chr<k>.txt and summary.tsv) and the upstream check table.
#
# Usage (from the PRS-MIX root):
#   Rscript inst/scripts/null_sim/23_panel_validation_tables.R <results_dir> <upstream_tsv> <out_dir>
#
#   <results_dir>   e.g. inst/scripts/results, holding panel_validation_fixed250kb/
#                   and panel_validation_eurcuts/
#   <upstream_tsv>  e.g. inst/scripts/results/null_sim_250kb_100rep/upstream_check.tsv

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "common.R"))
args <- ns_args("Rscript 23_panel_validation_tables.R <results_dir> <upstream_tsv> <out_dir>", 3)
results_dir <- args[[1]]
upstream_tsv <- args[[2]]
out_dir <- args[[3]]
program <- "inst/scripts/null_sim/23_panel_validation_tables.R"
dir.create(file.path(out_dir, "inputs"), recursive = TRUE, showWarnings = FALSE)
to_paper <- grepl("(^|/)paper/", normalizePath(out_dir))

sha <- function(path) sub("\\s.*$", "", system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE))
git <- function(...) system2("git", c(...), stdout = TRUE)

layouts <- c(fixed250kb = "Fixed 250 kb blocks", eurcuts = "PRS-CS EUR cut points")
pops <- c("EUR", "ASN", "AFR", "AMR", "SAS")
prscs_ref <- c(EUR = "EUR", ASN = "EAS", AFR = "AFR")
chrs <- 1:22

# ---- Freeze inputs -------------------------------------------------------
freeze <- function(src, dst_rel) {
  dst <- file.path(out_dir, "inputs", dst_rel)
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
  file.copy(src, dst, overwrite = TRUE)
  stopifnot(sha(src) == sha(dst))
  list(path = dst, source = src, sha256 = sha(dst),
       source_commit = c(git("log", "-1", "--format=%H", "--", src), NA_character_)[1],
       source_uncommitted = length(git("status", "--porcelain", "--", src)) > 0)
}
rel <- c(unlist(lapply(names(layouts), function(l) file.path(paste0("panel_validation_", l), c(sprintf("chr%d.txt", chrs), "summary.tsv")))))
inputs <- c(lapply(rel, function(f) freeze(file.path(results_dir, f), f)), list(freeze(upstream_tsv, "upstream_check.tsv")))
names(inputs) <- c(rel, "upstream_check.tsv")
inputs_dirty <- any(vapply(inputs, `[[`, TRUE, "source_uncommitted"))
prog_dirty <- length(git("status", "--porcelain", "--", program, "inst/scripts/null_sim/common.R")) > 0
if (to_paper && (inputs_dirty || prog_dirty)) {
  stop("Writing under paper/ needs committed inputs and program; commit them first.", call. = FALSE)
}

# ---- Parse the per-chromosome reports ------------------------------------
num <- function(s) as.numeric(gsub(",", "", s))
reports <- list()
for (l in names(layouts)) for (k in chrs) {
  txt <- readLines(inputs[[file.path(paste0("panel_validation_", l), sprintf("chr%d.txt", k))]]$path, warn = FALSE)
  fails <- sum(grepl("\\[FAIL\\]", txt))
  overall <- grepl("Overall structural and matrix checks: PASS", txt)
  for (p in pops) {
    st <- regmatches(txt, regexec(sprintf("^\\s+\\[PASS\\] %s: (\\d+) blocks list the same ([0-9,]+) SNPs", p), txt))
    st <- st[lengths(st) > 0]
    ag <- regmatches(txt, regexec(sprintf("\\[INFO\\] %s vs PRS-CS 1KG \\w+: r correlation ([0-9.]+), RMSE ([0-9.]+) over ([0-9,]+) shared within-block pairs \\(([0-9,]+) SNPs in both\\)", p), txt))
    ag <- ag[lengths(ag) > 0]
    pc <- regmatches(txt, regexec(sprintf("^\\s+%s\\s+([0-9,]+)\\s+([0-9.]+)\\s+([0-9.]+)\\s+([0-9.]+)%%", p), txt))
    pc <- pc[lengths(pc) > 0]
    reports[[length(reports) + 1]] <- data.frame(
      layout = l, chr = k, pop = p, checks_pass = any(overall) && fails == 0,
      blocks = if (length(st)) num(st[[1]][2]) else NA, snps = if (length(st)) num(st[[1]][3]) else NA,
      cor = if (length(ag)) num(ag[[1]][2]) else NA, rmse = if (length(ag)) num(ag[[1]][3]) else NA,
      pairs = if (length(ag)) num(ag[[1]][4]) else NA, snps_shared = if (length(ag)) num(ag[[1]][5]) else NA,
      n = if (length(pc)) num(pc[[1]][2]) else NA, r2_neighbour = if (length(pc)) num(pc[[1]][3]) else NA,
      maf_median = if (length(pc)) num(pc[[1]][4]) else NA, maf_lt05 = if (length(pc)) num(pc[[1]][5]) else NA)
  }
}
rep_df <- do.call(rbind, reports)
# Cross-check against each layout's summary.tsv (written by the same validation runs).
for (l in names(layouts)) {
  sm <- ns_read_tsv(inputs[[file.path(paste0("panel_validation_", l), "summary.tsv")]]$path)
  for (p in names(prscs_ref)) {
    a <- rep_df$cor[rep_df$layout == l & rep_df$pop == p][match(sm$chr, chrs)]
    stopifnot(isTRUE(all.equal(a, sm[[sprintf("%s_vs_PRScs_%s_cor", p, prscs_ref[[p]])]], tolerance = 1e-9)))
  }
  stopifnot(all(sm$checks == "PASS"))
}

# ---- Statistics ------------------------------------------------------------
cells <- list()
put <- function(table, layout, row, column, value) {
  cells[[length(cells) + 1]] <<- data.frame(table = table, layout = layout, row = row, column = column, value = value)
}
for (l in names(layouts)) for (p in pops) {
  x <- rep_df[rep_df$layout == l & rep_df$pop == p, ]
  put("checks", l, p, "chromosomes_pass", sum(x$checks_pass))
  put("checks", l, p, "blocks", sum(x$blocks))
  put("checks", l, p, "snps", sum(x$snps))
  if (p %in% names(prscs_ref)) {
    put("checks", l, p, "cor_median", stats::median(x$cor)); put("checks", l, p, "cor_min", min(x$cor)); put("checks", l, p, "cor_max", max(x$cor))
    put("checks", l, p, "rmse_median", stats::median(x$rmse)); put("checks", l, p, "rmse_min", min(x$rmse)); put("checks", l, p, "rmse_max", max(x$rmse))
    put("checks", l, p, "pairs", sum(x$pairs))
  }
}
# S1b: per superpopulation, 250 kb panels, chromosome values weighted by SNP count.
for (p in pops) {
  x <- rep_df[rep_df$layout == "fixed250kb" & rep_df$pop == p, ]
  stopifnot(length(unique(x$n)) == 1)
  wm <- function(v) sum(v * x$snps) / sum(x$snps)
  put("samples", "fixed250kb", p, "n", x$n[1])
  put("samples", "fixed250kb", p, "r2_neighbour", wm(x$r2_neighbour))
  put("samples", "fixed250kb", p, "maf_median", wm(x$maf_median))
  put("samples", "fixed250kb", p, "maf_lt05", wm(x$maf_lt05))
}
up <- ns_read_tsv(inputs[["upstream_check.tsv"]]$path)
for (kd in unique(up$kind)) for (ph in unique(up$phi)) {
  x <- up[up$kind == kd & up$phi == ph, ]
  row <- sprintf("%s|%s", kd, ph)
  put("upstream", "", row, "files", nrow(x))
  put("upstream", "", row, "identical", sum(x$bytes == "identical"))
  put("upstream", "", row, "max_abs_diff", max(x$max_abs_diff))
  put("upstream", "", row, "snp_mismatches", sum(x$snp_mismatches))
  put("upstream", "", row, "replicates", length(unique(x$rep)))
}
cells <- do.call(rbind, cells)
ns_write_tsv(cells, file.path(out_dir, "tables.tsv"))
ns_write_tsv(rep_df, file.path(out_dir, "per_chromosome.tsv"))

# ---- Rendering -------------------------------------------------------------
tx <- function(tex, t, h) if (tex) t else h
v <- function(table, layout, row, column) { r <- cells$value[cells$table == table & cells$layout == layout & cells$row == row & cells$column == column]; if (length(r)) r else NA }
big <- function(x, tex) formatC(x, format = "d", big.mark = ",")
f3 <- function(x) formatC(x, format = "f", digits = 3)
dash <- function(tex) tx(tex, "--", "&ndash;")
render <- function(heads, groups, tex) {
  nc <- length(heads)
  if (tex) {
    c(sprintf("%% Generated by %s; do not edit.", program), sprintf("\\begin{tabular}{l%s}", strrep("r", nc - 1)), "\\toprule",
      paste(paste(heads, collapse = " & "), "\\\\"), "\\midrule",
      unlist(lapply(seq_along(groups), function(g) c(
        if (g > 1) "\\addlinespace",
        if (nzchar(names(groups)[g])) sprintf("\\multicolumn{%d}{l}{\\textbf{%s}} \\\\", nc, names(groups)[g]),
        vapply(groups[[g]], function(r) paste(paste(r, collapse = " & "), "\\\\"), "")))),
      "\\bottomrule", "\\end{tabular}")
  } else {
    c("<table>", "<thead>", sprintf("<tr class=\"top\">%s</tr>", paste(sprintf("<th>%s</th>", heads), collapse = "")), "</thead>",
      unlist(lapply(seq_along(groups), function(g) c(
        "<tbody>", if (nzchar(names(groups)[g])) sprintf("<tr class=\"block\"><th colspan=\"%d\" scope=\"rowgroup\">%s</th></tr>", nc, names(groups)[g]),
        vapply(groups[[g]], function(r) sprintf("<tr>%s</tr>", paste(sprintf("<td>%s</td>", r), collapse = "")), ""), "</tbody>"))),
      "</table>")
  }
}
build <- function(tex) {
  out <- list()
  pct <- tx(tex, "\\%", "%")
  rng <- function(l, p, m) if (is.na(v("checks", l, p, paste0(m, "_median")))) dash(tex) else
    sprintf("%s [%s, %s]", f3(v("checks", l, p, paste0(m, "_median"))), f3(v("checks", l, p, paste0(m, "_min"))), f3(v("checks", l, p, paste0(m, "_max"))))
  heads <- c("Superpopulation", "Chromosomes passing", "Blocks", "SNPs", "PRS-CS 1KG panel",
             "Correlation of $r$", "RMSE of $r$", "Shared SNP pairs")
  if (!tex) heads[6:7] <- c("Correlation of <i>r</i>", "RMSE of <i>r</i>")
  groups <- lapply(names(layouts), function(l) lapply(pops, function(p) c(
    p, sprintf("%d of 22", v("checks", l, p, "chromosomes_pass")), big(v("checks", l, p, "blocks")), big(v("checks", l, p, "snps")),
    if (p %in% names(prscs_ref)) prscs_ref[[p]] else tx(tex, "none", "none"),
    rng(l, p, "cor"), rng(l, p, "rmse"),
    if (p %in% names(prscs_ref)) big(v("checks", l, p, "pairs")) else dash(tex))))
  names(groups) <- unname(layouts)
  out$panel_checks <- render(heads, groups, tex)
  heads <- c("Superpopulation", "Reference individuals",
             tx(tex, "Mean $r^2$, neighbouring SNPs", "Mean <i>r</i><sup>2</sup>, neighbouring SNPs"), "Median MAF",
             paste0("SNPs with MAF < 0.05 (", pct, ")"))
  if (tex) heads[5] <- paste0("SNPs with MAF $<$ 0.05 (", pct, ")")
  out$panel_samples <- render(heads, list(" " = lapply(pops, function(p) c(
    p, big(v("samples", "fixed250kb", p, "n")), f3(v("samples", "fixed250kb", p, "r2_neighbour")),
    f3(v("samples", "fixed250kb", p, "maf_median")), formatC(v("samples", "fixed250kb", p, "maf_lt05"), format = "f", digits = 1)))), tex)
  out$panel_samples <- sub("\\\\multicolumn\\{5\\}\\{l\\}\\{\\\\textbf\\{ \\}\\} \\\\\\\\", "", out$panel_samples)
  rows <- unique(cells$row[cells$table == "upstream"])
  heads <- c("GWAS", tx(tex, "PRS-CS $\\phi$", "PRS-CS <i>&phi;</i>"), "Replicates", "Files compared", "Byte-identical", tx(tex, "Max $|$difference$|$", "Max |difference|"))
  out$upstream <- render(heads, list(" " = lapply(rows, function(r) {
    kp <- strsplit(r, "|", fixed = TRUE)[[1]]
    c(if (kp[1] == "null") "Null" else "Non-null", if (kp[2] == "auto") "learned" else format(as.numeric(kp[2]), scientific = FALSE), v("upstream", "", r, "replicates"), v("upstream", "", r, "files"),
      v("upstream", "", r, "identical"), format(v("upstream", "", r, "max_abs_diff")))
  })), tex)
  out$upstream <- sub("\\\\multicolumn\\{6\\}\\{l\\}\\{\\\\textbf\\{ \\}\\} \\\\\\\\", "", out$upstream)
  out
}
tex_tabs <- build(TRUE)
for (nm in names(tex_tabs)) writeLines(tex_tabs[[nm]][nzchar(tex_tabs[[nm]])], file.path(out_dir, sprintf("table_%s.tex", nm)))
html_tabs <- build(FALSE)
html_tabs <- lapply(html_tabs, function(t) t[!grepl(">\\s<", t) | !grepl("rowgroup", t)])

titles <- c(panel_checks = "Table S1a. Panel checks and agreement with PRS-CS's published panels",
            panel_samples = "Table S1b. Reference individuals and LD features by superpopulation",
            upstream = "Table S1c. PRS-MIX with a single panel reproduces PRS-CS exactly")
captions <- c(
  panel_checks = paste("Genome-wide PRS-MIX panels (chromosomes 1&ndash;22, HapMap3 SNPs) built from the 33KG reference panel for",
                       "five superpopulations, at both block layouts. Chromosomes passing: every panel lists the same SNPs and",
                       "alleles in the same order, and every block's LD matrix is finite, symmetric, has a unit diagonal and",
                       "|<i>r</i>| &le; 1, and is positive definite. Agreement: correlation and root-mean-square difference between",
                       "our within-block <i>r</i> and that of PRS-CS's published 1000 Genomes panel for the matching population, over",
                       "SNP pairs in the same block in both panels, after aligning alleles; median and range over the 22 chromosomes.",
                       "PRS-CS publishes no AMR or SAS panel."),
  panel_samples = paste("Fixed 250 kb panels. Mean <i>r</i><sup>2</sup> between neighbouring SNPs, median minor allele frequency and the",
                        "share of SNPs with MAF &lt; 0.05, averaged over chromosomes weighted by SNP count. AFR shows the weaker",
                        "short-range LD expected of African-ancestry populations."),
  upstream = paste("Chromosome 22, 250 kb panels, EUR+ASN 70:30 GWAS, replicates 1&ndash;10. PRS-MIX's modified PRS-CS script",
                   "(<code>PRScs_mix.py</code>) given a single panel, against unmodified PRS-CS (<code>PRScs.py</code>, commit 5330390)",
                   "with the same panel, GWAS and seed, for both single panels. Files compared: posterior SNP-effect files."))
sans <- "-apple-system, \"Helvetica Neue\", Helvetica, Arial, sans-serif"
head_commit <- git("rev-parse", "HEAD")
html <- c("<!doctype html>", "<html lang=\"en\">", "<head>", "<meta charset=\"utf-8\">",
          "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
          sprintf("<title>%s tables</title>", basename(normalizePath(out_dir))),
          sprintf("<!-- Generated by %s; do not edit. -->", program), "<style>",
          "body { margin: 0; padding: 40px 24px 56px; background: #fff; color: #111; font: 15px/1.45 Georgia, \"Times New Roman\", Times, serif; }",
          "main { max-width: 1180px; margin: 0 auto; display: flex; flex-direction: column; gap: 44px; }",
          "h1 { font-size: 22px; margin: 0; } h2 { font-size: 17px; margin: 0 0 6px; }",
          "p.cap { margin: 0 0 12px; color: #333; font-size: 14px; max-width: 980px; }",
          ".wrap { overflow-x: auto; }",
          sprintf("table { border-collapse: collapse; font-family: %s; font-size: 13px; font-variant-numeric: tabular-nums; }", sans),
          "th, td { padding: 4px 12px; text-align: right; white-space: nowrap; } th:first-child, td:first-child { text-align: left; }",
          "thead tr.top th { border-top: 2px solid #111; border-bottom: 1px solid #111; font-weight: 600; }",
          "tr.block th { text-align: left; padding-top: 12px; font-weight: 700; }",
          "tbody:last-child tr:last-child td { border-bottom: 2px solid #111; }",
          sprintf("p.src { color: #666; font: 12px %s; }", sans),
          "</style>", "</head>", "<body>", "<main>",
          "<h1>Reference-panel validation: supplementary tables</h1>",
          unlist(lapply(names(html_tabs), function(nm) c("<section>", sprintf("<h2>%s</h2>", titles[[nm]]),
                                                         sprintf("<p class=\"cap\">%s</p>", captions[[nm]]),
                                                         "<div class=\"wrap\">", html_tabs[[nm]], "</div>", "</section>"))),
          sprintf("<p class=\"src\">Generated by %s at PRS-MIX %s%s.</p>", program, substr(head_commit, 1, 7),
                  if (inputs_dirty || prog_dirty) " with uncommitted program or inputs" else ""),
          "</main>", "</body>", "</html>")
writeLines(html, file.path(out_dir, "tables.html"))

# ---- Manifest ------------------------------------------------------------
out_files <- c(sprintf("table_%s.tex", names(tex_tabs)), "tables.html", "tables.tsv", "per_chromosome.tsv")
manifest <- list(
  study = basename(normalizePath(out_dir)), program = program,
  command = paste("Rscript", program, paste(shQuote(args), collapse = " ")),
  asset_commit = head_commit, program_uncommitted = prog_dirty, inputs_uncommitted = inputs_dirty,
  experiment_commit = unique(stats::na.omit(vapply(inputs, `[[`, "", "source_commit"))),
  r_version = R.version.string,
  inputs = unname(lapply(inputs, function(x) x[c("path", "source", "source_commit", "sha256")])),
  outputs = lapply(stats::setNames(out_files, out_files), function(f) list(path = file.path(out_dir, f), sha256 = sha(file.path(out_dir, f)))),
  settings = list(layouts = unname(layouts), superpopulations = pops, prscs_reference = as.list(prscs_ref), chromosomes = "1-22",
                  validator = "inst/scripts/validate_prscs_panels.R",
                  upstream = "PRScs_mix.py with one --ref_dir against PRScs.py (PRS-CS 5330390), same inputs and seed"),
  notes = "Validation and upstream runs were made on RedHawk (cluster/jobs/validate and null_sim/upstream_check.sbatch, not version-controlled).")
jsonlite::write_json(manifest, file.path(out_dir, "manifest.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
message(sprintf("Wrote %s (%d cells) at PRS-MIX %s%s", out_dir, nrow(cells), substr(head_commit, 1, 7),
                if (inputs_dirty || prog_dirty) " with uncommitted program or inputs" else ""))
