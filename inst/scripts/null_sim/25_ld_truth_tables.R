#!/usr/bin/env Rscript
# Supplementary table for the true-LD comparison: chromosome-wide summaries of
# the per-block scans behind the block-75 figures (21_ld_truth_figures.R), in
# the paper/<study> layout (LaTeX body, HTML preview, tables.tsv, frozen
# inputs, manifest).
#
#   table_ld_chr.tex   For each pair and GWAS composition: the number of blocks,
#                      the median over blocks of the mean |difference| from the
#                      true LD for PRS-LDmix, the EUR panel and the second panel,
#                      and the number of blocks where PRS-LDmix is closer to the
#                      true LD than the EUR panel.
#   tables.html, tables.tsv, inputs/, manifest.json
#
# The true LD is the discovery LD blended at the GWAS composition; no PRS-CS
# output and no replicates. Blocks are those in <pair>_block_scan.tsv (at least
# two analysis SNPs).
#
# Usage (from the PRS-LDmix root):
#   Rscript inst/scripts/null_sim/25_ld_truth_tables.R <scan_dir> <out_dir>
#
#   <scan_dir>  e.g. paper/null_sim_ld_truth_250kb/figures/inputs, holding
#               <pair>_block_scan.tsv and <pair>_params.tsv

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "common.R"))
args <- ns_args("Rscript 25_ld_truth_tables.R <scan_dir> <out_dir>", 2)
scan_dir <- args[[1]]
out_dir <- args[[2]]
program <- "inst/scripts/null_sim/25_ld_truth_tables.R"
dir.create(file.path(out_dir, "inputs"), recursive = TRUE, showWarnings = FALSE)
to_paper <- grepl("(^|/)paper/", normalizePath(out_dir))

sha <- function(path) sub("\\s.*$", "", system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE))
git <- function(...) system2("git", c(...), stdout = TRUE)

pairs <- list(
  list(id = "eur_asn", A = "EUR", B = "ASN", label = "EUR and ASN"),
  list(id = "eur_afr", A = "EUR", B = "AFR", label = "EUR and AFR"),
  list(id = "eur_sas", A = "EUR", B = "SAS", label = "EUR and SAS"))
weights <- c(0.9, 0.7, 0.5)

# ---- Freeze inputs -------------------------------------------------------
freeze <- function(src) {
  dst <- file.path(out_dir, "inputs", basename(src))
  file.copy(src, dst, overwrite = TRUE)
  stopifnot(sha(src) == sha(dst))
  list(path = dst, source = src, sha256 = sha(dst),
       source_commit = c(git("log", "-1", "--format=%H", "--", src), NA_character_)[1],
       source_uncommitted = length(git("status", "--porcelain", "--", src)) > 0)
}
inputs <- list()
for (p in pairs) for (f in c("block_scan", "params")) {
  inputs[[paste0(p$id, "/", f)]] <- freeze(file.path(scan_dir, sprintf("%s_%s.tsv", p$id, f)))
}
inputs_dirty <- any(vapply(inputs, `[[`, TRUE, "source_uncommitted"))
prog_dirty <- length(git("status", "--porcelain", "--", program, "inst/scripts/null_sim/common.R")) > 0
if (to_paper && (inputs_dirty || prog_dirty)) {
  stop("Writing under paper/ needs committed inputs and program; commit them first.", call. = FALSE)
}

# ---- Statistics ------------------------------------------------------------
cells <- list()
put <- function(p, w, column, value) {
  cells[[length(cells) + 1]] <<- data.frame(setting = sprintf("%s, %d:%d", p$label, round(100 * w), round(100 * (1 - w))),
                                            pair = p$id, w = w, column = column, value = value)
}
for (p in pairs) {
  scan <- ns_read_tsv(inputs[[paste0(p$id, "/block_scan")]]$path)
  # Scans extracted before the rename to PRS-LDmix name the column mad_prsmix.
  names(scan)[names(scan) == "mad_prsmix"] <- "mad_prsldmix"
  for (w in weights) {
    s <- scan[abs(scan$w - w) < 1e-9, ]
    stopifnot(nrow(s) > 0, !anyDuplicated(s$block))
    put(p, w, "blocks", nrow(s))
    put(p, w, "median_mad_prsldmix", stats::median(s$mad_prsldmix))
    put(p, w, "median_mad_a_only", stats::median(s$mad_a_only))
    put(p, w, "median_mad_b_only", stats::median(s$mad_b_only))
    put(p, w, "prsldmix_closer_than_a", sum(s$mad_prsldmix < s$mad_a_only))
  }
}
cells <- do.call(rbind, cells)
ns_write_tsv(cells, file.path(out_dir, "tables.tsv"))

# ---- Rendering -------------------------------------------------------------
tx <- function(tex, t, h) if (tex) t else h
val <- function(p, w, column) cells$value[cells$pair == p$id & abs(cells$w - w) < 1e-9 & cells$column == column]
build <- function(tex) {
  heads <- c(tx(tex, "GWAS ancestry mix", "GWAS ancestry mix"), "Blocks", "PRS-LDmix", "EUR panel", "Second panel",
             "PRS-LDmix closer than EUR panel")
  rows <- lapply(pairs, function(p) lapply(weights, function(w) c(
    tx(tex, sprintf("%d\\%% EUR, %d\\%% %s", round(100 * w), round(100 * (1 - w)), p$B),
       sprintf("%d%% EUR, %d%% %s", round(100 * w), round(100 * (1 - w)), p$B)),
    formatC(val(p, w, "blocks"), format = "d"),
    formatC(val(p, w, "median_mad_prsldmix"), format = "f", digits = 3),
    formatC(val(p, w, "median_mad_a_only"), format = "f", digits = 3),
    formatC(val(p, w, "median_mad_b_only"), format = "f", digits = 3),
    sprintf("%d of %d", as.integer(val(p, w, "prsldmix_closer_than_a")), as.integer(val(p, w, "blocks"))))))
  names(rows) <- vapply(pairs, `[[`, "", "label")
  nc <- length(heads)
  if (tex) {
    c(sprintf("%% Generated by %s; do not edit.", program), "\\begin{tabular}{lrrrrr}", "\\toprule",
      " &  & \\multicolumn{3}{c}{Median over blocks of the mean $|$difference$|$ from the true LD} &  \\\\",
      "\\cmidrule(lr){3-5}",
      paste(paste(heads, collapse = " & "), "\\\\"), "\\midrule",
      unlist(lapply(seq_along(rows), function(g) c(
        if (g > 1) "\\addlinespace",
        sprintf("\\multicolumn{%d}{l}{\\textbf{%s}} \\\\", nc, names(rows)[g]),
        vapply(rows[[g]], function(r) paste(paste(r, collapse = " & "), "\\\\"), "")))),
      "\\bottomrule", "\\end{tabular}")
  } else {
    c("<table>", "<thead>",
      "<tr class=\"top\"><th></th><th></th><th colspan=\"3\">Median over blocks of the mean |difference| from the true LD</th><th></th></tr>",
      sprintf("<tr class=\"top\">%s</tr>", paste(sprintf("<th>%s</th>", heads), collapse = "")), "</thead>",
      unlist(lapply(seq_along(rows), function(g) c(
        "<tbody>", sprintf("<tr class=\"block\"><th colspan=\"%d\" scope=\"rowgroup\">%s</th></tr>", nc, names(rows)[g]),
        vapply(rows[[g]], function(r) sprintf("<tr>%s</tr>", paste(sprintf("<td>%s</td>", r), collapse = "")), ""), "</tbody>"))),
      "</table>")
  }
}
writeLines(build(TRUE), file.path(out_dir, "table_ld_chr.tex"))
caption <- paste("Chromosome 22, fixed 250 kb blocks with at least two analysis SNPs. For each block, the mean |difference| over",
                 "off-diagonal SNP pairs between the true LD (the discovery LD blended at the GWAS composition) and the LD each",
                 "analysis assumes: PRS-LDmix's blend of the reference panels at the same composition, the EUR panel, or the",
                 "second ancestry's panel. Medians over blocks, and the number of blocks where PRS-LDmix is closer to the true LD than",
                 "the EUR panel. No PRS-CS output and no replicates.")
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
          "</style>", "</head>", "<body>", "<main>", "<h1>True-LD comparison: supplementary table</h1>",
          "<section>", "<h2>Chromosome-wide LD difference from the true GWAS LD</h2>",
          sprintf("<p class=\"cap\">%s</p>", caption), "<div class=\"wrap\">", build(FALSE), "</div>", "</section>",
          sprintf("<p class=\"src\">Generated by %s at PRS-LDmix %s%s.</p>",
                  program, substr(head_commit, 1, 7), if (inputs_dirty || prog_dirty) " with uncommitted program or inputs" else ""),
          "</main>", "</body>", "</html>")
writeLines(html, file.path(out_dir, "tables.html"))

out_files <- c("table_ld_chr.tex", "tables.html", "tables.tsv")
manifest <- list(
  study = basename(normalizePath(out_dir)), program = program,
  command = paste("Rscript", program, paste(shQuote(args), collapse = " ")),
  asset_commit = head_commit, program_uncommitted = prog_dirty, inputs_uncommitted = inputs_dirty,
  experiment_commit = unique(stats::na.omit(vapply(inputs, `[[`, "", "source_commit"))),
  r_version = R.version.string,
  inputs = unname(lapply(inputs, function(x) x[c("path", "source", "source_commit", "sha256")])),
  outputs = lapply(stats::setNames(out_files, out_files), function(f) list(path = file.path(out_dir, f), sha256 = sha(file.path(out_dir, f)))),
  settings = list(chromosome = 22, block_layout = "fixed 250 kb", pairs = vapply(pairs, `[[`, "", "id"), compositions_A = weights,
                  blocks = "blocks in <pair>_block_scan.tsv (at least two analysis SNPs)",
                  summary = "median over blocks of mean |true - candidate| over off-diagonal SNP pairs; count of blocks with mad_prsldmix < mad_a_only"),
  notes = "Inputs are the frozen per-block scans of paper/null_sim_ld_truth_250kb/figures (RedHawk extracts from 20_ld_truth_extract.R).")
jsonlite::write_json(manifest, file.path(out_dir, "manifest.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
message(sprintf("Wrote %s (%d cells) at PRS-LDmix %s%s", out_dir, nrow(cells), substr(head_commit, 1, 7),
                if (inputs_dirty || prog_dirty) " with uncommitted program or inputs" else ""))
