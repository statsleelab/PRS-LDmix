#!/usr/bin/env Rscript
# Supplementary tables for the second-chromosome check: the chr2 runs of
# EUR+ASN 70:30 and EUR+AFR 50:50 next to the same settings on chr22, in the
# paper/<study> layout (LaTeX bodies, HTML preview, tables.tsv, frozen inputs,
# manifest).
#
#   table_chr2_accuracy.tex   S9a. r^2 of PRS-MIX and its paired difference
#                             from each single panel, chr22 against chr2.
#   table_chr2_overlap.tex    S9b. Top-10% overlap of PRS-MIX with each single
#                             panel minus the seed-rerun overlap, chr22 against
#                             chr2.
#   table_chr2_settings.tex   S9c. SNPs, causal SNPs and signal variance per run.
#   tables.html, tables.tsv, inputs/, manifest.json
#
# chr2 signal variance is scaled to the SNP count (v_g = 0.01 x chr2 SNPs /
# chr22 SNPs) so per-SNP effects match chr22. EUR+AFR is also run at chr22's
# v_g = 0.01: AFR's 205-person discovery sample adds chance correlations that
# grow with total v_g, so the two chr2 EUR+AFR runs bracket that effect.
# LaTeX and HTML show 250 kb blocks; tables.tsv also has chr2 at PRS-CS's EUR
# cut points. Cells are means over replicates with Monte Carlo SEs, paired
# within replicate, as in 18_sweep_paper_tables.R.
#
# Usage (from the PRS-MIX root):
#   Rscript inst/scripts/null_sim/24_chr2_tables.R <chr2_dir> <sweep_dir> <out_dir>
#
#   <chr2_dir>   e.g. inst/scripts/results/chr2_check, holding
#                <run>/<layout>/metrics/nonnull_metrics.tsv, <run>/nonnull/params.tsv
#                and <pair>/inputs/geno/analysis_snps.tsv
#   <sweep_dir>  e.g. inst/scripts/results/sim250_sweep_100rep (chr22)

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "common.R"))
args <- ns_args("Rscript 24_chr2_tables.R <chr2_dir> <sweep_dir> <out_dir>", 3)
chr2_dir <- args[[1]]
sweep_dir <- args[[2]]
out_dir <- args[[3]]
program <- "inst/scripts/null_sim/24_chr2_tables.R"
dir.create(file.path(out_dir, "inputs"), recursive = TRUE, showWarnings = FALSE)
to_paper <- grepl("(^|/)paper/", normalizePath(out_dir))

sha <- function(path) sub("\\s.*$", "", system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE))
git <- function(...) system2("git", c(...), stdout = TRUE)

# Each comparison: the chr22 sweep setting and the chr2 run(s).
runs <- list(
  list(id = "eur_asn_w0.7", A = "EUR", B = "ASN", w = 0.7, pair = "eur_asn", label = "EUR + ASN, 70:30",
       chr2 = c(scaled = "eur_asn_w0.7")),
  list(id = "eur_afr_w0.5", A = "EUR", B = "AFR", w = 0.5, pair = "eur_afr", label = "EUR + AFR, 50:50",
       chr2 = c(scaled = "eur_afr_w0.5", vg001 = "eur_afr_w0.5_vg0.01")))
chr2_layouts <- c(fixed250kb = "250 kb", eurcuts = "EUR cut points")

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
inputs <- list()
add <- function(key, src, dst) inputs[[key]] <<- freeze(src, dst)
for (r in runs) {
  add(paste0("chr22/", r$id, "/metrics"), file.path(sweep_dir, r$pair, sprintf("w%.1f", r$w), "metrics", "nonnull_metrics.tsv"),
      file.path("chr22", sprintf("%s_nonnull_metrics.tsv", r$id)))
  add(paste0("chr22/", r$id, "/params"), file.path(sweep_dir, r$pair, sprintf("w%.1f", r$w), "metrics", "nonnull_params.tsv"),
      file.path("chr22", sprintf("%s_nonnull_params.tsv", r$id)))
  add(paste0("chr22/", r$pair, "/snps"), file.path(sweep_dir, r$pair, "analysis_snps.tsv"), file.path("chr22", sprintf("%s_analysis_snps.tsv", r$pair)))
  add(paste0("chr2/", r$pair, "/snps"), file.path(chr2_dir, r$pair, "inputs", "geno", "analysis_snps.tsv"),
      file.path("chr2", sprintf("%s_analysis_snps.tsv", r$pair)))
  for (run in r$chr2) {
    add(paste0("chr2/", run, "/params"), file.path(chr2_dir, run, "nonnull", "params.tsv"), file.path("chr2", sprintf("%s_nonnull_params.tsv", run)))
    for (l in names(chr2_layouts)) add(paste0("chr2/", run, "/", l), file.path(chr2_dir, run, l, "metrics", "nonnull_metrics.tsv"),
                                       file.path("chr2", sprintf("%s_%s_nonnull_metrics.tsv", run, l)))
  }
}
inputs_dirty <- any(vapply(inputs, `[[`, TRUE, "source_uncommitted"))
prog_dirty <- length(git("status", "--porcelain", "--", program, "inst/scripts/null_sim/common.R")) > 0
if (to_paper && (inputs_dirty || prog_dirty)) {
  stop("Writing under paper/ needs committed inputs and program; commit them first.", call. = FALSE)
}
inp <- function(key) ns_read_tsv(inputs[[key]]$path)
prm <- function(key) { x <- inp(key); stats::setNames(x$value, x$param) }

# ---- Statistics ------------------------------------------------------------
cells <- list()
put <- function(table, label, column, chrom, target, x, layout = "fixed250kb") {
  x <- x[is.finite(x)]
  n <- length(x)
  cells[[length(cells) + 1]] <<- data.frame(table = table, setting = label, column = column, chrom = chrom, layout = layout,
                                            target = target, n = n, mean = if (n) mean(x) else NA_real_,
                                            se = if (n > 1) stats::sd(x) / sqrt(n) else NA_real_)
}
n_reps <- integer(0)
summarise <- function(x, r, chrom, column, layout) {
  x <- x[x$phi_mode == "auto", ]
  mix <- sprintf("mix_w%.1f", r$w)
  sa <- paste0("single_", r$A); sb <- paste0("single_", r$B)
  reps <- sort(Reduce(intersect, lapply(c(mix, sa, sb, paste0(sa, "_seed2"), paste0(sb, "_seed2")), function(a) x$rep[x$arm == a])))
  n_reps[[sprintf("%s|%s|%s|%s", r$id, chrom, column, layout)]] <<- length(reps)
  v <- function(arm, col) { y <- x[x$arm == arm, ]; y[[col]][match(reps, y$rep)] }
  for (tg in c(r$A, r$B)) {
    r2 <- function(a) v(a, paste0("r2_", tg))
    put("accuracy", r$label, column, chrom, tg, r2(mix), layout)
    cells[[length(cells)]]$column <<- paste0(column, "|r2_mix")
    put("accuracy", r$label, paste0(column, "|d_A"), chrom, tg, r2(mix) - r2(sa), layout)
    put("accuracy", r$label, paste0(column, "|d_B"), chrom, tg, r2(mix) - r2(sb), layout)
    for (ref in c("A", "B")) {
      rp <- r[[ref]]
      col <- sprintf("top10_%s_%s", rp, tg)
      put("overlap", r$label, paste0(column, "|gap_", ref), chrom, tg, 100 * (v(mix, col) - v(paste0("single_", rp, "_seed2"), col)), layout)
    }
  }
}
for (r in runs) {
  summarise(inp(paste0("chr22/", r$id, "/metrics")), r, "chr22", "chr22", "fixed250kb")
  for (k in names(r$chr2)) for (l in names(chr2_layouts)) {
    summarise(inp(paste0("chr2/", r$chr2[[k]], "/", l)), r, "chr2", paste0("chr2_", k), l)
  }
}
cells <- do.call(rbind, cells)
cells$quantity <- sub(".*\\|", "", cells$column)
cells$column <- sub("\\|.*", "", cells$column)
if (to_paper && any(n_reps < 100)) stop(sprintf("Fewest replicates is %d (%s).", min(n_reps), names(which.min(n_reps))), call. = FALSE)
ns_write_tsv(cells, file.path(out_dir, "tables.tsv"))

settings <- list()
for (r in runs) {
  p22 <- prm(paste0("chr22/", r$id, "/params"))
  n22 <- nrow(inp(paste0("chr22/", r$pair, "/snps")))
  settings[[length(settings) + 1]] <- data.frame(setting = r$label, run = "chr22", snps = n22, n_causal = p22[["n_causal"]], v_g = p22[["v_g"]])
  n2 <- nrow(inp(paste0("chr2/", r$pair, "/snps")))
  for (k in names(r$chr2)) {
    p2 <- prm(paste0("chr2/", r$chr2[[k]], "/params"))
    settings[[length(settings) + 1]] <- data.frame(setting = r$label, run = paste0("chr2_", k), snps = n2, n_causal = p2[["n_causal"]], v_g = p2[["v_g"]])
  }
}
settings <- do.call(rbind, settings)
ns_write_tsv(settings, file.path(out_dir, "settings.tsv"))

# ---- Rendering -------------------------------------------------------------
tx <- function(tex, t, h) if (tex) t else h
fmt <- function(x, d, signed, tex) {
  if (!is.finite(x)) return(tx(tex, "--", "&ndash;"))
  s <- formatC(abs(x), format = "f", digits = d)
  if (as.numeric(s) == 0) return(s)
  paste0(if (x < 0) tx(tex, "$-$", "&minus;") else if (signed) "+" else "", s)
}
cell <- function(table, label, column, target, quantity, d, signed, tex) {
  r <- cells[cells$table == table & cells$setting == label & cells$column == column & cells$layout == "fixed250kb" &
               cells$target == target & cells$quantity == quantity, ]
  if (nrow(r) == 0) return(tx(tex, "--", "&ndash;"))
  unit <- 10^-d
  se <- if (!is.finite(r$se)) "" else if (r$se < unit / 2) sprintf(" (%s%s)", tx(tex, "$<$", "&lt;"), formatC(unit, format = "f", digits = d))
        else sprintf(" (%s)", formatC(r$se, format = "f", digits = d))
  paste0(fmt(r$mean, d, signed, tex), se)
}
col_label <- function(column, tex) switch(column, chr22 = "chr22", chr2_scaled = tx(tex, "chr2, scaled $v_g$", "chr2, scaled <i>v</i><sub>g</sub>"),
                                          chr2_vg001 = tx(tex, "chr2, $v_g$ = 0.01", "chr2, <i>v</i><sub>g</sub> = 0.01"))
render <- function(heads, groups, tex) {
  nc <- length(heads)
  if (tex) {
    c(sprintf("%% Generated by %s; do not edit.", program), sprintf("\\begin{tabular}{l%s}", strrep("r", nc - 1)), "\\toprule",
      paste(paste(heads, collapse = " & "), "\\\\"), "\\midrule",
      unlist(lapply(seq_along(groups), function(g) c(
        if (g > 1) "\\addlinespace",
        sprintf("\\multicolumn{%d}{l}{\\textbf{%s}} \\\\", nc, names(groups)[g]),
        vapply(groups[[g]], function(r) paste(paste(r, collapse = " & "), "\\\\"), "")))),
      "\\bottomrule", "\\end{tabular}")
  } else {
    c("<table>", "<thead>", sprintf("<tr class=\"top\">%s</tr>", paste(sprintf("<th>%s</th>", heads), collapse = "")), "</thead>",
      unlist(lapply(seq_along(groups), function(g) c(
        "<tbody>", sprintf("<tr class=\"block\"><th colspan=\"%d\" scope=\"rowgroup\">%s</th></tr>", nc, names(groups)[g]),
        vapply(groups[[g]], function(r) sprintf("<tr>%s</tr>", paste(sprintf("<td>%s</td>", r), collapse = "")), ""), "</tbody>"))),
      "</table>")
  }
}
r2 <- function(tex) tx(tex, "$r^2$", "<i>r</i><sup>2</sup>")
build <- function(tex) {
  out <- list()
  groups <- list(); groups2 <- list()
  for (r in runs) {
    cols <- c("chr22", paste0("chr2_", names(r$chr2)))
    groups[[r$label]] <- unlist(lapply(cols, function(cn) lapply(c(r$A, r$B), function(tg) c(
      sprintf("%s, scored in %s", col_label(cn, tex), tg),
      cell("accuracy", r$label, cn, tg, "r2_mix", 3, FALSE, tex),
      cell("accuracy", r$label, cn, tg, "d_A", 3, TRUE, tex),
      cell("accuracy", r$label, cn, tg, "d_B", 3, TRUE, tex)))), recursive = FALSE)
    groups2[[r$label]] <- unlist(lapply(cols, function(cn) lapply(c(r$A, r$B), function(tg) c(
      sprintf("%s, scored in %s", col_label(cn, tex), tg),
      cell("overlap", r$label, cn, tg, "gap_A", 1, TRUE, tex),
      cell("overlap", r$label, cn, tg, "gap_B", 1, TRUE, tex)))), recursive = FALSE)
  }
  out$chr2_accuracy <- render(c("Chromosome, scored in", paste("PRS-MIX", r2(tex)), "PRS-MIX minus first panel", "PRS-MIX minus second panel"), groups, tex)
  out$chr2_overlap <- render(c("Chromosome, scored in", "Compared with the first panel", "Compared with the second panel"), groups2, tex)
  sg <- list()
  for (r in runs) {
    s <- settings[settings$setting == r$label, ]
    sg[[r$label]] <- lapply(seq_len(nrow(s)), function(i) c(col_label(s$run[i], tex), formatC(s$snps[i], format = "d", big.mark = ","),
                                                            formatC(s$n_causal[i], format = "d", big.mark = ","), format(signif(s$v_g[i], 3))))
  }
  out$chr2_settings <- render(c("Run", "Analysis SNPs", "Causal SNPs", tx(tex, "$v_g$", "<i>v</i><sub>g</sub>")), sg, tex)
  out
}
tex_tabs <- build(TRUE)
for (nm in names(tex_tabs)) writeLines(tex_tabs[[nm]], file.path(out_dir, sprintf("table_%s.tex", nm)))
html_tabs <- build(FALSE)
titles <- c(chr2_accuracy = "Table S9a. Accuracy on chromosome 22 and chromosome 2",
            chr2_overlap = "Table S9b. Change in the top 10% beyond rerun randomness, chromosome 22 and chromosome 2",
            chr2_settings = "Table S9c. Second-chromosome run settings")
cap_common <- paste("Fixed 250 kb blocks, 100 replicates per run, PRS-CS defaults with <i>&phi;</i> learned. First panel: EUR;",
                    "second panel: ASN or AFR. Means over replicates with Monte Carlo standard errors in parentheses; differences are paired",
                    "within replicate. Scaled <i>v</i><sub>g</sub>: signal variance scaled to the SNP count so per-SNP effects match chromosome 22.")
captions <- c(
  chr2_accuracy = paste("<i>r</i><sup>2</sup> of PRS-MIX with the true genetic value, and PRS-MIX minus each single panel.", cap_common),
  chr2_overlap = paste("PRS-MIX's top-10% overlap with each single panel minus the overlap of that panel's rerun with a new seed,",
                       "in percentage points; more negative means PRS-MIX moves the top 10% more than sampler randomness.", cap_common),
  chr2_settings = paste("Chromosome 22 runs are the sweep settings; chromosome 2 runs use the same individuals. EUR + AFR is also run on",
                        "chromosome 2 at chromosome 22's <i>v</i><sub>g</sub> = 0.01: chance correlations in the 205-person AFR discovery",
                        "sample inflate the GWAS signal in proportion to total <i>v</i><sub>g</sub>, so the two chromosome 2 runs bracket that effect."))
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
          "</style>", "</head>", "<body>", "<main>", "<h1>Second-chromosome check: supplementary tables</h1>",
          unlist(lapply(names(html_tabs), function(nm) c("<section>", sprintf("<h2>%s</h2>", titles[[nm]]),
                                                         sprintf("<p class=\"cap\">%s</p>", captions[[nm]]),
                                                         "<div class=\"wrap\">", html_tabs[[nm]], "</div>", "</section>"))),
          sprintf("<p class=\"src\">Generated by %s at PRS-MIX %s%s. tables.tsv also has chromosome 2 at PRS-CS's EUR cut points.</p>",
                  program, substr(head_commit, 1, 7), if (inputs_dirty || prog_dirty) " with uncommitted program or inputs" else ""),
          "</main>", "</body>", "</html>")
writeLines(html, file.path(out_dir, "tables.html"))

out_files <- c(sprintf("table_%s.tex", names(tex_tabs)), "tables.html", "tables.tsv", "settings.tsv")
manifest <- list(
  study = basename(normalizePath(out_dir)), program = program,
  command = paste("Rscript", program, paste(shQuote(args), collapse = " ")),
  asset_commit = head_commit, program_uncommitted = prog_dirty, inputs_uncommitted = inputs_dirty,
  experiment_commit = unique(stats::na.omit(vapply(inputs, `[[`, "", "source_commit"))),
  r_version = R.version.string,
  inputs = unname(lapply(inputs, function(x) x[c("path", "source", "source_commit", "sha256")])),
  outputs = lapply(stats::setNames(out_files, out_files), function(f) list(path = file.path(out_dir, f), sha256 = sha(file.path(out_dir, f)))),
  settings = list(comparisons = unname(vapply(runs, `[[`, "", "label")), replicates_used = as.list(n_reps), phi_mode = "auto",
                  chr2_vg = "scaled: 0.01 x chr2 analysis SNPs / chr22 analysis SNPs of the pair; EUR+AFR also at 0.01",
                  uncertainty = "mean over replicates with Monte Carlo SE; paired differences use the SE of the per-replicate difference"),
  notes = "chr2 runs were made on RedHawk with cluster/jobs/chr2_check (not version-controlled), on the chr22 sweep's sample split.")
jsonlite::write_json(manifest, file.path(out_dir, "manifest.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
message(sprintf("Wrote %s (%d cells) at PRS-MIX %s%s", out_dir, nrow(cells), substr(head_commit, 1, 7),
                if (inputs_dirty || prog_dirty) " with uncommitted program or inputs" else ""))
