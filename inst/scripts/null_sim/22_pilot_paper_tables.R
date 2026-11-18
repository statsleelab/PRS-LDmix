#!/usr/bin/env Rscript
# Supplementary tables from the prespecified layout pilot (LAYOUT_PILOT.md) and
# the 100-replicate chain check, in the paper/<study> layout: LaTeX bodies, an
# HTML preview, every number in tables.tsv, frozen inputs and a manifest.
#
#   table_chain_pilot.tex      S6a. Chain length in the pilot (replicates 1-5):
#                              r^2 at the defaults, 4,000 minus default and
#                              8,000 minus 4,000 iterations, the seed-rerun
#                              spread and the learned phi, per arm and target.
#   table_chain_check.tex      S6b. Chain check at 100 replicates (EUR+ASN 70:30):
#                              r^2 change at 4,000 iterations per arm, and the
#                              PRS-MIX vs single-panel contrasts at both lengths.
#   table_discovery.tex        S7.  PRS-MIX minus each single panel with reference
#                              panels and with discovery-built panels (panels
#                              from the GWAS's own individuals).
#   table_layouts.tex          S8a. r^2 by block layout (250 kb, 1 Mb, EUR cut
#                              points) and the PRS-MIX contrasts at each layout.
#   table_diagnostics.tex      S8b. Sampler-free diagnostics: share of the
#                              generating LD that falls between blocks, and the
#                              within-block error of each panel.
#   tables.html, tables.tsv, inputs/, manifest.json
#
# Intervals are 95% paired t intervals over replicates (as prespecified in
# LAYOUT_PILOT.md), not adjusted for multiplicity. LaTeX and HTML show the 250 kb
# layout where a table is per layout; tables.tsv has every layout.
#
# Usage (from the PRS-MIX root):
#   Rscript inst/scripts/null_sim/22_pilot_paper_tables.R <pilot_root> <chain_dir> <out_dir>
#
#   <pilot_root>  e.g. inst/scripts/results/layout_pilot, holding
#                 <setting>/runs/<layout>/metrics/nonnull_metrics.tsv and
#                 <setting>/diagnostics/{boundary,mismatch}.tsv
#   <chain_dir>   e.g. inst/scripts/results/chain_check/eur_asn_w0.7, holding
#                 metrics/nonnull_metrics.tsv

script_dir <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))
source(file.path(script_dir, "common.R"))
args <- ns_args("Rscript 22_pilot_paper_tables.R <pilot_root> <chain_dir> <out_dir>", 3)
pilot_root <- args[[1]]
chain_dir <- args[[2]]
out_dir <- args[[3]]
program <- "inst/scripts/null_sim/22_pilot_paper_tables.R"
dir.create(file.path(out_dir, "inputs"), recursive = TRUE, showWarnings = FALSE)
to_paper <- grepl("(^|/)paper/", normalizePath(out_dir))

sha <- function(path) sub("\\s.*$", "", system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE))
git <- function(...) system2("git", c(...), stdout = TRUE)

pilot <- list(eur_afr_w0.5 = list(A = "EUR", B = "AFR", w = 0.5, label = "EUR + AFR, 50:50"),
              eur_asn_w0.5 = list(A = "EUR", B = "ASN", w = 0.5, label = "EUR + ASN, 50:50"))
layouts <- c(fixed250kb = "250 kb", fixed1000kb = "1 Mb", eurcuts = "EUR cut points")

# ---- Freeze inputs -------------------------------------------------------
srcs <- c(unlist(lapply(names(pilot), function(s) c(
  file.path(s, "runs", names(layouts), "metrics", "nonnull_metrics.tsv"),
  file.path(s, "diagnostics", c("boundary.tsv", "mismatch.tsv"))))))
freeze <- function(src, dst_rel) {
  dst <- file.path(out_dir, "inputs", dst_rel)
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
  file.copy(src, dst, overwrite = TRUE)
  stopifnot(sha(src) == sha(dst))
  list(path = dst, source = src, sha256 = sha(dst),
       source_commit = c(git("log", "-1", "--format=%H", "--", src), NA_character_)[1],
       source_uncommitted = length(git("status", "--porcelain", "--", src)) > 0)
}
inputs <- c(lapply(srcs, function(f) freeze(file.path(pilot_root, f), file.path("layout_pilot", f))),
            list(freeze(file.path(chain_dir, "metrics", "nonnull_metrics.tsv"), "chain_check_eur_asn_w0.7_nonnull_metrics.tsv")))
names(inputs) <- c(file.path("layout_pilot", srcs), "chain_check")
inputs_dirty <- any(vapply(inputs, `[[`, TRUE, "source_uncommitted"))
prog_dirty <- length(git("status", "--porcelain", "--", program, "inst/scripts/null_sim/common.R")) > 0
if (to_paper && (inputs_dirty || prog_dirty)) {
  stop("Writing under paper/ needs committed inputs and program; commit them first.", call. = FALSE)
}
inp <- function(key) ns_read_tsv(inputs[[key]]$path)

# ---- Statistics ------------------------------------------------------------
cells <- list()
put <- function(table, group, row, target, column, x, layout = "fixed250kb") {
  x <- x[is.finite(x)]
  n <- length(x)
  m <- if (n) mean(x) else NA_real_
  h <- if (n > 1) stats::qt(0.975, n - 1) * stats::sd(x) / sqrt(n) else NA_real_
  cells[[length(cells) + 1]] <<- data.frame(table = table, layout = layout, group = group, row = row, target = target,
                                            column = column, n = n, mean = m, ci_lo = m - h, ci_hi = m + h)
}
# Per-replicate value of metric col for arm, at replicates reps.
val <- function(x, arm, col, reps) { y <- x[x$arm == arm, ]; y[[col]][match(reps, y$rep)] }
arm_rows <- function(s) c(single_A = paste(s$A, "panel"), single_B = paste(s$B, "panel"), mix = "PRS-MIX")
arm_name <- function(key, s, prefix = "") switch(key, single_A = paste0(prefix, "single_", s$A),
                                                 single_B = paste0(prefix, "single_", s$B),
                                                 mix = sprintf("%smix_w%.1f", prefix, s$w))

met <- list()
for (s_id in names(pilot)) for (l in names(layouts)) {
  x <- inp(file.path("layout_pilot", s_id, "runs", l, "metrics", "nonnull_metrics.tsv"))
  met[[s_id]][[l]] <- x[x$phi_mode == "auto", ]
}

# S6a: chain length in the pilot.
for (s_id in names(pilot)) for (l in names(layouts)) {
  s <- pilot[[s_id]]; x <- met[[s_id]][[l]]
  reps5 <- sort(unique(x$rep[x$arm == arm_name("mix", s, "long_")]))
  reps20 <- sort(unique(x$rep[x$arm == paste0("single_", s$A, "_seed2")]))
  has_x <- any(grepl("^xlong_", x$arm))
  for (k in names(arm_rows(s))) {
    a0 <- arm_name(k, s); a1 <- arm_name(k, s, "long_"); a2 <- arm_name(k, s, "xlong_")
    for (tg in c(s$A, s$B)) {
      cn <- paste0("r2_", tg)
      put("chain_pilot", s$label, arm_rows(s)[[k]], tg, "r2_default", val(x, a0, cn, reps5), l)
      put("chain_pilot", s$label, arm_rows(s)[[k]], tg, "d_long", val(x, a1, cn, reps5) - val(x, a0, cn, reps5), l)
      if (has_x) put("chain_pilot", s$label, arm_rows(s)[[k]], tg, "d_xlong", val(x, a2, cn, reps5) - val(x, a1, cn, reps5), l)
      # Spread from sampler randomness alone: SD over replicates 1-20 of the
      # first panel's seed rerun minus its main run (the same for every arm).
      put("chain_pilot", s$label, arm_rows(s)[[k]], tg, "seed_sd",
          stats::sd(val(x, paste0("single_", s$A, "_seed2"), cn, reps20) - val(x, paste0("single_", s$A), cn, reps20)), l)
    }
    put("chain_pilot", s$label, arm_rows(s)[[k]], "", "phi_default", val(x, a0, "phi", reps5), l)
    put("chain_pilot", s$label, arm_rows(s)[[k]], "", "phi_long", val(x, a1, "phi", reps5), l)
    if (has_x) put("chain_pilot", s$label, arm_rows(s)[[k]], "", "phi_xlong", val(x, a2, "phi", reps5), l)
  }
}

# S6b: chain check at 100 replicates.
cc <- inp("chain_check"); cc <- cc[cc$phi_mode == "auto", ]
s_cc <- list(A = "EUR", B = "ASN", w = 0.7, label = "EUR + ASN, 70:30")
reps_cc <- sort(Reduce(intersect, lapply(c("single_EUR", "single_ASN", "mix_w0.7", "long_single_EUR", "long_single_ASN", "long_mix_w0.7"),
                                         function(a) cc$rep[cc$arm == a])))
for (k in names(arm_rows(s_cc))) for (tg in c("EUR", "ASN")) {
  cn <- paste0("r2_", tg)
  a0 <- arm_name(k, s_cc); a1 <- arm_name(k, s_cc, "long_")
  put("chain_check", s_cc$label, arm_rows(s_cc)[[k]], tg, "r2_default", val(cc, a0, cn, reps_cc))
  put("chain_check", s_cc$label, arm_rows(s_cc)[[k]], tg, "r2_long", val(cc, a1, cn, reps_cc))
  put("chain_check", s_cc$label, arm_rows(s_cc)[[k]], tg, "d_long", val(cc, a1, cn, reps_cc) - val(cc, a0, cn, reps_cc))
}
for (k in c("single_A", "single_B")) for (tg in c("EUR", "ASN")) {
  cn <- paste0("r2_", tg)
  row <- sprintf("PRS-MIX minus %s", arm_rows(s_cc)[[k]])
  d0 <- val(cc, "mix_w0.7", cn, reps_cc) - val(cc, arm_name(k, s_cc), cn, reps_cc)
  d1 <- val(cc, "long_mix_w0.7", cn, reps_cc) - val(cc, arm_name(k, s_cc, "long_"), cn, reps_cc)
  put("chain_contrast", s_cc$label, row, tg, "contrast_default", d0)
  put("chain_contrast", s_cc$label, row, tg, "contrast_long", d1)
  put("chain_contrast", s_cc$label, row, tg, "contrast_change", d1 - d0)
}

# S7: discovery-built against reference panels; S8a: layouts.
for (s_id in names(pilot)) for (l in names(layouts)) {
  s <- pilot[[s_id]]; x <- met[[s_id]][[l]]
  reps <- sort(Reduce(intersect, lapply(c(arm_name("mix", s), arm_name("single_A", s), arm_name("single_B", s),
                                          arm_name("mix", s, "disc_"), arm_name("single_A", s, "disc_"), arm_name("single_B", s, "disc_")),
                                        function(a) x$rep[x$arm == a])))
  for (tg in c(s$A, s$B)) {
    cn <- paste0("r2_", tg)
    for (src in c("", "disc_")) for (k in c("single_A", "single_B")) {
      put("discovery", s$label, sprintf("PRS-MIX minus %s", arm_rows(s)[[k]]), tg, paste0(if (src == "") "reference" else "discovery"),
          val(x, arm_name("mix", s, src), cn, reps) - val(x, arm_name(k, s, src), cn, reps), l)
    }
    for (k in names(arm_rows(s))) {
      put("discovery_gain", s$label, arm_rows(s)[[k]], tg, "disc_minus_ref",
          val(x, arm_name(k, s, "disc_"), cn, reps) - val(x, arm_name(k, s), cn, reps), l)
    }
  }
}
for (s_id in names(pilot)) {
  s <- pilot[[s_id]]
  x0 <- met[[s_id]][["fixed250kb"]]
  for (tg in c(s$A, s$B)) {
    cn <- paste0("r2_", tg)
    for (k in names(arm_rows(s))) {
      reps <- sort(Reduce(intersect, lapply(names(layouts), function(l) met[[s_id]][[l]]$rep[met[[s_id]][[l]]$arm == arm_name(k, s)])))
      r250 <- val(x0, arm_name(k, s), cn, reps)
      put("layouts", s$label, arm_rows(s)[[k]], tg, "r2_250kb", r250)
      put("layouts", s$label, arm_rows(s)[[k]], tg, "d_1mb", val(met[[s_id]][["fixed1000kb"]], arm_name(k, s), cn, reps) - r250)
      put("layouts", s$label, arm_rows(s)[[k]], tg, "d_eurcuts", val(met[[s_id]][["eurcuts"]], arm_name(k, s), cn, reps) - r250)
    }
  }
}

# S8b: sampler-free diagnostics (one value per setting; no intervals).
diag <- list()
for (s_id in names(pilot)) {
  s <- pilot[[s_id]]
  bd <- inp(file.path("layout_pilot", s_id, "diagnostics", "boundary.tsv"))
  mi <- inp(file.path("layout_pilot", s_id, "diagnostics", "mismatch.tsv"))
  for (i in seq_len(nrow(bd))) diag[[length(diag) + 1]] <- data.frame(
    table = "diagnostics", layout = bd$layout[i], group = s$label, row = sprintf("%s LD (discovery n = %d)", bd$pop[i], bd$n_discovery[i]),
    target = "", column = "cross_share", n = 1, mean = 100 * bd$cross_share_adj[i], ci_lo = NA, ci_hi = NA)
  for (i in seq_len(nrow(mi))) diag[[length(diag) + 1]] <- data.frame(
    table = "diagnostics", layout = mi$layout[i], group = s$label, row = sprintf("%s LD", mi$ld[i]),
    target = "", column = paste0("rel_err_", mi$source[i]), n = 1, mean = 100 * mi$rel_err[i], ci_lo = NA, ci_hi = NA)
}
cells <- rbind(do.call(rbind, cells), do.call(rbind, diag))
ns_write_tsv(cells, file.path(out_dir, "tables.tsv"))

# ---- Rendering -------------------------------------------------------------
tx <- function(tex, t, h) if (tex) t else h
fmt <- function(x, d, signed, tex) {
  if (!is.finite(x)) return(tx(tex, "--", "&ndash;"))
  s <- formatC(abs(x), format = "f", digits = d)
  if (as.numeric(s) == 0) return(s)
  paste0(if (x < 0) tx(tex, "$-$", "&minus;") else if (signed) "+" else "", s)
}
cell <- function(table, group, row, target, column, d, signed = TRUE, ci = TRUE, tex = TRUE, layout = "fixed250kb") {
  r <- cells[cells$table == table & cells$layout == layout & cells$group == group & cells$row == row &
               cells$target == target & cells$column == column, ]
  if (nrow(r) == 0) return(tx(tex, "--", "&ndash;"))
  if (!ci || !is.finite(r$ci_lo)) return(fmt(r$mean, d, signed, tex))
  sprintf("%s [%s, %s]", fmt(r$mean, d, signed, tex), fmt(r$ci_lo, d, signed, tex), fmt(r$ci_hi, d, signed, tex))
}
# A table: header labels plus groups of rows; each row is a label and its cell strings.
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
phi <- function(tex) tx(tex, "$\\phi$", "<i>&phi;</i>")
build <- function(tex) {
  out <- list()
  # S6a
  heads <- c("Arm, scored in", paste(r2(tex), "at defaults"), "4,000 minus defaults", "8,000 minus 4,000", "Seed-rerun SD",
             paste(phi(tex), "(defaults, 4,000, 8,000)"))
  groups <- list()
  for (s_id in names(pilot)) {
    s <- pilot[[s_id]]
    groups[[s$label]] <- unlist(lapply(names(arm_rows(s)), function(k) lapply(c(s$A, s$B), function(tg) {
      a <- arm_rows(s)[[k]]
      ph <- vapply(c("phi_default", "phi_long", "phi_xlong"), function(cn) {
        r <- cells[cells$table == "chain_pilot" & cells$layout == "fixed250kb" & cells$group == s$label & cells$row == a & cells$column == cn, ]
        if (nrow(r)) formatC(signif(r$mean, 2), format = "fg", digits = 2) else tx(tex, "--", "&ndash;")
      }, "")
      c(sprintf("%s, %s", a, tg),
        cell("chain_pilot", s$label, a, tg, "r2_default", 3, FALSE, FALSE, tex),
        cell("chain_pilot", s$label, a, tg, "d_long", 3, tex = tex),
        cell("chain_pilot", s$label, a, tg, "d_xlong", 3, tex = tex),
        cell("chain_pilot", s$label, a, tg, "seed_sd", 3, FALSE, FALSE, tex),
        if (tg == s$A) paste(ph, collapse = ", ") else "")
    })), recursive = FALSE)
  }
  out$chain_pilot <- render(heads, groups, tex)
  # S6b
  heads <- c("", paste(r2(tex), "at defaults"), paste(r2(tex), "at 4,000"), tx(tex, "Change [95\\% CI]", "Change [95% CI]"))
  rows <- unlist(lapply(names(arm_rows(s_cc)), function(k) lapply(c("EUR", "ASN"), function(tg) {
    a <- arm_rows(s_cc)[[k]]
    c(sprintf("%s, scored in %s", a, tg), cell("chain_check", s_cc$label, a, tg, "r2_default", 3, FALSE, FALSE, tex),
      cell("chain_check", s_cc$label, a, tg, "r2_long", 3, FALSE, FALSE, tex), cell("chain_check", s_cc$label, a, tg, "d_long", 4, tex = tex))
  })), recursive = FALSE)
  crow <- unlist(lapply(c("single_A", "single_B"), function(k) lapply(c("EUR", "ASN"), function(tg) {
    a <- sprintf("PRS-MIX minus %s", arm_rows(s_cc)[[k]])
    c(sprintf("%s, scored in %s", a, tg), cell("chain_contrast", s_cc$label, a, tg, "contrast_default", 4, tex = tex),
      cell("chain_contrast", s_cc$label, a, tg, "contrast_long", 4, tex = tex), cell("chain_contrast", s_cc$label, a, tg, "contrast_change", 4, tex = tex))
  })), recursive = FALSE)
  out$chain_check <- render(heads, setNames(list(rows, crow), c(paste(r2(tex), "of each arm"), "Paired contrasts")), tex)
  # S7
  heads <- c("Contrast, scored in", "Reference panels", "Discovery-built panels")
  groups <- list()
  for (s_id in names(pilot)) {
    s <- pilot[[s_id]]
    groups[[s$label]] <- c(
      unlist(lapply(c("single_A", "single_B"), function(k) lapply(c(s$A, s$B), function(tg) {
        a <- sprintf("PRS-MIX minus %s", arm_rows(s)[[k]])
        c(sprintf("%s, %s", a, tg), cell("discovery", s$label, a, tg, "reference", 3, tex = tex), cell("discovery", s$label, a, tg, "discovery", 3, tex = tex))
      })), recursive = FALSE),
      unlist(lapply(names(arm_rows(s)), function(k) lapply(c(s$A, s$B), function(tg) {
        a <- arm_rows(s)[[k]]
        c(sprintf("%s: discovery minus reference, %s", a, tg), "", cell("discovery_gain", s$label, a, tg, "disc_minus_ref", 3, tex = tex))
      })), recursive = FALSE))
  }
  out$discovery <- render(heads, groups, tex)
  # S8a
  heads <- c("Arm, scored in", paste(r2(tex), "at 250 kb"), "1 Mb minus 250 kb", "EUR cut points minus 250 kb")
  groups <- list()
  for (s_id in names(pilot)) {
    s <- pilot[[s_id]]
    groups[[s$label]] <- unlist(lapply(names(arm_rows(s)), function(k) lapply(c(s$A, s$B), function(tg) {
      a <- arm_rows(s)[[k]]
      c(sprintf("%s, %s", a, tg), cell("layouts", s$label, a, tg, "r2_250kb", 3, FALSE, FALSE, tex),
        cell("layouts", s$label, a, tg, "d_1mb", 3, tex = tex), cell("layouts", s$label, a, tg, "d_eurcuts", 3, tex = tex))
    })), recursive = FALSE)
  }
  out$layouts <- render(heads, groups, tex)
  # S8b
  pct <- tx(tex, "\\%", "%")
  heads <- c("LD", paste0("Between-block share at 250 kb, 1 Mb, EUR cuts (", pct, ")"),
             paste0("Within-block error, reference panel (", pct, ")"), paste0("Within-block error, discovery panel (", pct, ")"))
  groups <- list()
  for (s_id in names(pilot)) {
    s <- pilot[[s_id]]
    groups[[s$label]] <- lapply(c(s$A, s$B), function(p) {
      d <- cells[cells$table == "diagnostics" & cells$group == s$label, ]
      g <- function(col, l) { r <- d[d$column == col & d$layout == l & startsWith(d$row, paste(p, "LD")), ]; if (nrow(r)) formatC(r$mean[1], format = "f", digits = 1) else tx(tex, "--", "&ndash;") }
      three <- function(col) paste(vapply(names(layouts), function(l) g(col, l), ""), collapse = ", ")
      c(sprintf("%s LD", p), three("cross_share"), three("rel_err_reference"), three("rel_err_discovery"))
    })
  }
  out$diagnostics <- render(heads, groups, tex)
  out
}
tex_tabs <- build(TRUE)
for (nm in names(tex_tabs)) writeLines(tex_tabs[[nm]], file.path(out_dir, sprintf("table_%s.tex", nm)))
html_tabs <- build(FALSE)

titles <- c(chain_pilot = "Table S6a. Sampler chain length in the layout pilot",
            chain_check = "Table S6b. Chain-length check at 100 replicates",
            discovery = "Table S7. Reference panels against panels built from the GWAS's own individuals",
            layouts = "Table S8a. Accuracy by LD block layout",
            diagnostics = "Table S8b. How much LD each layout and panel misses")
captions <- c(
  chain_pilot = paste("Layout pilot, chromosome 22, 250 kb blocks, replicates 1&ndash;5. PRS-CS at its defaults (1,000 iterations, 500 burn-in)",
                      "against 4,000 iterations (2,000 burn-in) and, for EUR + ASN only, 8,000 iterations (4,000 burn-in); thinning 5 throughout.",
                      "Differences are paired within replicate, with 95% <i>t</i> intervals. Seed-rerun SD: the SD over replicates 1&ndash;20 of the",
                      "first panel's rerun with a new seed minus its main run, the spread from sampler randomness alone.",
                      "<i>&phi;</i>: the learned global shrinkage, mean over replicates, to two significant figures."),
  chain_check = paste("EUR + ASN 70:30, chromosome 22, 250 kb blocks, 100 replicates: PRS-CS at its defaults against 4,000 iterations",
                      "(2,000 burn-in, thinning 5). Paired differences with 95% <i>t</i> intervals. The change column of the contrasts",
                      "is the 4,000-iteration contrast minus the default contrast."),
  discovery = paste("Layout pilot, chromosome 22, 250 kb blocks, replicates 1&ndash;20, GWAS mix 50:50. <i>r</i><sup>2</sup> with the true",
                    "genetic value, paired differences with 95% <i>t</i> intervals. Reference panels are built from the reference individuals,",
                    "as in the main analysis; discovery-built panels from the individuals the GWAS was simulated on, so they differ from the",
                    "LD the GWAS carries only by the panel estimator. Discovery-built panels are a diagnostic, not a usable method."),
  layouts = paste("Layout pilot, chromosome 22, replicates 1&ndash;20, GWAS mix 50:50, reference panels. <i>r</i><sup>2</sup> with the true",
                  "genetic value at fixed 250 kb blocks, and the paired change at fixed 1 Mb blocks and at PRS-CS's EUR LD-block cut points,",
                  "with 95% <i>t</i> intervals. Every layout uses the same SNPs."),
  diagnostics = paste("Sampler-free, chromosome 22. Between-block share: of the generating LD's squared correlations between SNPs",
                      "within 3 Mb, the percentage between blocks, which block-diagonal LD discards (with 1/<i>n</i> sampling noise",
                      "removed). Within-block error: the panel's squared off-diagonal error within blocks against the generating LD, as a",
                      "percentage of the generating LD's own off-diagonal sum of squares. Values for 250 kb, 1 Mb and EUR cut points."))
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
          sprintf("<h1>%s</h1>", "Layout pilot and chain-length checks: supplementary tables"),
          unlist(lapply(names(html_tabs), function(nm) c("<section>", sprintf("<h2>%s</h2>", titles[[nm]]),
                                                         sprintf("<p class=\"cap\">%s</p>", captions[[nm]]),
                                                         "<div class=\"wrap\">", html_tabs[[nm]], "</div>", "</section>"))),
          sprintf("<p class=\"src\">Generated by %s at PRS-MIX %s%s.</p>", program, substr(head_commit, 1, 7),
                  if (inputs_dirty || prog_dirty) " with uncommitted program or inputs" else ""),
          "</main>", "</body>", "</html>")
writeLines(html, file.path(out_dir, "tables.html"))

# ---- Manifest ------------------------------------------------------------
out_files <- c(sprintf("table_%s.tex", names(tex_tabs)), "tables.html", "tables.tsv")
manifest <- list(
  study = basename(normalizePath(out_dir)), program = program,
  command = paste("Rscript", program, paste(shQuote(args), collapse = " ")),
  asset_commit = head_commit, program_uncommitted = prog_dirty, inputs_uncommitted = inputs_dirty,
  experiment_commit = unique(stats::na.omit(vapply(inputs, `[[`, "", "source_commit"))),
  r_version = R.version.string,
  inputs = unname(lapply(inputs, function(x) x[c("path", "source", "source_commit", "sha256")])),
  outputs = lapply(stats::setNames(out_files, out_files), function(f) list(path = file.path(out_dir, f), sha256 = sha(file.path(out_dir, f)))),
  settings = list(
    design = "inst/scripts/null_sim/LAYOUT_PILOT.md (prespecified 2026-11-05; deviation 2026-11-09: 8,000-iteration chains for EUR+ASN)",
    pilot_settings = unname(vapply(pilot, `[[`, "", "label")), pilot_replicates = "1-20; chain arms 1-5",
    chain_check = "EUR+ASN 70:30, 250 kb, replicates with every default and 4,000-iteration arm",
    chain_check_replicates = length(reps_cc),
    intervals = "95% paired t intervals over replicates, unadjusted", phi_mode = "auto",
    sampler = "defaults 1,000/500/5; long 4,000/2,000/5; xlong 8,000/4,000/5 (iterations/burn-in/thinning)"),
  notes = "Pilot and chain-check runs were made on RedHawk (cluster/jobs/layout_pilot, not version-controlled).")
jsonlite::write_json(manifest, file.path(out_dir, "manifest.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA)
message(sprintf("Wrote %s (%d cells) at PRS-MIX %s%s", out_dir, nrow(cells), substr(head_commit, 1, 7),
                if (inputs_dirty || prog_dirty) " with uncommitted program or inputs" else ""))
