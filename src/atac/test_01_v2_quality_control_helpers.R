## Extract the pure-function definitions from the script and exercise them on
## synthetic data with known structure. No ArchR / no data needed.
set.seed(1)

# Locate the script: next to this test file, or in the working directory.
script_candidates <- c(
  "01_v2_quality_control.R",
  file.path(dirname(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "01_v2_quality_control.R"),
  "src/atac/01_v2_quality_control.R"
)
script_path <- script_candidates[which(file.exists(script_candidates))[1L]]
if (is.na(script_path)) stop("Cannot find 01_v2_quality_control.R next to this test file.")
cat("testing against:", script_path, "\n")
exprs <- parse(script_path)
want <- c("is_group_constant", "collapse_to_group", "assoc_test",
          "first_unit", "safe_max", "safe_at_max",
          "gap_threshold", "half_min_nonzero", "bin_to_numeric",
          "strip_timepoint", "neglog10")
env <- new.env()
got <- character(0)
for (e in exprs) {
  if (is.call(e) && length(e) >= 3L &&
      (identical(e[[1]], as.name("<-")) || identical(e[[1]], as.name("=")))) {
    nm <- if (is.name(e[[2]])) as.character(e[[2]]) else ""
    if (length(nm) == 1L && nm %in% want) { eval(e, envir = env); got <- c(got, nm) }
  }
}
cat("extracted:", paste(sort(got), collapse = ", "), "\n")
stopifnot(setequal(got, want))
for (n in want) assign(n, get(n, envir = env))

ok <- function(label, cond) cat(sprintf("%-58s %s\n", label, if (isTRUE(cond)) "PASS" else "** FAIL **"))

## ---- strip_timepoint / bin_to_numeric --------------------------------------
ok("strip_timepoint collapses 555d0/555d2 to one donor",
   strip_timepoint("28205-0555d0") == strip_timepoint("28205-0555d2"))
ok("strip_timepoint leaves plain ids alone",
   strip_timepoint("HIP043") == "HIP043" && strip_timepoint("9313454") == "9313454")
ok("bin_to_numeric midpoints and open-ended bins",
   isTRUE(all.equal(bin_to_numeric(c("30-39", "80+", "45", NA)), c(34.5, 85, 45, NA))))

## ---- gap_threshold round-trips on the raw scale ----------------------------
x <- c(runif(20, 1e-6, 5e-6), runif(20, 1e-3, 5e-3))   # two clear clusters
eps <- half_min_nonzero(x)
thr <- gap_threshold(x, eps)
ok("gap_threshold separates the two clusters on the raw scale",
   thr > 5e-6 && thr < 1e-3)

## ---- is_group_constant ------------------------------------------------------
g  <- rep(paste0("D", 1:10), each = 3)
xc <- rep(letters[1:10], each = 3)          # donor-constant
xv <- rep(1:3, times = 10)                  # varies within donor
ok("is_group_constant TRUE for a donor-constant covariate", is_group_constant(xc, g))
ok("is_group_constant FALSE for a within-donor covariate", !is_group_constant(xv, g))
ok("is_group_constant ignores NAs", is_group_constant(replace(xc, c(2, 5), NA), g))

## ---- collapse_to_group ------------------------------------------------------
y  <- as.numeric(seq_along(g))
cl <- collapse_to_group(y, factor(xc), factor(g))
ok("collapse_to_group returns one row per donor", cl$n == 10 && length(cl$y) == 10)
ok("collapse_to_group averages y within donor",
   isTRUE(all.equal(cl$y[1], mean(y[g == "D1"]))))
ok("collapse_to_group keeps x a factor", is.factor(cl$x) && nlevels(cl$x) == 10)

## ---- assoc_test: basic dispatch --------------------------------------------
donor  <- rep(paste0("D", 1:8), each = 3)              # 24 samples, 8 donors
cohort <- factor(rep(c("ctl", "exp"), each = 12))      # a DONOR property
yy     <- rnorm(24)
res_donor <- assoc_test(yy, cohort, donor)
ok("assoc_test picks the donor unit for a donor-constant covariate",
   identical(res_donor$unit, "donor"))
ok("assoc_test reports n = number of donors, not samples",
   res_donor$n == 8 && res_donor$n_group == 8)

## ---- THE CENTRAL CLAIM: calibration under the null --------------------------
# 24 samples from 8 donors (3 each). Donor identity drives the coordinate;
# cohort is assigned at the donor level and has NO effect. A correct test
# rejects at ~5%. A sample-level OLS that treats the 3 samples of a donor as
# independent should reject far more often. This is the whole reason the script
# groups on donor, so it is worth checking rather than asserting.
nsim <- 1000
naive_rej <- donor_rej <- logical(nsim)
for (s in seq_len(nsim)) {
  de <- setNames(rnorm(8, sd = 1), paste0("D", 1:8))   # donor random intercepts
  y  <- de[donor] + rnorm(24, sd = 0.1)                # high intra-donor correlation
  naive_rej[s] <- anova(lm(y ~ 1), lm(y ~ cohort))[["Pr(>F)"]][2] < 0.05
  donor_rej[s] <- isTRUE(assoc_test(y, cohort, donor)$p < 0.05)
}
cat(sprintf("   null rejection rate at alpha = 0.05 over %d simulations:\n", nsim))
cat(sprintf("     naive sample-level OLS : %.3f   (should be ~0.05)\n", mean(naive_rej)))
cat(sprintf("     donor-level assoc_test : %.3f   (should be ~0.05)\n", mean(donor_rej)))
ok("naive sample-level OLS is anti-conservative (rejects > 15% under the null)",
   mean(naive_rej) > 0.15)
ok("donor-level assoc_test is calibrated (rejects 2-9% under the null)",
   mean(donor_rej) > 0.02 && mean(donor_rej) < 0.09)

## ---- assoc_test: a real donor-level effect is still detected ---------------
y2 <- ifelse(cohort == "exp", 3, 0) + rnorm(24, sd = 0.3)
ok("a genuine cohort effect survives the donor collapse",
   assoc_test(y2, cohort, donor)$p < 0.05)

## ---- assoc_test: adjusted R^2 penalises degrees of freedom -----------------
# a 4-level factor on 8 donors explains a lot of pure noise by construction
noise_f <- factor(rep(c("a", "b", "c", "d"), each = 6))
rs <- assoc_test(rnorm(24), noise_f, donor)
cat(sprintf("   noise factor: raw R2 = %.2f, adjusted R2 = %.2f\n", rs$r2, rs$adj_r2))
ok("raw R2 is inflated for a multi-level factor on few donors", rs$r2 > 0.15)
ok("adjusted R2 discounts it", rs$adj_r2 < rs$r2)

# a factor with one level per donor leaves zero residual df: the test is not
# defined and must come back NA rather than a spurious R2 of 1
ok("saturated factor -> NA (no residual df)",
   is.na(assoc_test(rnorm(24), factor(donor), donor)$p))

## ---- assoc_test: cohort-adjusted partial test -------------------------------
# age is perfectly aliased with cohort; beyond cohort it explains nothing
age <- ifelse(cohort == "exp", 60, 30)
rp  <- assoc_test(y2, age, donor, cohort)
ok("partial test of a covariate aliased with cohort returns NA/non-significant",
   is.na(rp$p) || rp$p > 0.05)

# a covariate that varies across donors WITHIN cohort is still detectable
age2 <- as.numeric(setNames(c(20, 30, 40, 50, 22, 32, 42, 52), paste0("D", 1:8))[donor])
y3   <- 0.1 * age2 + rnorm(24, sd = 0.2)
ok("partial test detects a real within-cohort donor covariate",
   assoc_test(y3, age2, donor, cohort)$p < 0.05)

## ---- assoc_test: degenerate inputs return NA rather than erroring ----------
ok("constant covariate -> NA", is.na(assoc_test(yy, rep(1, 24), donor)$p))
ok("single-level factor -> NA", is.na(assoc_test(yy, factor(rep("a", 24)), donor)$p))
ok("all-NA covariate -> NA", is.na(assoc_test(yy, rep(NA_real_, 24), donor)$p))
ok("too few observations -> NA", is.na(assoc_test(yy[1:3], cohort[1:3], donor[1:3])$p))

## ---- dispatch to the LMM branch --------------------------------------------
within_donor <- rep(1:3, 8)
ok("a within-donor covariate is NOT sent to the donor collapse",
   !is_group_constant(within_donor, donor))

## ---- summary helpers tolerate all-NA groups --------------------------------
ok("first_unit on all-NA", is.na(first_unit(c(NA_character_, NA_character_))))
ok("first_unit picks the first non-NA", first_unit(c(NA, "donor", "x")) == "donor")
ok("safe_max on all-NA", is.na(safe_max(c(NA_real_, NA_real_))))
ok("safe_at_max on all-NA", is.na(safe_at_max(c("a", "b"), c(NA_real_, NA_real_))))
ok("safe_at_max picks the right label",
   safe_at_max(c("Dim1", "Dim2", "Dim3"), c(0.1, 0.9, NA)) == "Dim2")

## ---- neglog10 caps instead of producing Inf --------------------------------
p_cap <- 50   # defined immediately above neglog10() in the script
ok("neglog10 caps an underflowed p at the cap rather than Inf",
   is.finite(neglog10(0)) && neglog10(0) == 50)
