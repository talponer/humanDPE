# Generate resampled data sets

Rscript --vanilla ./scripts/hg19_bootstrap_split.R

# Phase of 1 of extended algorithm - collections of motifs
# (executed in parallel on multiple cores)

Rscript --vanilla ./scripts/hg19_bootstrap_phase1.R  0 >/dev/null 2>&1 &
Rscript --vanilla ./scripts/hg19_bootstrap_phase1.R  1 >/dev/null 2>&1 &
Rscript --vanilla ./scripts/hg19_bootstrap_phase1.R  2 >/dev/null 2>&1 &
Rscript --vanilla ./scripts/hg19_bootstrap_phase1.R  3 >/dev/null 2>&1 &
Rscript --vanilla ./scripts/hg19_bootstrap_phase1.R  4 >/dev/null 2>&1 &
Rscript --vanilla ./scripts/hg19_bootstrap_phase1.R  5 >/dev/null 2>&1 &
Rscript --vanilla ./scripts/hg19_bootstrap_phase1.R  6 >/dev/null 2>&1 &
Rscript --vanilla ./scripts/hg19_bootstrap_phase1.R  7 >/dev/null 2>&1 &
Rscript --vanilla ./scripts/hg19_bootstrap_phase1.R  8 >/dev/null 2>&1 &
Rscript --vanilla ./scripts/hg19_bootstrap_phase1.R  9 >/dev/null 2>&1 &
Rscript --vanilla ./scripts/hg19_bootstrap_phase1.R 10 >/dev/null 2>&1 &

# Phase 2 of extended algoritm

Rscript --vanilla ./scripts/hg19_bootstrap_phase2.R

# Correlation plot

Rscript --vanilla figS2.R 

