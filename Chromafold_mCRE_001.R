#!/usr/bin/env Rscript
## Chromafold_mCRE_001.R
## Elliott Ferris
## Feb 23 2026
##
## Classify Chromafold-predicted enhancer-promoter contacts as single-gene CREs
## (sCREs) or multi-gene CREs (mCREs).
##
## Input:  A gzipped BEDPE file produced by intersecting Chromafold predictions
##         (bedtools pairtopair -type either) with a promoter BED file (EPD1kb_v6).
##         The file name must follow the pattern:
##           chromafold_<CellType>_<FeedingState>_top10_EPD1kb.bedpe.gz
##
## Output: A single tab-delimited line printed to stdout:
##           CellType  FeedingState  N_contacts  N_sCRE  N_mCRE
##
## Usage:
##   Rscript Chromafold_mCRE_001.R chromafold_Microglia_Fed_top10_EPD1kb.bedpe.gz
##
## To accumulate results across conditions:
##   for f in chromafold_*_EPD1kb.bedpe.gz; do
##     Rscript Chromafold_mCRE_001.R "$f"
##   done >> ChromafoldCRE_001.txt

library(limma)  # for strsplit2()

## --- Parse command-line arguments -------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript Chromafold_mCRE_001.R <bedpe.gz>")
}

cell_type     <- strsplit2(args[1], split = "_")[2]
feeding_state <- strsplit2(args[1], split = "_")[3]

## --- Read BEDPE file --------------------------------------------------------
contacts <- read.csv(gzfile(args[1]), sep = "\t", header = FALSE,
                     stringsAsFactors = FALSE)

if (nrow(contacts) == 0) {
  cat(paste(cell_type, feeding_state, 0, 0, 0, sep = "\t"), "\n")
  quit(save = "no")
}

colnames(contacts) <- c(
  "chr1", "start1", "end1",
  "chr2", "start2", "end2",
  "ChromafoldScore",
  "chr_promo", "start_promo", "end_promo",
  "promoter_name", "score", "strand",
  "EPD_start", "EPD_end"
)

contacts$GeneName      <- strsplit2(contacts$promoter_name, split = "_")[, 1]
contacts$ContactString <- paste(contacts$chr1, contacts$start1,
                                contacts$start2, sep = "_")

## --- Determine which bin contains the annotated promoter (PromoSide) --------
left_hit <- contacts$chr_promo == contacts$chr1 &
  ((contacts$start_promo <= contacts$end1 & contacts$start_promo >= contacts$start1) |
   (contacts$end_promo   <= contacts$end1 & contacts$end_promo   >= contacts$start1))

right_hit <- contacts$chr_promo == contacts$chr2 &
  ((contacts$start_promo <= contacts$end2 & contacts$start_promo >= contacts$start2) |
   (contacts$end_promo   <  contacts$end2 & contacts$end_promo   >= contacts$start2))

contacts$PromoSide <- ifelse(
  left_hit & right_hit, "Both",
  ifelse(left_hit, "Left",
         ifelse(right_hit, "Right", NA_character_))
)

## --- Check if BOTH bins overlap a promoter (PromoPromo) ---------------------
promos <- unique(contacts[, c("chr_promo", "start_promo", "end_promo")])

bin_has_promo <- function(chrom, start, end) {
  mapply(function(ch, s, e) {
    p <- promos[promos$chr_promo == ch, ]
    any((p$start_promo <= e & p$start_promo >= s) |
        (p$end_promo   <= e & p$end_promo   >= s))
  }, chrom, start, end, USE.NAMES = FALSE)
}

contacts$PromoPromo <-
  bin_has_promo(contacts$chr1, contacts$start1, contacts$end1) &
  bin_has_promo(contacts$chr2, contacts$start2, contacts$end2)

contacts$ContactCategory <- contacts$PromoSide
contacts$ContactCategory[contacts$PromoPromo] <- "Both"

## --- Keep only enhancer-to-promoter contacts --------------------------------
enh_promo <- unique(contacts[!contacts$PromoPromo, ])

## --- Identify the cis-regulatory bin (the non-promoter side) ----------------
enh_promo$CisStart <- NA
enh_promo$CisEnd   <- NA

left_idx  <- enh_promo$ContactCategory == "Left"
right_idx <- enh_promo$ContactCategory == "Right"

enh_promo$CisStart[left_idx]  <- enh_promo$start2[left_idx]
enh_promo$CisEnd[left_idx]    <- enh_promo$end2[left_idx]
enh_promo$CisStart[right_idx] <- enh_promo$start1[right_idx]
enh_promo$CisEnd[right_idx]   <- enh_promo$end1[right_idx]

enh_promo$CisBin <- paste0(enh_promo$chr1, "_", enh_promo$CisStart)

## --- Summarise each cis-regulatory bin --------------------------------------
cis_summary <- do.call(rbind, lapply(
  split(enh_promo, enh_promo$CisBin),
  function(df) {
    data.frame(
      CisBin      = df$CisBin[1],
      N_Promoters = length(unique(df$promoter_name)),
      N_Genes     = length(unique(df$GeneName)),
      GeneNames   = paste(unique(df$GeneName), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
))

cis_summary$CREclass <- ifelse(cis_summary$N_Genes > 1, "mCRE", "sCRE")

n_sCRE <- sum(cis_summary$CREclass == "sCRE")
n_mCRE <- sum(cis_summary$CREclass == "mCRE")

## --- Output -----------------------------------------------------------------
cat(paste(cell_type, feeding_state, nrow(enh_promo), n_sCRE, n_mCRE, sep = "\t"),
    "\n")
