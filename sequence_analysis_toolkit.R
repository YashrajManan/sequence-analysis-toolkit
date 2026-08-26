# Sequence Analysis Toolkit - R implementation
#
# Companion to the Python version (sequence_analysis_toolkit.ipynb): the same DNA
# sequence analysis, written from first principles in R and cross-checked against
# the seqinr package. Worked example: the human insulin (INS) mRNA.
#
# Data: NCBI Entrez API (RefSeq accession NM_000207), fetched via rentrez.
# Run:  open in RStudio -> Session -> Set Working Directory -> To Source File
#       Location, then Source. First run installs rentrez + seqinr (~1 min).

## ---- Setup ----
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
if (!requireNamespace("rentrez", quietly = TRUE)) install.packages("rentrez")
if (!requireNamespace("seqinr",  quietly = TRUE)) install.packages("seqinr")
library(rentrez)
library(seqinr)
dir.create("data_R",    showWarnings = FALSE)
dir.create("results_R", showWarnings = FALSE)

## ---- Fetch from the NCBI Entrez API ----
accession <- "NM_000207"   # human insulin (INS) mRNA
fasta_txt <- entrez_fetch(db = "nuccore", id = accession, rettype = "fasta", retmode = "text")
writeLines(fasta_txt, "data_R/sequence.fasta")

## ---- Parse FASTA (custom parser) ----
parse_fasta <- function(path) {
  lines <- readLines(path); seqs <- list(); current_id <- NULL
  for (line in lines) {
    line <- trimws(line)
    if (nchar(line) == 0) next
    if (startsWith(line, ">")) {
      current_id <- strsplit(substring(line, 2), " ")[[1]][1]  # id = header, first word
      seqs[[current_id]] <- ""
    } else {
      seqs[[current_id]] <- paste0(seqs[[current_id]], line)
    }
  }
  seqs
}
records <- parse_fasta("data_R/sequence.fasta")
seq <- records[[1]]

## ---- GC content ----
gc_content <- function(s) {
  chars <- strsplit(toupper(s), "")[[1]]
  if (length(chars) == 0) return(0)
  (sum(chars == "G") + sum(chars == "C")) / length(chars)
}

## ---- GC skew (non-overlapping windows) ----
gc_skew <- function(s, window = 50) {
  s <- toupper(s); n <- nchar(s); skews <- c()
  for (i in seq(1, n - window + 1, by = window)) {
    v <- strsplit(substr(s, i, i + window - 1), "")[[1]]
    g <- sum(v == "G"); c <- sum(v == "C")
    skews <- c(skews, if (g + c == 0) 0 else (g - c) / (g + c))
  }
  skews
}

## ---- k-mer counts ----
kmer_counts <- function(s, k = 3) {
  s <- toupper(s)
  table(sapply(1:(nchar(s) - k + 1), function(i) substr(s, i, i + k - 1)))
}

## ---- Reverse complement ----
reverse_complement <- function(s) {
  comp  <- c(A = "T", T = "A", G = "C", C = "G", N = "N")   # named vector = lookup table
  chars <- strsplit(toupper(s), "")[[1]]
  paste(rev(comp[chars]), collapse = "")
}

## ---- Open reading frames (three forward frames) ----
stops <- c("TAA", "TAG", "TGA")
find_orfs <- function(s, min_len = 60) {
  s <- toupper(s); n <- nchar(s); orfs <- c()
  for (frame in 0:2) {
    i <- 1 + frame
    while (i <= n - 2) {
      if (substr(s, i, i + 2) == "ATG") {
        j <- i
        while (j <= n - 2) {
          if (substr(s, j, j + 2) %in% stops) {
            orf <- substr(s, i, j + 2)
            if (nchar(orf) >= min_len) orfs <- c(orfs, orf)
            i <- j; break
          }
          j <- j + 3
        }
      }
      i <- i + 3
    }
  }
  orfs
}

## ---- Translation (codon table -> protein) ----
CODON <- c(
 TTT="F",TTC="F",TTA="L",TTG="L",CTT="L",CTC="L",CTA="L",CTG="L",
 ATT="I",ATC="I",ATA="I",ATG="M",GTT="V",GTC="V",GTA="V",GTG="V",
 TCT="S",TCC="S",TCA="S",TCG="S",CCT="P",CCC="P",CCA="P",CCG="P",
 ACT="T",ACC="T",ACA="T",ACG="T",GCT="A",GCC="A",GCA="A",GCG="A",
 TAT="Y",TAC="Y",TAA="*",TAG="*",CAT="H",CAC="H",CAA="Q",CAG="Q",
 AAT="N",AAC="N",AAA="K",AAG="K",GAT="D",GAC="D",GAA="E",GAG="E",
 TGT="C",TGC="C",TGA="*",TGG="W",CGT="R",CGC="R",CGA="R",CGG="R",
 AGT="S",AGC="S",AGA="R",AGG="R",GGT="G",GGC="G",GGA="G",GGG="G")
translate_dna <- function(dna) {
  dna <- toupper(dna); protein <- ""
  for (i in seq(1, nchar(dna) - 2, by = 3)) {
    aa <- CODON[substr(dna, i, i + 2)]           # single [] -> NA if codon not found
    protein <- paste0(protein, if (is.na(aa)) "X" else aa)
  }
  protein
}

## ---- Run + cross-check against seqinr ----
record_id <- names(records)[1]
cat(record_id, "| length", nchar(seq), "\n")
cat("GC content:", round(gc_content(seq), 4), "\n")
stopifnot(round(gc_content(seq), 3) == round(GC(strsplit(tolower(seq), "")[[1]]), 3))

skews  <- gc_skew(seq, 50)
counts <- kmer_counts(seq, 3)
stopifnot(reverse_complement(reverse_complement(seq)) == seq)   # involution holds

orfs    <- find_orfs(seq, 60)
longest <- orfs[[which.max(nchar(orfs))]]
protein <- translate_dna(longest)
stopifnot(protein == paste(translate(strsplit(tolower(longest), "")[[1]]), collapse = ""))
cat("longest ORF:", nchar(longest), "bp ->", protein, "\n")

## ---- Save summary CSV ----
summary_df <- data.frame(
  seq_id         = record_id,
  length         = nchar(seq),
  gc_content     = round(gc_content(seq), 4),
  n_orfs         = length(orfs),
  longest_orf_bp = max(nchar(orfs)),
  top_3mer       = names(sort(counts, decreasing = TRUE))[1]
)
write.csv(summary_df, "results_R/summary_R.csv", row.names = FALSE)
print(summary_df)

## ---- Figure: GC skew ----
png("results_R/gc_skew_R.png", width = 1400, height = 800, res = 200)
plot(skews, type = "o", xlab = "window index (50 bp)",
     ylab = "GC skew (G-C)/(G+C)", main = paste("GC skew across", accession))
abline(h = 0, col = "grey")
dev.off()
