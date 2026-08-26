# Sequence Analysis Toolkit

A small toolkit for basic DNA sequence analysis, written from first principles in **Python and R** and
validated against established libraries. I use the human insulin (*INS*) mRNA as a worked example, fetched
live from NCBI, and compute its base composition, coding potential, and protein product.

![GC skew across NM_000207](results_py/gc_skew.png)

## 1. Research question

For a real human gene, can I compute the standard sequence properties — GC content, GC skew, k-mer
frequencies, reverse complement, open reading frames (ORFs), and protein translation — from scratch, and
confirm they agree with reference implementations?

## 2. Biological background

DNA is a string over A, C, G, T, but each operation has a biological meaning: GC content reflects base
composition (and affects melting temperature); GC skew relates to strand asymmetry; ORFs mark candidate
protein-coding regions (ATG → in-frame stop); and translation converts codons to amino acids. Insulin is a
good test case because its mRNA contains a single, well-characterised coding sequence (preproinsulin).

## 3. Dataset

- **Source:** NCBI Entrez API (RefSeq accession **`NM_000207`**, human insulin mRNA).
- **Access:** fetched programmatically in code — `Bio.Entrez` (Python) / `rentrez` (R) — not a manual download.
- **Size:** one transcript, 465 nucleotides.

## 4. Methods

Each function is implemented by hand and cross-checked against a library:

| Task | From scratch | Cross-check |
|------|--------------|-------------|
| Parsing | custom FASTA parser | Biopython `SeqIO` / R `readLines` |
| GC content, skew, k-mers | base Python / base R | seqinr (R) |
| Reverse complement | dictionary map + reverse | Biopython / seqinr |
| Translation | codon-table lookup | Biopython `Seq.translate` / seqinr `translate` |

Python was the primary language (Biopython makes sequence handling clean); the analysis is reproduced in R
to confirm the results are implementation-independent.

## 5. Analysis pipeline

1. Fetch the mRNA from the Entrez API and cache it locally.
2. Parse the FASTA (custom parser, verified against the library).
3. Compute GC content and windowed GC skew.
4. Count k-mers.
5. Compute the reverse complement (checked with the involution `revcomp(revcomp(x)) == x`).
6. Find ORFs in the three forward frames.
7. Translate the longest ORF and compare to the library translation.
8. Write a summary CSV and the GC-skew figure.

## 6. Results

| metric | value |
|--------|-------|
| length | 465 bp |
| GC content | 0.639 |
| distinct 3-mers | 60 (most frequent: `CTG`) |
| ORFs (forward, ≥60 bp) | 1 |
| longest ORF | 333 bp |

Full table in [`results_py/summary.csv`](results_py/summary.csv).

## 7. Figures

`results_py/gc_skew.png` — GC skew across the transcript in 50 bp windows (shown above).

## 8. Limitations

- ORF detection covers the **three forward frames only** — reverse-strand ORFs are not searched.
- This is a **single sequence**, so nothing here generalises.
- GC skew is really a genome / replication-origin–scale signal; on a 465 bp mRNA it is included for
  completeness rather than biological meaning.

## 9. Biological interpretation

The single 333 bp ORF translates to a 110-residue protein beginning `MALWMRLLPLL...` and ending in a stop —
this is **preproinsulin** (signal peptide, B-chain, C-peptide, A-chain). The from-scratch translation matches
Biopython and seqinr exactly, and the ~0.64 GC content is consistent with a human coding transcript. So the
toolkit correctly recovers a real, annotated protein from raw mRNA.

## Language comparison (Python vs R)

The toolkit is implemented twice — **Python** (Biopython) and **R** (seqinr) — on the same NCBI record.
The results are **identical**: GC content 0.6387, one 333 bp forward ORF, top 3-mer `CTG`, and the same
110-residue preproinsulin translation (each cross-checked against its own library). Matching outputs from
two independent implementations is a reproducibility check — the biology doesn't depend on one toolchain.
The differences are idiomatic, not numerical: R is 1-indexed with inclusive ranges, uses named vectors as
lookup tables, and needed care with operator precedence in the GC-skew formula `(g - c) / (g + c)`.

## 10. Reproduce

```bash
# Python (Google Colab or local)
pip install biopython matplotlib
# open sequence_analysis_toolkit.ipynb and run all cells

# R (RStudio)
install.packages(c("rentrez", "seqinr"))
# open sequence_analysis_toolkit.R in RStudio and Source it
```

Data is pulled live from NCBI (accession `NM_000207`); change the accession to run the toolkit on any other
gene.

## 11. Abstract

*Implemented a DNA sequence-analysis toolkit in Python and R that fetches gene records from the NCBI Entrez
API and computes GC content, GC skew, k-mer spectra, reverse complements, ORFs, and protein translation from
first principles, validated against Biopython and seqinr. Applied to the human insulin mRNA, it recovers the
preproinsulin ORF and reproduces the reference translation exactly. Demonstrates API-based data access,
biological-sequence handling, and reproducible, cross-validated code.*

## Tech

`Python` (Biopython, matplotlib) · `R` (rentrez, seqinr) · NCBI Entrez API

## Files

```
sequence_analysis_toolkit.ipynb   # Python implementation
sequence_analysis_toolkit.R       # R implementation (RStudio)
results_py/summary.csv            # Python: per-sequence metrics
results_py/gc_skew.png            # Python: GC-skew figure
results_R/summary_R.csv           # R: per-sequence metrics
results_R/gc_skew_R.png           # R: GC-skew figure
```
