Developed for the Cambridge iDEC 2026 team, see [the wiki](https://idec-teams.github.io/2026_Cambridge/) for the full project. 

## Sequencing data analyser - `sequencing_data_analyzerV3.py` 

Can be found in [`sequencing_data_analyzer/`](sequencing_data_analyzer/).
Processes raw FASTQ reads from the degenerate libraries (PrCAHS1, SSB, NEXT) into quality-filtered DNA and amino-acid variant counts. 

### Pipeline summary

1. **Parse** FASTQ reads (ID, sequence, quality).
2. **Find the left anchor** with Edlib semi-global alignment (≤16 differences). If it isn't found, retry on the reverse complement. Reads with no anchor in either orientation go to a second pass.
3. **Extract the variable region**, a fixed length after the anchor (57 bp for PrCAHS1-1, 126 bp for SSB, 156 bp for NEXT). Reads too short are discarded.
4. **Check the right-hand constant region** (60 nt, ≤16 mismatches). Reads that fail are excluded because the boundaries can't be trusted.
5. **Rescue pass.** Anchor-failed reads are re-aligned to the reference using a larger tolerance (≤32) and IUPAC-aware matching. The CIGAR string locates the region. Reads with an indel inside the variable region are not rescued.
6. **Compare to the degenerate design** using IUPAC codes at each position.
7. **Quality-aware mismatch handling:** Q ≥ 30 mismatch → read discarded. Q < 30 at a single-nucleotide position → corrected to the expected base. Q < 30 at a degenerate position → read discarded.
8. **Translate** to amino acids and drop reads with premature stop codons. Count each distinct variant.
9. **Optional abundance-based clustering** (greedy; merges a variant into a cluster if the cluster is ≥3× more abundant and Hamming distance ≤2). Downstream analysis uses the raw, unclustered variants.

Full step-by-step logic: [Report which can be accessed from the Wiki](https://idec-teams.github.io/2026_Cambridge/)

### Usage

```bash
python sequencing_data_analyzerV3.py <input.fastq> \
    <dna_variants.csv> <aa_variants.csv> <valid_var.txt> \
    <discarded_at_find_anchor.csv> <discarded_degen_mis.csv>
```

Example:

```bash
python sequencing_data_analyzerV3.py VRZQTP_1_sample_1.fastq \
    cahs_dna_variants.csv cahs_aa_variants.csv valid_var.txt \
    discarded_var_at_find_anchor.csv discarded_degen_mis.csv
```

Currently the anchors and sequences are hardcoded, so they will need to be changed to process other sequences. 

| Argument | Contents |
|---|---|
| `dna_variants.csv` | Reads that passed the anchor checks or were rescued (`read_id`, `strand`, `nt`, `qual`, `source` = direct/rescued). Written **before** the degenerate-sequence filter, so it also includes reads discarded later. |
| `aa_variants.csv` | Raw, unclustered amino-acid variants and their read counts. This is what downstream analysis uses. |
| `valid_var.txt` | Abundance-clustered variants with the percentage of cleanly called reads each cluster accounts for. |
| `discarded_at_find_anchor.csv` | Reads rejected in the first pass (`no_left_anchor`, `short_window`, `indel_in_right_window`). Includes reads that were later rescued. |
| `discarded_degen_mis.csv` | Reads discarded at the degenerate-sequence check (`length`, `confident_mismatch`, `ambiguous_low_q`). |

A read-level QC summary and the top 10 clustered variants are printed to the terminal.

### Requirements
Python 3.7+ and:
```bash
pip install edlib tqdm biopython
```

## `analyzerV4`: read QC and alignment tool 

Can be found in [`sequence_alignment/`](sequence_alignment/).
Quality-control and alignment tool for raw FASTQ reads. We used it to check for **mispriming** by looking at read lengths, 5'/3' primer motifs and alignments to the designed sequence. It was also meant for manually re-analysing reads discarded by the main sequencing analyser (see "sequencing_data_analyzer"), but that was never integrated due to time constraints. 

### What it does
1. Parses the FASTQ and writes per-read stats (length, mean Phred quality) to a CSV.
2. Writes global metrics to a text file: read-length distribution, 5' motif and 3' junction counts, 5'/3' motif combinations, and summed quality by position.
3. Classifies each read by primer orientation and length:
   - starts with the sense primer: kept as is
   - starts with the antisense primer: reverse-complemented first
   - neither: flagged as unknown
4. Aligns each read to the desired sequence using Levenshtein distance and outputs the alignment to a text file. 

### Usage
```bash
python analyzerV4.py input.fastq output.csv metrics.txt alignment.txt
```
The script then prompts for the **desired sequence** (paste the reference sequence that the unknown sequence will be aligned to).

Currently the primers are hard-coded, so they will need to be changed for other sequences.

| Output | Contents |
|---|---|
| `output.csv` | read_id, sequence, quality, length, mean quality |
| `metrics.txt` | global read-length, motif and quality metrics |
| `alignment.txt` | per-read classification and alignments |

### Requirements
`tqdm`, plus the local modules `fastq_parser.py` and `distance.py` (`l_matrix`). 