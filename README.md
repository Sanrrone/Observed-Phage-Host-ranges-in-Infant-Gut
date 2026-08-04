# Observed phage–host ranges in the infant gut

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21218558.svg)](https://doi.org/10.5281/zenodo.21218558)

Sequence data and host-assignment tables supporting the study of phage–bacteria
interactions in the early-life gut microbiome. This repository contains the phage
contigs analysed in the paper, the per-sample bacterial genomes (MAGs) they were
mapped against, and the resulting observed phage–host associations across three
reference sets.

All identifiers have been anonymised: infants are labelled `B1 … B41` and each
biological sample is named `B<n>_<timepoint>` (e.g. `B1_6_months`). No file name,
sequence header, or table field contains the original participant IDs.

---

## Dataset at a glance

| Item                                               | Count                       |
| -------------------------------------------------- | --------------------------- |
| Infants                                            | 41                          |
| Sample-timepoints (1, 6 and 12 months)             | 55                          |
| Faecal metagenomes                                 | 1,366                       |
| Phage contigs (≥ 10 kb)                            | 13,235                      |
| Bacterial genomes / MAGs                           | 1,688 (≈ 31 per sample)     |
| Host assignments — HumGut DB (HGDB)                | 8,764 hits · 2,236 contigs  |
| Host assignments — self assemblies (SBA)           | 5,597 hits · 5,595 contigs  |
| Host assignments — infant-cohort assemblies (GIBA) | 22,406 hits · 7,043 contigs |

Timepoint breakdown: 4 samples at 1 month, 21 at 6 months, 30 at 12 months.

---

## Repository structure

```
.
├── phage_contigs/          # 55 FASTA files, one per sample (B<n>_<age>_phages.fna)
│                           #   → 13,235 phage contigs > 10 kb in total
├── bact_contigs/           # 55 sample folders (B<n>_<age>/)
│   └── B<n>_<age>/
│       ├── genomes/        #   host bacterial genomes / MAGs (<taxid>_<species>.fna)
│       └── B<n>_<age>_checkm2.tsv   # CheckM2 completeness/contamination report
├── host_ranges/            # observed phage → host tables (one per reference set)
│   ├── host_hgdb.tsv       #   mapping vs the HumGut reference database
│   ├── host_sba.tsv        #   mapping vs same-sample bacterial assemblies
│   └── host_giba.tsv       #   mapping vs the pooled infant-cohort assemblies
├── bacterial_taxonomy_profiles.tsv   # per-metagenome bacterial taxonomy (MetaPhlAn4)
├── slurm_scripts/          # the SLURM pipeline used to generate the data
│   └── supp_files/         # anonymised supplementary tables (see below)
├── LICENSE                 # MIT
└── README.md
```

---

## File formats

### `phage_contigs/B<n>_<age>_phages.fna`

Standard nucleotide FASTA. Each record is a phage contig ≥ 10 kb. Header is the
contig ID, e.g.

```
>k141_11970
>k141_11987_phage171_86355_120697_42_0.74      # predicted prophage region on a host contig
```

### `bact_contigs/B<n>_<age>/genomes/<taxid>_<species>.fna`

Nucleotide FASTA for each recovered host genome / MAG, named by NCBI taxid and
species. Sequence headers encode `sample__taxid_species__contig`, e.g.
`>B1_6_months__1128111_Veillonella_atypica__k141_22086`. Genome quality
(completeness, contamination, size, N50, …) is reported per sample in the
accompanying `B<n>_<age>_checkm2.tsv` (CheckM2 output).

### `host_ranges/*.tsv`

Tab-separated observed phage–host associations, one row per alignment. Columns:

| Column                       | Description                                                                                                                                     |
| ---------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`                     | Anonymised sample (`B<n>_<age>`) the phage contig comes from                                                                                    |
| `contig`                     | Phage contig ID                                                                                                                                 |
| `qlen`                       | Phage contig length (bp)                                                                                                                        |
| `HumGut_name`                | *(host_hgdb.tsv only)* matched HumGut genome                                                                                                    |
| `ref`                        | Matched host reference. For SBA/GIBA this encodes `donor_sample__taxid_species__host_contig`; for HGDB it is the HumGut/kraken reference string |
| `sstart`, `send`             | Alignment coordinates on the host reference                                                                                                     |
| `host_genus`, `host_species` | Taxonomic assignment of the matched host                                                                                                        |

The three tables differ only in the reference set the phage contigs were mapped
against:

- **HGDB** — the [HumGut](https://doi.org/10.1186/s40168-021-01114-w) reference
  collection of human-gut prokaryotic genomes.
- **SBA** (self bacterial assemblies) — bacterial genomes assembled from the
  **same** sample as the phage; `sample` and the donor in `ref` always match.
- **GIBA** (gut infant bacterial assemblies) — bacterial genomes assembled across
  the **whole infant cohort**, so a phage from one infant can be matched to a host
  MAG recovered from another infant (the donor is given in `ref`).

### `bacterial_taxonomy_profiles.tsv`

Per-metagenome bacterial taxonomic relative-abundance profiles (MetaPhlAn4), long
format — one row per taxon per metagenome. Columns: `sample`, the seven rank
columns (`kingdom … species`, `NA` where unresolved) and `rel_abu` (relative
abundance, %). Keyed on the same anonymised `sample` IDs (`B<n>_<k>`) as
`slurm_scripts/supp_files/supp_table1.tsv`, so the two join directly; all 1,366
study metagenomes are covered.

---

## Methods (summary)

Phage contigs were identified from per-sample metagenomic assemblies (including
predicted prophage regions) and retained at ≥ 10 kb. Host ranges were determined
by direct nucleotide mapping of each phage contig against the three reference sets
above, keeping alignments passing the identity and coverage thresholds described
in the manuscript. Bacterial host genomes were assembled and binned per sample and
their quality assessed with CheckM2. The full processing pipeline is provided in
`slurm_scripts/` (numbered `0_*` → `13_*`).

## Containers and pipelines

All analyses were executed within a Singularity (SIF) container to ensure
reproducibility across computing environments.

- **Bacteria mining pipeline** — container:
  <https://datacloud.helsinki.fi/index.php/s/4SX7wmZBttpnWRg>
- **Phage mining pipeline** — container:
  <https://datacloud.helsinki.fi/index.php/s/JfFByXAwgYbBkQ9>

More information about the pipelines can be obtained from
[Sanrrone/Gutbusters](https://github.com/Sanrrone/Gutbusters).

## Supplementary metadata (`slurm_scripts/supp_files/`)

- `supp_table1.tsv` — per-metagenome metadata (sample, day, age, feeding/antibiotic
  protocol, delivery mode); 1,366 rows.
- `allowed_10k.tsv` — the list of the 13,235 phage contigs (≥ 10 kb) per sample.
- `hep_sa.txt`, `hep.txt` — anonymised sample and infant identifier lists.

---

## Citation

If you use these data, please cite the associated article:

> *Phage–bacteria interactions in the infant gut.* **Microbiology Spectrum** (2026).

*(Article DOI https://doi.org/10.1128/spectrum.00361-26)*

The dataset itself is archived on Zenodo; please also cite it via its DOI
[10.5281/zenodo.21218558](https://doi.org/10.5281/zenodo.21218558) (see
`CITATION.cff` for the full metadata).

## License

Released under the [MIT License](LICENSE).
