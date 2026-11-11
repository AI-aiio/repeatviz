# repeatviz

Chromosome-layout heatmaps of **tandem** and **multi-copy** repeats for draft-assembly
curation: *where are the repeats, what kind are they, and are they where the contigs break?*

`repeatviz` is meant for the stage right after automatic scaffolding (e.g. synscaf / YaHS)
and before manual curation. Every chromosome is drawn as a bar of contigs (contig boundaries
marked). `plot` writes **four PDFs** (`--maps` to choose) (portrait, ~15 chromosomes per page, then UNPLACED pages):

| PDF | heatmap shows |
|---|---|
| `repeatviz_pdf/<prefix>.all.pdf` | **All repeats**: the larger of k-mers occurring ≥10× in the assembly or tandem arrays |
| `repeatviz_pdf/<prefix>.tandem.pdf` | **Tandem repeats, short units** (<100 bp, incl. telomere) |
| `repeatviz_pdf/<prefix>.long.pdf` | **Tandem repeats, long units** (≥100 bp, satellites) |
| `repeatviz_pdf/<prefix>.dup.pdf` | **Duplications**: k-mers occurring 2–9× (segmental duplications, recent paralogs, collapsed/alt haplotype copies) – the other main reason HiFi contigs break |

Each chromosome row:

| element | what it shows | source |
|---|---|---|
| **contig bar = repeat heatmap** | per 100 kb window (500 kb moving average for display). Diverging colours: **blue = repeat-poor**, grey = genome median (tandem maps: at least 2 %), **red = repeat-rich**, dark red = 100 % | meryl + ULTRA / k-mer long pass |
| contig extent | double arrow ◀ name ▶ above each contig, white gap + thin black outline at boundaries | layout |
| **depth wave** *(optional, `--depth`)* | same height as the bar; read depth relative to the genome median (dotted = median, top = 3×); red fill ≥ `--collapse`, blue fill ≤0.5×; where depth ≥ `--depth-warn` (2×) a red ▼ marks it with the peak as a small red label, e.g. `high: max 110x` (peaks on both sides of a junction are merged) | mosdepth / any bedGraph |
| **yellow flag** | telomere array ≥300 bp (also *interstitial* ones); the pennant points away from the chromosome | TTAGGG/CCCTAA motif density (v0.3) |

Unit lengths, low-copy (≥2×) fractions and per-class tandem fractions are not drawn but are in
`windows.tsv.gz` / `arrays.tsv`. **Alerts** (windows mostly high-copy with little tandem
annotation: an SD / interspersed-repeat cluster, or an array the long pass could not resolve) are
always written to `alerts.tsv`; `--alerts` also draws them as small orange marks above the bar.

Outputs are a multi-page PDF (+ one PNG per page; unplaced contigs on their own pages) and tables:

| file | content |
|---|---|
| `arrays.tsv` | every tandem array: contig, start, end (0-based, half-open), length, period, copies, class, consensus unit |
| `units.fa` | consensus unit of every array ≥ 10 kb (e.g. to look up in Dfam / RepeatMasker) |
| `windows.tsv.gz` | the per-window track values |
| `<prefix>.contig_ends.tsv` | for **each contig end** (default 200 kb): role (`chr_end` / `junction` / `unplaced`), bp per tandem class, top unit, multi-copy fractions, depth flags and a `category` (`telomere`, `satellite`, `short-tandem`, `high-copy`, `low-copy/SD`, `collapse`, or **`no repeat`**) |
| `<prefix>.alerts.tsv` | merged alert windows |
| `<prefix>.arrays.layout.tsv` | `arrays.tsv` + chromosome coordinates in the layout |

## How it works

`repeatviz run` (compute; needs `meryl`, `meryl-lookup`, `ultra`, `tantan`):

1. `meryl count k=21` on the assembly, `greater-than` 1 and 9 → `meryl-lookup -bed-runs` → bp per window covered by ≥2× / ≥10× k-mers.
2. **ULTRA, short pass** (`-p 100`): micro/minisatellites and telomeres. By default (`--short-mode tantan`)
   ULTRA only sees the windows that [tantan](https://gitlab.com/mcfrith/tantan) (a fast tandem-repeat
   masker, `-w 100`) flags, padded by 100 bp and merged – typically a few % of a mammalian genome, so
   several times faster than the whole assembly, while keeping nearly all of the bp of ULTRA arrays ≥ 500 bp and
   the telomere arrays (what is lost: weak low-complexity calls of a few hundred bp). `--short-mode kmer` uses exact
   k-mer recurrence instead (no extra tool, but misses diverged arrays); `--short-mode full` = v0.1.
3. **Long pass** only on *long, dense* high-copy regions (≥10 kb, ≥80 % covered by ≥10× k-mers, + padding),
   chosen with `--long-mode`:
   - **`kmer` (default, light).** For every position, the distance to the next occurrence of the same
     k-mer (k=16) is computed (one sort per region). Inside a tandem array with unit P most distances are
     P (or 2P, 3P where a copy is mutated). Each 5 kb tile whose distance mode explains ≥30 % of it is
     tandem; tiles of similar period are merged into arrays and trimmed to the first/last supported
     position. There is **no maximum period** (40 kb units are found like 171 bp ones). Other,
     non-harmonic distance peaks are reported as `period2` / `period3` (monomer of a HOR, mixed arrays).
     Then **ULTRA is run on one window** (max(10 kb, 4 × period)) from the middle of every array
     ≥ 50 kb with period ≤ `--p-long`, with `-p` just above the period, only to get a consensus unit
     (`--no-sample-ultra` skips it; units are then one raw copy from the middle of the array).
     On a mammalian genome this takes a few minutes. Compared with ULTRA `-p 2000`, most k-mer array bp
     get the same period or a harmonic of it; units > 2 kb are found only by the k-mer pass.
   - **`ultra` (detailed, slow).** ULTRA `-p 2000` on the whole regions (v0.1 behaviour): exact
     boundaries and a consensus for every array, but hours on a mammalian genome
     (much longer on a small workstation) and no units > 2 kb.
   - **`none`.** Short pass only.
4. ULTRA reports long arrays in ~10–50 kb pieces; adjacent calls with a similar period (±2 %) are merged into one array. Telomeres (v0.3): the assembly is scanned for exact TTAGGG / CCCTAA (like `tidk search`); hits ≤18 bp apart are chained, and a run ≥100 bp that is ≥50 % exact motif is a telomere run (`telomere_runs.tsv`). Any array ≥50 % covered by such runs is labelled `telomere` whatever ULTRA's period (6/12/18) or variant consensus (TTGGGG-type, …); runs no array covers are added (`method = telo_scan`). Up to v0.2 only a 1×/2× TTAGGG consensus counted, which missed many chromosome ends that tidk calls.

`arrays.tsv` columns added in v0.2: `method` (`ultra`, `kmer`, `kmer+ultra` = k-mer extent and period,
ULTRA consensus), `tandem_frac` (share of the array explained by its period; low = diverged or mixed),
`period2`, `period3`, `ultra_period` (period ULTRA found in the sampled window, also when it disagrees).

The number of ULTRA threads is chosen from `--memory` using `ultra --mem` (the long pass needs a few GB per thread). Each step's output is kept, so an interrupted run resumes.

**Run time (v0.2 defaults)**: for a ~3 Gb mammalian assembly on 8 threads, roughly **30 min** in total
and under 10 GB of memory; the ULTRA short pass and meryl take most of the time (v0.1 took hours).

`repeatviz plot` is pure Python (numpy / pandas / matplotlib) and only reads the files written by `run`.

## Install

Works on any Linux or macOS machine; an HPC cluster is **not** required.

```bash
conda env create -f environment.yml      # python + meryl + ULTRA + tantan from bioconda
conda activate repeatviz
pip install -e .
```

Without conda: `pip install -e .` for the Python part, and put `meryl`, `meryl-lookup`, `ultra` and `tantan`
on `PATH` — or point repeatviz to containers:

```bash
export REPEATVIZ_MERYL="apptainer exec /path/meryl.sif meryl"
export REPEATVIZ_MERYL_LOOKUP="apptainer exec /path/meryl.sif meryl-lookup"
export REPEATVIZ_ULTRA="apptainer exec /path/ultra.sif ultra"
export REPEATVIZ_TANTAN="apptainer exec /path/tantan.sif tantan"
```

(`--meryl-cmd`, `--meryl-lookup-cmd`, `--ultra-cmd`, `--tantan-cmd` do the same per run.) `scripts/run_repeatviz.slurm.example.sh` is a SLURM wrapper for clusters with apptainer images.

## Usage

```bash
repeatviz run  assembly.fa -o rv -t 16 -m 64
repeatviz plot rv -o rv/fig --edits edits.tsv --depth depth.regions.bed.gz     # curation plan
repeatviz plot rv -o rv/synscaf --agp synscaf.agp                              # or an AGP
repeatviz plot rv -o rv/raw                                                     # no layout: contigs >= 10 Mb as rows

# v0.2 defaults are the light ones (--short-mode tantan --long-mode kmer). Slower, more detail:
repeatviz run  assembly.fa -o rv -t 16 -m 64 --long-mode ultra                # ULTRA consensus for every long array
repeatviz run  assembly.fa -o rv -t 16 -m 64 --short-mode full                # ULTRA short pass on the whole assembly
```

Layouts: an **AGP** (objects whose name matches `--chrom-regex`, default `chr1`…/`scaf_chrX`, are
chromosomes; everything else is unplaced) or a curation **edits.tsv** (`object part contig start end
orient gap_after status note`, 1-based inclusive; objects `UNPLACED` / `REMOVED`). Chromosomes whose
plan still has `draft` rows are labelled in red. Contigs missing from the layout are listed as unplaced.

Demo (≈1 min): `python3 demo/make_demo.py demo/ && repeatviz run demo/demo.fa -o demo/out/run --window 10000 --region-min-len 5000 -t 4 -m 16 && repeatviz plot demo/out/run -o demo/out/demo --edits demo/demo.edits.tsv --depth demo/demo.depth.bed --end-flank 50000`

Tests: `pytest` (no external tools needed).

## Limitations

- **Interspersed repeats are not named.** LINE / SINE / LTR / DNA elements only appear as
  "high-copy" or "low-copy" k-mer signal; run RepeatMasker separately if you need classes.
- **Repeats missing from the assembly are not shown.** Satellite that the assembler left out
  (inside gaps) is invisible here; a comparison with read k-mers is planned as a separate output.
- **Long pass, `kmer` mode:** periods come from exact 16-mer repeats, so strongly diverged
  arrays (old satellites, < ~85 % identity between neighbouring copies) can be missed or split;
  a HOR is reported with its HOR length as `period` (monomer, if visible, in `period2`); boundaries
  are exact only to the first/last repeated k-mer, and mixed arrays are cut at 5 kb tile resolution.
  Two-copy tandem duplications (period = SD length) are reported as arrays with ~2 copies.
  Use `--long-mode ultra` on a cluster when a consensus for every array is needed.
- **Tandem units longer than `--p-long`** get no ULTRA consensus (raw copy instead); with
  `--long-mode ultra` / `none` they are not annotated at all, only flagged by an alert.
- **Segmental duplications:** the low-copy track says *that* a region is duplicated, not *where*
  the other copy is or in which orientation. (A self-alignment track is optional / planned.)
- **Non-repeat break causes** (heterozygous SVs that split haplotypes, coverage drops,
  assembler graph issues) are not diagnosed; such contig ends are reported as `no repeat`.
- **Collapsed repeats** are only visible with the optional depth track; without reads a collapsed
  array looks like a normal-copy region.
- The k-mer tracks are computed on the assembly as given: in a primary (pseudo-haploid)
  assembly with retained haplotigs, haplotig/primary pairs show up as low-copy.

## License

MIT License. See [LICENSE](LICENSE).
