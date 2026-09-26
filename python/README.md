# DiffDomain-Spectrum — Python implementation

*Enhanced detection of reorganized TADs using sparse aggregated single-cell
Hi-C contact maps.*

Detect whether a topologically associating domain (TAD) is significantly
**reorganized between two conditions**, from sparse aggregated single-cell
Hi-C (scHi-C) contact maps, using a random-matrix-theory (Wigner semicircle)
full-spectrum goodness-of-fit test on the distance-normalized log-ratio
contact matrix.

Given two aggregated contact maps built on the same genomic interval and
resolution, the tool

1. **extracts** the TAD-level contact matrix `M` of each reference TAD from both
   conditions (panel b),
2. **builds and normalizes** the difference matrix `M_diff = log(M1) - log(M2)`
   and standardizes each off-diagonal by its mean/sd to obtain `M_norm`
   (panel c),
3. **tests** the eigenvalue spectrum of `M_norm / sqrt(N)` against the Wigner
   semicircle law with a Monte-Carlo null distribution, then applies a
   Benjamini-Hochberg correction to call **reorganized TADs** (panel d).

![DiffDomain-Spectrum workflow](../assets/workflow.png)

**Figure 1 | DiffDomain-Spectrum workflow.**
**a** Raw scHi-C contact maps are aggregated per condition into one
condition-level contact map. **b** For each reference TAD, the corresponding
TAD-level contact matrices `M1` and `M2` are extracted from the two aggregated
maps over the same genomic interval. **c** The difference matrix
`M_diff = log(M1) - log(M2)` is computed and normalized off-diagonal-wise to
`M_norm`. **d** The eigenvalues of `M_norm / sqrt(N)` are compared with the
Wigner semicircle law through a full-spectrum goodness-of-fit statistic `Z_A`;
a Monte-Carlo p-value is obtained by sampling from the semicircle distribution
and corrected with the Benjamini-Hochberg procedure (`adjusted P <= 0.05`).

> **Scope.** Step **a** (aggregation of raw scHi-C contacts) is performed
> upstream by the single-cell Hi-C processing pipeline. **This tool starts from
> aggregated contact maps** (`.hic` / `.cool`) and implements steps **b-d**.
> Original vector figure: [`workflow.pdf`](../workflow.pdf).

> **Also available in R.** An R implementation of the same method — the
> `DiffDomainSpectrum` package — lives at the **root of this repository**
> (see the [top-level README](../README.md)). Install it with
> `devtools::install_github("FocusPaka/DiffDomain-Spectrum")`.

---

## Repository contents

| File | Purpose |
|---|---|
| `codes_run_spectrum.py` | command line entry point (docopt) |
| `functions.py` | core module: TAD-list loading, matrix extraction from `.hic` / `.cool`, normalization, the Wigner semicircle distribution and the parallel-Monte-Carlo full-spectrum test |
| `environment.yml` | conda environment definition |

---

## Dependencies

Handled by `environment.yml`:

- Python >= 3.9 (3.10 by default)
- `numpy`, `pandas`, `scipy`, `h5py`
- `cooler` — reads `.cool` / `.mcool`
- `hic-straw` (`import hicstraw`) — reads `.hic`
- `matplotlib`, `seaborn` — only for the `visualization` subcommand
- `docopt`, `joblib`, `statsmodels`

---

## Installation

**Important:** this code lives in the `python/` subdirectory of the repository,
so `cd` into it before creating the environment.

```bash
git clone https://github.com/FocusPaka/DiffDomain-Spectrum.git
cd DiffDomain-Spectrum/python
conda env create --name spectrum -f environment.yml
conda activate spectrum
```

Once activated, every dependency is available and the commands below can be run
directly. (Plain `pip install cooler hic-straw numpy pandas scipy h5py docopt
joblib statsmodels matplotlib seaborn` works too if you prefer not to use conda.)

---

## Usage

```
python codes_run_spectrum.py <command> ...

  dvsd one         <chr> <start> <end> <hic0> <hic1> [options]
  dvsd multiple    <hic0> <hic1> <bed>              [options]
  visualization    <chr> <start> <end> <hic0> <hic1> [options]
  adjustment       <method> <input> <output>         [options]
```

`<start>` / `<end>` are genomic coordinates in **bp**; `<reso>` is the bin size
in bp and must match the resolution stored in the Hi-C file.

### 1. Single domain — spectral test

```bash
python codes_run_spectrum.py dvsd one "X" 154425001 154700001 \
    sample_A.hic sample_B.hic --reso 25000 --hicnorm NONE
```

Result is printed to stdout and appended to `--ofile` if given.

### 2. Batched domains

`<bed>` is a tab-separated file whose **first line is a header** (it is skipped
unconditionally) followed by `chrom  start  end  ...` columns:

```bash
python codes_run_spectrum.py dvsd multiple \
    sample_A.cool sample_B.cool tads.bed \
    --reso 25000 --min_nbin 10 --chrn ALL --ncore 10 \
    --ofile results.txt
```

### 3. Visualization

Produces one PDF. The **lower triangle is `hic0`, the upper triangle is `hic1`**:

```bash
python codes_run_spectrum.py visualization "X" 154425001 154700001 \
    sample_A.hic sample_B.hic --reso 25000 --hicnorm NONE \
    --ofile region.pdf
```

### 4. Multiple-testing adjustment

Reads a `dvsd multiple` output file and appends Benjamini-Hochberg (or any
`statsmodels` method) adjusted p-values:

```bash
python codes_run_spectrum.py adjustment "fdr_bh" results.txt results_adj.txt
# add --filter true to keep only adjusted p <= 0.05
```

---

## Input formats

| Extension | Backend | Notes |
|---|---|---|
| `.hic` | `hicstraw` | the requested normalization and resolution must already exist in the file (straw only reads, never computes them) |
| `.cool` | `cooler` | raw counts (`balance=False`); the `--hicnorm` flag is **ignored** for cool input |
| `.mcool` | `cooler` | opened at `::resolutions/<reso>` |
| anything else | pandas | parsed as a 3-column sparse matrix (`bin_i`, `bin_j`, `count`) separated by `--sep` |

`--hicnorm` values for `.hic`: `NONE`, `VC`, `VC_SQRT`, `KR`, `SCALE`, … — use
`NONE` when the file carries no normalization weights.

---

## Output format

`dvsd one` / `dvsd multiple` write two blocks to `--ofile`:

1. the run options, prefixed with `#`;
2. one result line per domain:

```
chr  start  end  region  pvalue  bins  stat
```

where `bins` is the number of retained matrix bins and `stat` is the observed
test statistic `ZA_obs`.

> `adjustment` hard-codes `skiprows=26`, i.e. it assumes the option block above
> the results is exactly 26 lines. If you add or remove a command line option,
> update that number in `codes_run_spectrum.py`.

---

## Options

| Option | Default | Meaning |
|---|---|---|
| `--ofile` | stdout | output file path |
| `--hicnorm` | `KR` | normalization for `.hic` input (ignored for `.cool`) |
| `--chrn` | `ALL` | restrict `dvsd multiple` to one chromosome (without the `chr` prefix) |
| `--reso` | `10000` | bin size in bp |
| `--ncore` | `10` | parallel workers for `dvsd multiple` |
| `--min_nbin` | `10` | minimum number of retained bins for a domain to be tested |
| `--f` | `0.5` | sparsity filter: drop a bin if more than this fraction of its entries are NaN (the R implementation calls this `prop = 1 - f`) |
| `--N` | `10000` | Monte-Carlo replicates for the p-value |
| `--filter` | `false` | for `adjustment`: keep only adjusted p <= 0.05 |

---

## Notes

- The test is resolution- and bin-count sensitive: very sparse intervals are
  reported as `NaN` rows instead of raising an error.
- Running two datasets that are not on a common coordinate system / resolution
  will silently produce meaningless ratios — always verify `--reso` and the
  chromosome naming (`X` vs `chrX`).
- `.hic` and `.cool` inputs are handled by different code paths; results are
  not guaranteed to be bit-identical between the two backends.
- The Python and R implementations are numerically equivalent up to the five
  small differences listed in the [top-level README](../README.md#known-numerical-differences-between-the-two-implementations).

---

## Data availability

Large Hi-C files are **not** stored in this repository. Example data links go
here:

- TODO: dataset A name and download link
- TODO: dataset B name and download link

---

## Citation

If you use this tool, please cite:

- **DiffDomain-Spectrum: enhanced detection of reorganized TADs using sparse
  aggregated single-cell Hi-C contact maps** — TODO: add authors, journal /
  preprint server, year and DOI.
- Cooler: Abdennur N., Mirny L.A. (2020). *Cooler: scalable storage for Hi-C
  data and other genomically labeled arrays.* Bioinformatics 36(1):311–316.
  doi:10.1093/bioinformatics/btz540
- Straw: Durand N.C. et al. (2016). *Juicebox provides a visualization system
  for Hi-C contact maps with unlimited zoom.* Cell Systems 3(1):99–101.

---

## License

Released under the **GNU General Public License, version 2 (GPL-2)** — see
[`LICENSE`](../LICENSE). The Python implementation is a port of the R package in
this repository and carries the same license. If any code was adapted from
another project, that project's upstream license must be preserved as well.

---

## Contact

TODO: add name and email / homepage.
