# INSTALLATION GUIDE — kRISP-mER

* Use a Conda virtual environment
* Install dependencies from default channels
* Build **Jellyfish** from the official upstream GitHub
* Manually place **dna_jellyfish** into the environment's `site-packages`

This is the **authoritative installation procedure**.

---

# 1. Create Conda Environment

```bash
conda create -n krispmer python=3.8 -y
conda activate krispmer
```

---

# 2. Install Required Tools (samtools, bowtie2, Java)

These are available from default Conda channels. No special configuration needed.

```bash
conda install -y samtools bowtie2 openjdk
```

Verify:

```bash
samtools --version
bowtie2 --version
java -version
```

---

# 3. Install Python Dependencies

Install common scientific Python packages:

```bash
conda install -y numpy scipy pandas biopython scikit-learn
```

---

# 4. Build & Install Jellyfish (from official GitHub)

kRISP-mER requires the Jellyfish binary. Install it from upstream:

**Repository:** [https://github.com/gmarcais/Jellyfish](https://github.com/gmarcais/Jellyfish)

### 4.1 Clone the repository

```bash
cd ~
git clone https://github.com/gmarcais/Jellyfish.git
cd Jellyfish
```

### 4.2 Install build dependencies into the conda environment

```bash
conda install -y autoconf automake libtool make gxx_linux-64 pkg-config
```

### 4.3 Build Jellyfish into your Conda environment

```bash
autoreconf -i
./configure --prefix="$CONDA_PREFIX"
make -j 4
make install
```

### 4.4 Verify

```bash
jellyfish --version
```

---

# 5. Install kRISP-mER

Clone the repository:

```bash
cd ~
git clone https://github.com/mahmudhera/kRISP-mER.git
cd kRISP-mER
```

Install in editable mode:

```bash
python setup.py install
```

---

# 6. Copy dna_jellyfish (Required Step)

After building Jellyfish from the upstream repository, the Python binding
`dna_jellyfish` is available inside the Jellyfish source tree (typically
under a `python/` directory or generated during build).

kRISP-mER requires `import dna_jellyfish` to work. You must manually copy
this module into the active Conda environment's `site-packages`.

### 6.1 Locate Conda site-packages

```bash
python -c "import sysconfig; print(sysconfig.get_paths()['purelib'])"
```

Example output:

```
/home/user/miniconda3/envs/krispmer/lib/python3.8/site-packages
```

### 6.2 Locate dna_jellyfish inside Jellyfish build

From the Jellyfish repository you built earlier:

```bash
cd ~/Jellyfish
find . -name "dna_jellyfish*"
```

This should locate the `dna_jellyfish` directory or compiled module.

### 6.3 Copy dna_jellyfish into Conda environment

```bash
cp -r path/to/dna_jellyfish \
$(python -c "import sysconfig; print(sysconfig.get_paths()['purelib'])")/
```

Ensure the copied folder contains `__init__.py` or compiled extension files.

### 6.4 Verify import

```bash
python -c "import dna_jellyfish; print('dna_jellyfish OK')"
```

---

# 7. Final Verification

Final Verification

```bash
python -c "import dna_jellyfish, numpy, sklearn; print('Imports OK')"
jellyfish --version
krispmer --help
```

If all three commands work, installation is complete.

---

# 8. Example Usage

```bash
krispmer -nvr reads.fastq target.fasta out.csv 1
```

Ensure the input FASTQ and FASTA files exist.

---

# 9. Notes

* This guide intentionally avoids additional conda channels.
* Jellyfish is always installed from the upstream GitHub repository.
* dna_jellyfish must be copied manually because it is not available from pip or conda.

---



# 10. Rule Set 2 On-Target Scoring (Azimuth 2.0)

On-target scores come from Azimuth 2.0, the reference implementation of
Doench's Rule Set 2. It is included at `vendor/Azimuth-2.0`, so there
is nothing to install. Verify with:

```bash
python -m krispmer.rs2_azimuth
```

This scores the reference 30-mer `ACGCAGCGCCTCCCTCCACTCGGAAGGACT` and
should print 55, the value given in CRISPOR's own doctest.

Azimuth 2.0 is BSD-3-Clause, copyright 2015 Microsoft Research. Its
licence travels with it at `vendor/Azimuth-2.0/LICENSE.txt`. If you
report Rule Set 2 scores, cite Doench et al., *Nature Biotechnology*
34:184-191, 2016.

**Why the model inside `krispmer/on_target_scoring/` is not used.** That
directory holds `V3_model_nopos.pickle`, a Python 2 protocol 2 pickle
referencing `PresortBestSplitter`, a Cython class dropped from
scikit-learn after version 0.17. It cannot be unpickled on a current
scikit-learn, so that code path is dead. The vendored Azimuth copy carries
a Python 3 re-pickled model that loads.

**Why `rs2_azimuth.py` adds two paths to `sys.path`.** Azimuth's `util.py`
uses Python 2 implicit relative imports, which under Python 3 resolve only
when both the package directory and its parent are on the path. Appending
both is deliberate, not redundant.

**Using a different copy.** The scorer checks `$AZIMUTH_DIR` first if it is
set, then the vendored directory, then `~/crisporWebsite/bin/Azimuth-2.0`.

**If Azimuth cannot be imported**, kRISP-meR logs a warning and writes `-1`
in the `on_target_score` column rather than failing. A `-1` also appears
when no 30 bp context could be extracted, which happens for guides within
4 bp of the target boundary.

---

# 11. Output Columns

The output CSV carries a two-line header comment giving the estimated
target coverage, followed by nine columns:

| Column | Meaning |
|---|---|
| `tgt_in_plus` | the guide sequence, including PAM |
| `tgt_in_minus` | its reverse complement |
| `inverse_specificity` | E[x_G] / (E[x_t] x target_coverage) |
| `Expected_cuts_in_genome` | E[x_G] |
| `Expected_presence_in_target` | E[x_t] |
| `gc_content` | GC fraction of the guide |
| `strand` | `+` or `-` |
| `on_target_score` | Rule Set 2, 0-100 |
| `context_30bp` | the 30 bp window scored by Rule Set 2 |

`Expected_cuts_in_genome` is the quantity guides are ranked and thresholded
on. A candidate is accepted when it falls below 1.5.

`-1` appears in `inverse_specificity` when the target coverage or the
expected presence in target is zero, and in `on_target_score` as described
above. In both columns `-1` means "not computed", not a low score. Filter
these rows out rather than treating them as values.
