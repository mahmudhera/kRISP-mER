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

On-target scores are produced by the Azimuth 2.0 implementation of
Doench's Rule Set 2, as bundled with CRISPOR.

## Why not the bundled model

`krispmer/on_target_scoring/` contains `V3_model_nopos.pickle`. This is a
Python 2, protocol 2 pickle that references `PresortBestSplitter`, a Cython
class removed from scikit-learn after version 0.17. It cannot be unpickled
on a modern scikit-learn, so that code path is unusable. Azimuth 2.0 ships
its own model that loads correctly.

## Obtaining Azimuth

Azimuth 2.0 is distributed inside CRISPOR:

```bash
git clone https://github.com/maximilianh/crisporWebsite
```

The scorer looks for `bin/Azimuth-2.0` under, in order:

1. `$AZIMUTH_DIR`, if set
2. `~/crisporWebsite/bin/Azimuth-2.0`
3. `~/Crispor/crisporWebsite/bin/Azimuth-2.0`

If CRISPOR is cloned elsewhere, set the environment variable:

```bash
export AZIMUTH_DIR=/path/to/crisporWebsite/bin/Azimuth-2.0
```

## Why both paths go on sys.path

Azimuth's `util.py` uses Python 2 implicit relative imports. Under Python 3
these resolve only if both the package directory and its parent are on
`sys.path`. `rs2_azimuth.py` appends both; this is not redundant.

## Verification

Run the scorer directly:

```bash
python -m krispmer.rs2_azimuth
```

Expected output:

```
reference 30-mer: ACGCAGCGCCTCCCTCCACTCGGAAGGACT
score           : 55
expected        : 55  (CRISPOR doctest)
MATCH
```

A score of 55 confirms the integration matches CRISPOR's own expected
value for this sequence.

## Behaviour when Azimuth is absent

If Azimuth cannot be found or imported, kRISP-meR logs a warning and writes
`-1` in the `on_target_score` column. The run does not fail. A `-1` is also
written when no 30 bp context could be extracted, which happens for guides
closer than 4 bp to the target boundary.

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
