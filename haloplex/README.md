# HaloPlex Cardio panel — re-analysis pipeline

Task: **HaloPlex samples re-analysis and annotation with ANNOVAR.**

This directory adapts the WES pipeline in `../new_version/` for Agilent **HaloPlex**
targeted data from the 96-gene cardio panel. It is a separate directory rather than
a patch because three of the WES pipeline's steps are not merely suboptimal here —
they are actively wrong, and silently so.

---

## The panel

| | |
|---|---|
| Design name | `Cardio_Panel_96-genes` |
| Design ID | `62192-1551680838` |
| Chemistry | Agilent HaloPlex |
| Genome build | **hg19 / GRCh37** |
| Run date | 2018-05-17 |
| Finalized | 2019-03-04 |

### There are TWO designs on the NAS. Use the right one.

| Design | Regions | Total bp | refGene genes |
|---|---|---|---|
| **`22408-1387190209_Regions.bed`** ✅ | 2007 | **463 767** (464 kb) | **111** |
| `ZMF96cardio_1_*` ❌ | 1425 | 351 759 | **14** |

Both carry an "Agilent HaloPlex — ZMF96cardio" track name, and `ZMF96cardio_1`
sits one directory *deeper* (`Ulykbek_CardioPanel_96/Haloplex_merged_fastq_gz/
ZMF96cardio_1/`) than the correct one, which is at the top of
`Ulykbek_CardioPanel_96/`. It is very easy to grab the wrong one — this pipeline
did, on its first run.

**How the mistake was caught, and how to catch it again:** ask the *data*, not the
filenames. Take one aligned BAM and see where coverage actually is:

```bash
samtools depth sample.bam | awk '$3>=20' | wc -l          # bases at >=20x
```

A sample BAM has **818.7 kb** at ≥20×. Only **41.9%** of that falls inside
`ZMF96cardio_1`; the other **476 kb across 166 genes** — `MYBPC3`, `MYH7`, `PKP2`,
`DSP`, `KCNQ1`, `TNNT2`, `DES`, `CACNA1C` — lay outside it. Calling restricted to
that BED would have silently dropped the most clinically important
cardiomyopathy genes while producing perfectly healthy-looking output.

If the covered footprint and the BED disagree by more than a rounding error, the
BED is wrong. `runs.txt` is no help either: it names design `62192-1551680838`,
finalized 2019-03-04 — *after* these samples were sequenced in 2016–2018, so it
cannot be the design used.

`22408-1387190209` ships no separate amplicons track, so coverage QC uses the same
regions file as calling. Chromosomes touched: 1–8, 10–12, 14–20, 22, X.

**464 kb still drives every decision below.** It is about 1/150th of an exome —
small enough that BQSR and VQSR remain unusable, as argued in each section.

---

## What changed versus the WES pipeline, and why

### 1. `MarkDuplicates` — removed

**Why.** HaloPlex fragments DNA with restriction enzymes, not by random shearing.
Every read covering a given amplicon therefore starts and ends at the *same*
coordinates by construction.

Positional deduplication rests on the assumption that two reads sharing coordinates
are overwhelmingly likely to be PCR copies of one molecule. For randomly sheared
libraries that holds. For HaloPlex it is false by design — the collision is
guaranteed. Running `MarkDuplicates` typically flags 70–90% of reads and collapses
coverage to a fraction of what was sequenced.

**Generalisation worth remembering:** positional dedup is valid only when
fragmentation is *random*. Hybrid capture (SureSelect, SeqCap) and WGS — yes. Every
amplicon chemistry — HaloPlex, AmpliSeq, QIAseq, Fluidigm — no.

**Exception:** HaloPlex**HS** carries UMIs. That data *must* be deduplicated, but by
molecular barcode (Agilent AGeNT `LocatIt`), never by coordinate. These samples are
classic HaloPlex without UMIs, so there is nothing to deduplicate at all.

**How to verify on your own data** rather than taking this on faith:

```bash
samtools view sample.bam | awk '{print $3":"$4}' | sort | uniq -c | sort -rn | head -20
```

Thousands of reads at identical positions is the restriction-fragmentation
signature, not a bad library.

### 2. BQSR (`BaseRecalibrator` / `AnalyzeCovariates` / `PrintReads`) — removed

**Why.** The WES template restricts `BaseRecalibrator` with `-L ${target_region}`.
Over a 71 Mb exome that leaves plenty of known sites to fit covariate tables. Over
51.9 kb it leaves a few hundred — far below what GATK needs for a stable estimate.
The outcome is either an error or, worse, a confident but meaningless
recalibration applied to every base.

Running BQSR unrestricted does not rescue it either: panel libraries carry very few
off-target reads to widen the training set with.

Skipping BQSR on small panels is the standard choice. Base qualities from the
sequencer are used as-is.

### 3. VQSR (`VariantRecalibrator` / `ApplyRecalibration`) — removed, hard filtering kept

**Why.** GATK requires at least 30 exomes or one whole genome to train
`VariantRecalibrator`. The WES template already knows this — there is a comment
admitting it, and a `--maxGaussians 1` workaround that suppresses the "no data"
crash without making the model meaningful.

One sample over 51.9 kb yields on the order of tens of variants. No amount of tuning
makes a Gaussian mixture model trainable on that. GATK's own documented fallback for
small targets is hard filtering — which the WES template already ran downstream of
VQSR anyway. We keep that and drop the pretence.

**And joint genotyping does not change this**, contrary to the natural assumption
(including one I made earlier and am correcting here). 143 samples across 51.9 kb
give on the order of hundreds to ~1500 variant *sites* — VariantRecalibrator wants
thousands. The cohort is wide, not deep: joint calling adds samples per site, not
sites. Hard filtering stays, and it lives in `run_cohort.sh`.

### 4. `-L` target restriction with padding — added `-ip 50`

Calling is restricted to `Covered.bed`. Each region is padded by 50 bp because
amplicon ends are fixed: a variant sitting on a boundary would otherwise be missed
or called with truncated context.

### 5. Coverage QC over the panel — new step

A sample that failed enrichment still produces a perfectly valid-looking GVCF, just
with fewer calls — and once it is folded into the cohort, that failure is even
harder to spot. The new step reports mean depth and the fraction of target bases at
≥20× and ≥100× per sample; `run_cohort.sh` aggregates them and flags any sample
below 90% at 20×.

Note it measures over **`Amplicons.bed`** (what was physically amplified) while
calling uses **`Covered.bed`**. Mixing the two up is the classic HaloPlex QC error —
it produces plausible numbers that mean something else.

### 6. ANNOVAR annotation — new step

Run **once on the cohort VCF**, not 143 times per sample: same result, a fraction of
the work, one file to reason about. See below.

### 6a. Genotyping moved from per-sample to cohort

The WES pipeline emits a finished VCF per sample. This one emits a **GVCF** per
sample and genotypes the whole cohort in a second stage. Reasons in "Why joint
genotyping" below.

### 7. Smaller changes

- `RGPU` was hardcoded to `SureSelectV4` — a different capture chemistry, i.e. wrong
  metadata in every BAM header. Now configurable, set to `HaloPlex_ZMF96cardio_1`.
- Heap dropped from 64 G to 8 G. 64 G was sized for exome BAMs; here it only slows
  JVM startup.
- Reference paths moved out of the Python source into `settings.json`. The WES
  generator baked `/home/PublicData/...` into the code, so relocating a reference
  file meant editing Python.
- Thread count and interval padding are configurable instead of hardcoded.

---

## Sample-name delimiter — a trap worth spelling out

Fastq files are named:

```
Haloplex_18.i7_160606.R1.fastq.gz
```

The WES default `R1_fastq_delimiter` is `"_"`. Splitting that filename on `_` gives
`Haloplex` — **for every single file**. All 284 fastqs would collapse into one
"sample" and the generator would emit a single bogus script.

The delimiter must be `"."`, which yields `Haloplex_18`. This is set correctly in
`settings.json`, and `generate_haloplex_sh.py` now aborts loudly if the configured
delimiter collapses many fastq files into one sample.

---

## ANNOVAR

**Nothing needs downloading.** A complete `humandb` already exists on the Kilo NAS:

```
/mnt/nas/Kilo/ds1821p_III/annovar_20200608/annovar_src/humandb    2.0 TB, 188 files
```

It includes `refGene`, `avsnp150`/`avsnp151`, `clinvar_20240917`, `dbnsfp35a/42a/42c/47a`,
`gnomad211_exome`, `gnomad211_genome`, `exac03`, `intervar_20180118`, `dbscsnv11`,
`revel`, `mcap`, `cosmic70`, `cadd`/`cadd13`, `fathmm`, `eigen`, `gerp++`, and more.
The ANNOVAR Perl scripts sit next to it in `annovar_src/`.

### Cost of re-downloading, if it ever comes up

Measured throughput from `lbsb-hp-z800` to the ANNOVAR download server:
**≈1.7 MB/s**.

| Scope | Size | Time at 1.7 MB/s |
|---|---|---|
| The 8 databases this pipeline uses | ~70 GB | **~12 hours** |
| Adding `gnomad211_genome` | ~140 GB | ~24 hours |
| The full 2.0 TB mirror on Kilo | 2048 GB | **~14 days** |

The three single largest files are ~350 GB each (`fathmm`, `cadd`, `cadd13`) and
~258 GB (`eigen`). Re-downloading the full mirror is not a reasonable operation on
this link — use the NAS copy.

### Protocol chosen

```
refGene,avsnp151,clinvar_20240917,dbnsfp42a,gnomad211_exome,intervar_20180118,dbscsnv11,revel
g,f,f,f,f,f,f,f
```

Consequence, rsIDs, clinical assertions, in-silico predictors, population frequency,
ACMG-style classification, splice-site effects, missense pathogenicity. CADD, FATHMM
and Eigen are available on the NAS but are ~350 GB each and add little on a 96-gene
cardio panel — leave them out unless a specific question needs them.

---

## Conda environment

```bash
# once per machine — see the note below
conda config --system --remove channels defaults
conda config --system --set channel_priority strict

conda env create -f environment.yml
conda activate haloplex
```

> **Clear `defaults` first, or the solve fails.** Recent conda refuses to touch
> `repo.anaconda.com/pkgs/{main,r}` until their commercial Terms of Service are
> accepted, dying with `CondaToSNonInteractiveError`. Accepting those on an
> institution's behalf is not a call to make in passing, and it is unnecessary —
> every dependency here is in conda-forge or bioconda.
>
> Note `conda env create` has **no** `--override-channels` / `-c` flags (those
> belong to `conda create`, a different command), so the channel list has to come
> from the config and the yml. Verified on conda 26.5.3.

### Why GATK 4

The WES pipeline this was derived from uses GATK 3 (`-T HaplotypeCaller`, `-o`).
This one uses **GATK 4**, for three reasons:

1. **GATK 3 cannot be installed unattended.** bioconda ships it as a wrapper
   *without* the jar for licensing reasons; the jar must be downloaded from the
   Broad by hand and registered with `gatk3-register`. GATK 4 installs complete
   from bioconda.
2. **`GenomicsDBImport`.** Joint genotyping across 143 samples is what the GVCFs
   are for. GATK 4 has the tool built for it; GATK 3 would need `CombineGVCFs`
   run in hand-managed batches.
3. GATK 3 needs Java 8, which is itself EOL. GATK 4 runs on Java 17.

The CLI differs accordingly: `gatk HaplotypeCaller -O out` rather than
`gatk3 -T HaplotypeCaller -o out`.

### Checklist after creating the environment

1. **Verify every tool resolves:**
   ```bash
   for t in bwa samtools bcftools picard gatk vcf-sort bgzip tabix fastqc; do
     printf '%-12s %s\n' "$t" "$(command -v $t || echo MISSING)"
   done
   gatk --version
   ```
2. **Check the reference is indexed** — bwa indices and the `.fai`/`.dict` must sit
   next to the fasta:
   ```bash
   ls PublicData/hg19/ucsc.hg19.fasta.{amb,ann,bwt,pac,sa,fai}
   ```
   If `ucsc.hg19.dict` is missing: `picard CreateSequenceDictionary R=ucsc.hg19.fasta`
3. **Check ANNOVAR reachable:**
   ```bash
   ls /mnt/nas/Kilo/ds1821p_III/annovar_20200608/annovar_src/table_annovar.pl
   ```
   The NAS mounts are systemd automounts — the first `ls` may take a few seconds.

---

## Running

The pipeline runs in **two stages**.

**Stage 1 — per sample.** Aligns and emits one GVCF per sample. No genotypes yet.

```bash
conda activate haloplex
cd haloplex
python3 generate_haloplex_sh.py -j settings.json     # one .sh per sample
bash <script_dir>/Haloplex_18.sh                     # a single sample
ls <script_dir>/*.sh | xargs -P 4 -n 1 bash          # 4 samples at a time
```

**Stage 2 — the cohort, once.** Joint genotyping, filtering, annotation.

```bash
./run_cohort.sh settings.json
```

`run_cohort.sh` refuses to run on a partial cohort: joint genotyping a subset
produces different allele frequencies and different genotypes than the full set,
so a half-finished run must not be mistaken for a finished one. Override with
`COHORT_ALLOW_PARTIAL=1` only if you mean it.

The step-token/lock design is inherited unchanged from the WES pipeline: each step
writes a `token.<sample>.<step>` file on success and is skipped on re-run, and a
first/final lock pair prevents two machines running the same sample. Stage 2 is
likewise resumable — each artefact is skipped if it already exists.

### Why joint genotyping

Per-sample calling was the WES pipeline's model; this one defers genotyping to
the cohort. That buys three things:

1. **A squared-off matrix.** Every sample gets a genotype at every site the cohort
   varies at, so "homozygous reference" stops being indistinguishable from "no
   coverage here". For a cohort VCF that distinction is the whole point.
2. **Sensitivity at low coverage.** A weak signal in one sample is judged in the
   light of confident calls at the same site across the other 142.
3. **One consistent set of thresholds** applied across the cohort, rather than 143
   independent decisions.

It does **not** make VQSR usable — see above.

---

## Open questions — decide before the production run

1. **hg19 vs hg38.** The source project directory is named
   `HaloPlex_Dauren_0025_0026_170518_hg38`, but `runs.txt` states the design is
   hg19/GRCh37, the BED files are in hg19 coordinates, and both LabBandSB pipelines
   are hg19. The folder name looks like a leftover from an earlier hg38 attempt.
   This pipeline is hg19 throughout. Confirm that is what is wanted.

2. ~~**Single-sample vs joint calling.**~~ **Decided: joint.** Implemented as
   `-ERC GVCF` per sample → `GenomicsDBImport` → `GenotypeGVCFs`. Note the
   correction to an earlier claim of mine: joint calling does *not* make VQSR
   trainable at this panel size. Hard filtering remains.

3. ~~**Sample `Haloplex_336`.**~~ **Resolved — nothing is missing.**
   `Important_ChangeLog.txt` in the Lustre archive on mng001 records a 2018-09-03
   correction: `577 -> 527` and `336 -> 369`. The cohort VCF still carries the old,
   wrong numbers, which is exactly why `369` and `527` looked like extras on disk.
   **Consequence:** joining new results back to the old VCF by sample name will
   silently fail for those two. Decide once whether to keep the corrected names
   (recommended) or map back, and record the mapping with the output.

4. **Comparison against the previous analysis.** Timestamped BAMs
   (`..._Sorted__29Aug2018_..._reHeader__04Mar2019_...bam`) sit next to the panel on
   Kilo — that naming is Agilent SureCall's. There is a prior result set to compare
   the new calls against.
