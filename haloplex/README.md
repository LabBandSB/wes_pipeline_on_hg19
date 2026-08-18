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

BED files (exported from SureDesign, found on the Kilo NAS):

| File | Regions | Total bp | Used for |
|---|---|---|---|
| `ZMF96cardio_1_Covered.bed` | 60 | **51 895** (51.9 kb) | variant calling (`-L`) |
| `ZMF96cardio_1_Regions.bed` | 60 | 55 515 | reference |
| `ZMF96cardio_1_Amplicons.bed` | 1202 | 240 028 | coverage QC |
| `ZMF96cardio_1_AllTracks.bed` | 1430 | — | reference |

Chromosomes touched: chr1, 2, 3, 4, 7, 10, 12, 15, 17, 18, 19.

**51.9 kb is the number that drives almost every decision below.** It is about
1/1400th of an exome. Methods that are routine at exome scale run out of data here.

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

### 4. `-L` target restriction with padding — added `-ip 50`

Calling is restricted to `Covered.bed`. Each region is padded by 50 bp because
amplicon ends are fixed: a variant sitting on a boundary would otherwise be missed
or called with truncated context.

### 5. Coverage QC over the panel — new step

A sample that failed enrichment still produces a perfectly valid-looking VCF, just
with fewer calls. On a panel this small that is easy to miss. The new step reports
mean depth and the fraction of target bases at ≥20× and ≥100× per sample.

Note it measures over **`Amplicons.bed`** (what was physically amplified) while
calling uses **`Covered.bed`**. Mixing the two up is the classic HaloPlex QC error —
it produces plausible numbers that mean something else.

### 6. ANNOVAR annotation — new step

See below.

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
conda env create -f environment.yml
conda activate haloplex
```

### Checklist after creating the environment

1. **Java must be 8.** GATK 3 fails on newer JDKs with unhelpful errors.
   ```bash
   java -version        # expect "1.8.0_..."
   ```
2. **Register the GATK 3 jar.** bioconda ships a wrapper *without* the jar for
   licensing reasons. Download `GenomeAnalysisTK-3.8-1` from the Broad, then:
   ```bash
   gatk3-register /path/to/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
   gatk3 --version      # must print 3.8-1-0-gf15c1c3ef
   ```
   Without this the pipeline fails at the first HaplotypeCaller step.
3. **Verify every tool resolves:**
   ```bash
   for t in bwa samtools bcftools picard gatk3 vcf-concat vcf-sort bgzip tabix fastqc; do
     printf '%-12s %s\n' "$t" "$(command -v $t || echo MISSING)"
   done
   ```
4. **Check the reference is indexed** — bwa indices and the `.fai`/`.dict` must sit
   next to the fasta:
   ```bash
   ls PublicData/hg19/ucsc.hg19.fasta.{amb,ann,bwt,pac,sa,fai}
   ```
   If `ucsc.hg19.dict` is missing: `picard CreateSequenceDictionary R=ucsc.hg19.fasta`
5. **Check ANNOVAR reachable:**
   ```bash
   ls /mnt/nas/Kilo/ds1821p_III/annovar_20200608/annovar_src/table_annovar.pl
   ```
   The NAS mounts are systemd automounts — the first `ls` may take a few seconds.

---

## Running

```bash
cd haloplex
python3 generate_haloplex_sh.py -j settings.json     # writes one .sh per sample
bash <script_dir>/Haloplex_18.sh                     # single sample
ls <script_dir>/*.sh | xargs -P 4 -n 1 bash          # 4 samples in parallel
```

The step-token/lock design is inherited unchanged from the WES pipeline: each step
writes a `token.<sample>.<step>` file on success and is skipped on re-run, and a
first/final lock pair prevents two machines running the same sample.

---

## Open questions — decide before the production run

1. **hg19 vs hg38.** The source project directory is named
   `HaloPlex_Dauren_0025_0026_170518_hg38`, but `runs.txt` states the design is
   hg19/GRCh37, the BED files are in hg19 coordinates, and both LabBandSB pipelines
   are hg19. The folder name looks like a leftover from an earlier hg38 attempt.
   This pipeline is hg19 throughout. Confirm that is what is wanted.

2. **Single-sample vs joint calling.** Both LabBandSB pipelines are explicitly
   "single sample". With 143 samples, joint genotyping (`-ERC GVCF` →
   `CombineGVCFs` → `GenotypeGVCFs`) would give better sensitivity at low coverage
   and consistent reference/no-call distinction across the cohort — which matters
   for a cohort VCF. It would also make VQSR trainable again. This is an
   architectural change and is deliberately *not* in this version.

3. **Sample `Haloplex_336`.** Present in the cohort VCF header, no raw data found on
   any NAS. `Haloplex_39` and `Haloplex_577` were also missing from the main project
   directory but were found in the Kilo archive under different names
   (`Sample39_R[12].fastq.gz`, `577_S1_R[12]_001.fastq.gz`) and copied into
   `workflow/haloplex/from_kilo_archive/`. `336` may likewise be a renamed file.

4. **Comparison against the previous analysis.** Timestamped BAMs
   (`..._Sorted__29Aug2018_..._reHeader__04Mar2019_...bam`) sit next to the panel on
   Kilo — that naming is Agilent SureCall's. There is a prior result set to compare
   the new calls against.
