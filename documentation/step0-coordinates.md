# Project 1 - Step 0 Anchor: HBG1/HBG2 promoter coordinates + edit positions
# STATUS: VERIFIED GO (June 18 2026). NO-GO as-originally-drafted; corrections below applied.
# Arithmetic spine confirmed: both RefSeq TSS anchors correct; ClinVar-matched where a record exists;
# 4,924 bp homology invariant validated on GRCh38 and GRCh37.

## Build & strand convention (everything downstream trusts this)
- GRCh38: accession NC_000011.10, 1-based. GRCh37: accession NC_000011.9, 1-based.
- HBG1 & HBG2 are MINUS (reverse/complement) strand → promoters at numerically LARGER coords than CDS.
- Literature promoter positions (-198, -175, -113...) are numbered vs the experimentally defined
  transcription CAP SITE: first transcribed base = +1, immediately upstream = -1, NO zero.
  Conversion (reverse-strand genes): absolute genomic pos = TSS + |promoter offset|.
- DUAL ALLELE ORIENTATION IS MANDATORY. Literature uses promoter/SENSE-strand alleles
  ("-198 T>C"); genome browsers/VCFs use GENOMIC PLUS-strand alleles ("-198 A>G", the reverse
  complement). A single `ref_base` field is UNSAFE and guarantees a strand error. Store both.
- Canonical IDs = GENOMIC HGVS (e.g. NC_000011.10:g.5250055A>G). Current transcripts =
  NM_000559.3 (HBG1), NM_000184.3 (HBG2). Retain legacy .2 expressions ONLY as aliases.
  Do NOT use ClinVar-invalidated transcript expressions (e.g. c.-53-198T>C) as a coordinate;
  transcript HGVS counts from ATG incl. 5'UTR and is NOT promoter-offset arithmetic.

## REQUIRED SCHEMA (one row per gene per position)
genome_build | chrom | pos_1based | genomic_plus_ref | genomic_plus_alt | gene | gene_strand |
promoter_offset | promoter_sense_ref | promoter_sense_alt | editor_class | variant_db_status | mechanism_provenance
- Allele-orientation key:  sense T>C == genomic A>G ;  sense A>G == genomic T>C ;
  sense G>A == genomic C>T ;  sense C>T == genomic G>A.

## Gene anchors (RefSeq TSS: VERIFIED CORRECT. Use these, NOT Ensembl gene spans.)
- HBG1 (Gene 3047): NM_000559.3 / NG_000007.3. GRCh38 span chr11:5,248,269-5,249,857;
  TSS (first transcribed base, reverse strand) GRCh38 5,249,857 ; GRCh37 5,271,087.
- HBG2 (Gene 3048): NM_000184.3 / NG_000007.3. GRCh38 span chr11:5,253,188-5,254,781;
  TSS GRCh38 5,254,781 ; GRCh37 5,276,011.
- NCBI annotation release RS_2025_08 (records updated June 3 2026); Ensembl release 116 (June 9 2026).
  No GRCh38 TSS disagreement for the RefSeq/Ensembl-linked transcripts. Do NOT substitute broader
  Ensembl gene spans or alternate/legacy transcripts for the selected-transcript TSS.

## Proximal-promoter modeling window: -210 to -80 vs cap
- HBG1: GRCh38 chr11:5,249,937-5,250,067 ; HBG2: GRCh38 chr11:5,254,861-5,254,991
- Promoters are highly homologous → most guides co-target BOTH. Project every edit to both loci.

## IN-SCOPE base-editable positions  [genomic+ ref>alt | sense ref>alt | editor]
ABE class:
- -198 : HBG1 5,250,055 A>G (sense T>C) | HBG2 5,254,979 A>G (sense T>C) | de novo KLF1 motif
- -175 : HBG1 5,250,032 A>G (sense T>C) | HBG2 5,254,956 A>G (sense T>C) | TAL1 E-box (looping)
- -124 : HBG1 5,249,981 A>G (sense T>C) | HBG2 5,254,905 A>G (sense T>C) | synthetic; KLF1 (see note)
- -123 : HBG1 5,249,980 A>G (sense T>C) | HBG2 5,254,904 A>G (sense T>C) | synthetic; KLF1 (see note)
- -113 : HBG1 5,249,970 T>C (sense A>G) | HBG2 5,254,894 T>C (sense A>G) | de novo GATA1 motif
CBE class:
- -117 : HBG1 5,249,974 C>T (sense G>A) | HBG2 5,254,898 C>T (sense G>A) | distal BCL11A motif
- -114C>T : HBG1 5,249,971 G>A (sense C>T) | HBG2 5,254,895 G>A (sense C>T) | distal BCL11A motif
# ^ CORRECTED: HBG2 -114 = 5,254,895 (was tentatively 5,254,894, which is HBG2 -113). Fixed.

## variant_db_status  (evidentiary-tier labels: DO NOT overclaim)
- DB-CONFIRMED (ClinVar accession + rsID; gene-specific, do NOT cross-assign rsIDs between HBG1/HBG2):
    HBG1 -198  rs35710727  (RCV001814966 / VCV000015031)
    HBG1 -117  rs35378915  (RCV000016173 / VCV000015030)
    HBG1 -114  rs281860601 (RCV001814969 or RCV000016179 / VCV000015035)
    HBG2 -175  rs63750654  (RCV001814961 / VCV000014983)   <- HBG2-specific; NEVER annotate HBG1 -175
    HBG2 -114  rs34809449  (RCV000016130 / VCV000014990)   <- multiallelic; store allele-specific G>A, not just rsID
- variant_db_UNVERIFIED (experimentally demonstrated in CD34+ HSPCs, RefSeq-derived coordinate,
  but NO HBG1-specific ClinVar/dbSNP record found; usable as a target, NOT as a catalogued variant):
    HBG1 -175  (do NOT borrow HBG2 rs63750654)
    HBG1 -113  (no HBG1-specific record)
- SYNTHETIC SCREEN HITS (not clinical HPFH alleles): -123 and -124. KLF1 binding was shown for the
  -123/-124 DOUBLE substitution, not either single position alone (Ravi 2022).

## EXCLUDED (not standard ABE/CBE; keep as explicit negative cases)
- -110 A>C (sense): genomic T>G, NF-Y motif but a TRANSVERSION.  HBG1 ~5,249,967 / HBG2 ~5,254,891
- -114 C>A and -114 C>G: transversions.
- Delta-13 (~-114 to -102): a 13-nt deletion (Cas9-nuclease-modeled); base editors don't make it.

## PRIMARY ENGINEERING ANCHORS (build/demo the v0 on these first)
- Strongest EXACT-allele HSPC/HSC base-editing evidence: -175 (HBG2 DB-confirmed; HBG1 unverified),
  -113 (ABE8e, exact -113 A>G installed in patient CD34+ HSCs + in-vivo human-locus model), and
  the -123/-124 synthetic pair. Build/demo on these.
- -117 and -114: HSPC evidence is TILED / MOTIF-LEVEL (region included in base-editing screens /
  BCL11A-motif disruption), NOT a verified isolated exact-natural-allele HSC experiment. Label as such.

## BUILD-TIME INVARIANT (hard gate; failure must STOP the build)
- For equivalent HBG1/HBG2 promoter offsets:  hbg2_pos_grch38 - hbg1_pos_grch38 == 4924
  (TSS separation 5,254,781 - 5,249,857 = 4,924 on GRCh38; same on GRCh37.)
- This invariant catches row-shifts (it caught the -114/-113 off-by-one) but does NOT establish gene
  identity. NEVER infer gene by subtracting 4,924. Require absolute build-specific coord + explicit
  gene. For read processing: require a uniquely-assignable gene-specific flank, else label HBG1/HBG2_AMBIGUOUS.

## RESOLVED CITATIONS (real PMIDs; replace the earlier unverifiable tokens)
# Pattern: clinical report = variant's EXISTENCE; modern paper = MECHANISM (do not conflate; mechanisms postdate the clinical reports).
- -198 : clinical Tate/Wood/Weatherall 1986 PMID 2430647 | KLF1 mechanism Wienert 2017 PMID 28659276
         (+ competition model Doerfler 2021 PMID 34341563) | HSPC ABE Mayuranathan 2023 PMID 37400614
- -175 : clinical (HBG2) Surrey 1988 PMID 2449926 ; (HBG1, in cis w/ HBG2 -158) Coleman 1993 PMID 7679879
         | TAL1 E-box mechanism Wienert 2015 PMID 25971621 | HSPC ABE Mayuranathan 2023 PMID 37400614
- -123/-124 : Ravi 2022 PMID 35147495 (synthetic screen hit; KLF1 = double-mutant probe; HSPC-derived erythroid)
- -117 : clinical Waber 1986 PMID 2417646 | distal-CCAAT binding Superti-Furga 1988 PMID 3181130
         (+ Doerfler 2021 PMID 34341563) | region tiled Ravi 2022 PMID 35147495 ; motif BE Han 2023 PMID 37989316
- -114 : clinical (HBG2) Fucharoen 1990 PMID 1698280 ; (HBG1) Oner 1991 PMID 1704803
         | modern framework Doerfler 2021 PMID 34341563 | regional/motif HSPC BE Ravi 2022 / Han 2023
- -113 : clinical Amato 2014 PMID 23621512 | de novo GATA1 mechanism Martyn 2019 PMID 30617196
         (+ Doerfler 2021 PMID 34341563) | exact -113 A>G ABE8e Li 2022 PMID 36006707 ; also Mayuranathan 2023 PMID 37400614

## ===== APPLIED CORRECTIONS (from June 18 2026 verification) =====
## 1. HBG2 -114 set to 5,254,895 (GRCh37 5,276,125); 5,254,894 reserved for HBG2 -113.
## 2. Dual allele orientation stored (genomic_plus_* AND promoter_sense_*).
## 3. Transcripts updated to NM_000559.3 / NM_000184.3; .2 kept as legacy alias only.
## 4. Genomic HGVS canonical; ClinVar-invalidated c.-53-198T>C NOT used as a coordinate.
## 5. rs63750654 is HBG2 -175 only; never on HBG1 -175.
## 6. HBG1 -175 and HBG1 -113 marked variant_db_UNVERIFIED (experimentally demonstrated, not DB-confirmed).
## 7. -123/-124 marked SYNTHETIC screen hits; KLF1 attached to the double substitution.
## 8. -117/-114 HSPC evidence marked TILED/MOTIF-LEVEL.
## 9. Gene assignment explicit; ambiguous reads stay HBG1/HBG2_AMBIGUOUS (not assigned by guide name).
## 10. 4,924 bp invariant added as a hard build-stop gate.
