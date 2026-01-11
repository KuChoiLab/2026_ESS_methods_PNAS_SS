# 2026_ESS_methods_PNAS_SS
## [WES]
## Variant calling
**1. Detection of somatic SNVs and INDELs**

The somatic mutation calling pipeline was adapted based on sample type. For matched tumor-normal pairs, single-nucleotide variants (SNVs) were called with Mutect2, while small insertions and deletions (indels) were identified using a consensus of both MuTect2 and Strelka2.

For unmatched tumor samples, somatic variants were called using MuTect2 in tumor-only mode.

The following code is used to run GATK Mutect2.

```
# GATK version : 4.6.0.0

# Running GATK Mutect2
# Example of fasta file : hs38DH.fasta
# Example of a panel of normals (PON) file: 1000g_pon.hg38.vcf.gz

gatk Mutect2 \
        -R ${path_to_fasta_file} \
        -I ${path_to_tumor_bam} \
        -I ${path_to_normal_bam} \
        -tumor ${tumor_sample_name} \
        -normal ${normal_sample_name} \
        --panel-of-normals ${path_to_PON_file} \
        --f1r2-tar-gz ${path_to_output_directory}/${sample_pair_name}.f1r2.tar.gz \
        -O ${path_to_output_directory}/${sample_pair_name}.mutect2.vcf \
        -bamout ${path_to_output_directory}/${sample_pair_name}.mutect2.bam

# Running GATK LearnReadOrientationModel
gatk LearnReadOrientationModel \
	 -I ${path_to_output_directory}/${sample_pair_name}.f1r2.tar.gz \
	 -O ${path_to_output_directory}/${sample_pair_name}.read-orientation-model.tar.gz

# Contamination calculation and filtering
gatk GetPileupSummaries \
    -R ${path_to_fasta_file} \
    -I ${path_to_tumor_bam} \
    -O ${path_to_output_directory}/${tumor_sample_name}_getpileupsummaries.table

gatk GetPileupSummaries \
    -R ${path_to_fasta_file} \
    -I ${path_to_normal_bam} \
    -O ${path_to_output_directory}/${normal_sample_name}_getpileupsummaries.table

gatk CalculateContamination \
    -I ${path_to_output_directory}/${tumor_sample_name}_getpileupsummaries.table \
    -matched ${path_to_output_directory}/${normal_sample_name}_getpileupsummaries.table \
    -O ${path_to_output_directory}/${sample_pair_name}_calculatecontamination.table \
    -tumor-segmentation ${path_to_output_directory}/${sample_pair_name}_segments.table

gatk FilterMutectCalls \
    -V ${path_to_output_directory}/${sample_pair_name}.mutect2.vcf.gz \
    -R ${path_to_fasta_file} \
    --ob-priors ${path_to_output_directory}/${sample_pair_name}.read-orientation-model.tar.gz \
    --contamination-table ${path_to_output_directory}/${sample_pair_name}_calculatecontamination.table \
    --tumor-segmentation ${path_to_output_directory}/${sample_pair_name}_segments.table \
    --stats ${path_to_output_directory}/${sample_pair_name}.mutect2.vcf.gz.stats \
    -O ${path_to_output_directory}/${sample_pair_name}.mutect2.filtered.vcf
```

**2. Detection of germline SNVs and INDELs**

The following code is used to run GATK HaplotypeCaller.

```
gatk HaplotypeCaller \
   -R ${path_to_fasta_file} \
   -I ${path_to_normal_bam} \
   -O ${path_to_output_directory}/${sample_name}.g.vcf.gz \
   -ERC GVCF
```

**3. Filtering FFPE-induced artifacts**

The following code is used to run FFPolish.

```
# FFPolish version : 0.1.0
# ffpolish is available upon installation of FFPolish

# Filter FFPE artifacts from VCF
ffpolish filter \
  -o ${output_directory} \
  -p ${output_prefix} \
  ${reference_fasta} \
  ${vcf_gz} \
  ${ffpe_tumor_bam}
```

**4. Detection of CNVs**

The CNV calling pipeline was adapted based on the sample type. For matched tumor–normal pairs, copy number variations were inferred using both FACETS and Sequenza. 

For unmatched tumor samples, CNVs were identified using CNVkit.

The following code is used to run FACETS.
   
```
# FACETS version : 0.6.2
# snp-pileup and cnv_facets.R are available upon installation of the FACETS.

facets/inst/extcode/snp-pileup \
${path_to_output_directory}/${sample_pair_name}.csv.gz \
${path_to_normal_bam} \
${path_to_tumor_bam}

Rscript facets/bin/cnv_facets.R -p ${path_to_output_directory}/${sample_pair_name}.csv.gz -o ${path_to_output_directory}/${output_file}

```
The following code is used to run sequenza.

```
# Sequenza version: 3.9.0
# sequenza-utils is available upon installation of sequenza

sequenza-utils bam2seqz \
  -gc /path/to/hg38.gc50Base.wig.gz \
  --fasta /path/to/reference/hg38.fasta \
  -n ${path_to_normal_bam} \
  -t ${path_to_tumor_bam} \
  -o ${path_to_output_directory}/${sample_pair_name}.seqz.gz

sequenza-utils seqz_binning \
  --seqz ${path_to_output_directory}/${sample_pair_name}.seqz.gz \
  -w 50 \
  -o ${path_to_output_directory}/${sample_pair_name}.binned.seqz.gz
```

The following code is used to run CNVkit.

```
# CNVkit version: 0.9.10
# cnvkit.py is available upon installation of CNVkit
# Tumor-only mode (without matched normal)
cnvkit.py batch ${tumor_bam} \
  -n \
  --targets ${baits_bed} \
  --fasta ${reference_fasta} \
  --access ${access_bed} \
  --output-reference ${output_reference_cnn} \
  --output-dir ${output_directory}
```

## Downstream analysis

**1. Identifying mutational signatures**

The following code is used to run deconstructSigs.

```
# deconstructSigs version: 
# Convert mutation data to deconstructSigs input format
sigs.input <- mut.to.sigs.input(
  mut.ref = ${mutation_dataframe},
  sample.id = "Sample",
  chr = "chr",
  pos = "pos",
  ref = "ref",
  alt = "alt"
)

# Determine signature contributions for a sample
output <- whichSignatures(
  tumor.ref = sigs.input,
  signatures.ref = signatures.nature2013,
  sample.id = ${sample_id},
  contexts.needed = TRUE,
  tri.counts.method = 'default'
)

# Visualize results
plotSignatures(output, sub = ${sample_name})
makePie(output, sub = ${sample_name})
```

**2. Identifying driver genes**

The following code is used to run MutPanning.

```
# MutPanning version: 
java -Xmx8G -classpath ${mutpanning_dir}/commons-math3-3.6.1.jar:${mutpanning_dir}/jdistlib-0.4.5-bin.jar:${mutpanning_dir} MutPanning \
  "${root_directory}" \
  "${maf_file}" \
  "${sample_annotation_file}" \
  "${hg19_directory}"
```

**3. Identifying recurrent CNV**

The following code is used to run GISTIC2.

```
# GISTIC2 version: 8.3
gistic2 \
  -b ${output_directory} \
  -seg ${segmentation_file} \
  -refgene ${refgene_mat_file} \
  -genegistic  \
  -smallmem  \
  -broad  \
  -brlen  \
  -conf  \
  -armpeel  \
  -savegene  \
  -gcm extreme
```

## [RNA-seq]

## Identifying fusions

**1. Fusion calling**

Fusion calling was performed using both STAR-Fusion and Arriba to identify gene fusion events from RNA-seq data.

The following code is used to run STAR-Fusion.

```
# STAR-Fusion version: 1.15.0
STAR-Fusion \
  --genome_lib_dir ${ctat_genome_lib_dir} \
  --left_fq ${read1_fastq} \
  --right_fq ${read2_fastq} \
  --output_dir ${output_directory} \
  --CPU ${threads}
```

The following code is used to run Arriba.

```
# Arriba version: 2.5.0
arriba \
  -x /dev/stdin \
  -o ${output_fusions_tsv} \
  -O ${output_discarded_tsv} \
  -a ${assembly_fasta} \
  -g ${annotation_gtf} \
  -b ${blacklist_tsv} \
  -k ${known_fusions_tsv} \
  -t ${known_fusions_tsv} \
  -p ${protein_domains_gff3}

# Visualize fusions 
Rscript draw_fusions.R \
  --fusions=${output_fusions_tsv} \
  --alignments=${aligned_bam} \
  --output=${output_pdf} \
  --annotation=${annotation_gtf} \
  --cytobands=${cytobands_tsv} \
  --proteinDomains=${protein_domains_gff3}
```

## Transcriptomic profiling of ESS

**1. Differential Expression analysis**

The following code is used to run DESeq2.

```
# DESeq2 version: 1.38.3
# Create DESeq2 dataset from count matrix
dds <- DESeqDataSetFromMatrix(
  countData = ${count_matrix},
  colData = ${sample_metadata},
  design = ~ ${design_formula}
)

# Run differential expression analysis
dds <- DESeq(dds)

# Extract results
res <- results(dds, name = ${coefficient_name})

# Or with log fold change shrinkage (recommended)
res <- lfcShrink(dds, coef = ${coefficient_name}, type = "apeglm")

```

**2. Gene set enrichment analysis**

The following code is used to run fgsea.

```
# fgsea version: 1.38.3
# Run fgsea with ranked gene list
fgseaRes <- fgsea(
  pathways = ${pathway_list},
  stats = ${ranked_stats},
  minSize = 15,
  maxSize = 500
)
```

## [WGS]

**1. Identifying Complex Structural Variations**

The following code is used to run JaBbA.

```
# JaBbA version: 1.0
jba ${junctions_file} ${coverage_file} \
  --seg ${segmentation_rds} \
  --blacklist.junctions ${blacklist_junctions} \
  --field ${coverage_field} \
  --ploidy ${ploidy} \
  --purity ${purity} \
  --outdir ${output_directory} \
  --name ${sample_name}
```

## Reference

1. McKenna A, Hanna M, Banks E, Sivachenko A, Cibulskis K, Kernytsky A, Garimella K, Altshuler D, Gabriel S, Daly M, DePristo MA. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res. 2010;20(9):1297–1303.
2. Kim S, Scheffler K, Halpern AL, et al. Strelka2: fast and accurate calling of germline and somatic variants. Nat Methods. 2018;15(8):591–594.
3. Shen R, Seshan VE. FACETS: allele-specific copy number and clonal heterogeneity analysis tool for high-throughput DNA sequencing. Nucleic Acids Res. 2016;44(16):e131.
4. Favero F, Joshi T, Marquard AM, et al. Sequenza: allele-specific copy number and mutation profiles from tumor sequencing data. Ann Oncol. 2015;26(1):64-70.
5. Talevich E, Shain AH, Botton T, Bastian BC. CNVkit: Genome-Wide Copy Number Detection and Visualization from Targeted DNA Sequencing. PLoS Comput Biol. 2016;12(4):e1004873.
6. Rosenthal R, McGranahan N, Herrero J, Taylor BS, Swanton C. DeconstructSigs: delineating mutational processes in single tumors distinguishes DNA repair deficiencies and patterns of carcinoma evolution. Genome Biol. 2016;17:31.
7. Dietlein F, Weghorn D, Taylor-Weiner A, et al. Identification of cancer driver genes based on nucleotide context. Nat Genet. 2020;52(2):208-218.
8. Mermel CH, Schumacher SE, Hill B, Meyerson ML, Beroukhim R, Getz G. GISTIC2.0 facilitates sensitive and confident localization of the targets of focal somatic copy-number alteration in human cancers. Genome Biol. 2011;12(4):R41.
9. Haas BJ, Dobin A, Li B, Stransky N, Pochet N, Regev A. Accuracy assessment of fusion transcript detection via read-mapping and de novo fusion transcript assembly-based methods. Genome Biol. 2019;20(1):213.
10. Uhrig S, Ellermann J, Walther T, et al. Accurate and efficient detection of gene fusions from RNA sequencing data. Genome Res. 2021;31(3):448-460.
11. Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 2014;15(12):550.
12. Korotkevich G, Sukhov V, Sergushichev A. Fast gene set enrichment analysis. bioRxiv. Published online June 20, 2016.
13. Hadi K, Yao X, Behr JM, et al. Distinct Classes of Complex Structural Variation Uncovered across Thousands of Cancer Genome Graphs. Cell. 2020;183(1):197-210.e32.

