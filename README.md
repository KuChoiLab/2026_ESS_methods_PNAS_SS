# 2026_ESS_methods_PNAS_SS
## WES
### Variant calling
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
# FFPolish version : 
# ffpolish is available upon installation of FFPolish

# Filter FFPE artifacts from VCF
ffpolish filter \
  -o ${output_directory} \
  -p ${output_prefix} \
  ${reference_fasta} \
  ${vcf_gz} \
  ${ffpe_tumor_bam}
```

**4. Detection of somatic CNVs**

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
