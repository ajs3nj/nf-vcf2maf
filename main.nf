#!/usr/bin/env nextflow

nextflow.enable.dsl=2


// Process for downloading files from Synapse
process SYNAPSE_GET {

  tag "${meta.synapse_id}"

  container "sagebionetworks/synapsepythonclient:v2.6.0"

  secret "SYNAPSE_AUTH_TOKEN"

  input:
  tuple val(meta), val(synapse_id)

  output:
  tuple val(meta), path('*')

  script:
  """
  synapse get ${synapse_id}

  shopt -s nullglob
  for f in *\\ *; do mv "\${f}" "\${f// /_}"; done
  """

}


// Process for decompressing and extracting the VEP cache tarball
process EXTRACT_TAR_GZ {

  container "sagebionetworks/vcf2maf:107.2"

  input:
  path vep_tarball

  output:
  path "vep_data"

  script:
  """
  mkdir -p vep_data/
  tar -zxf ${vep_tarball} -C vep_data/
  """

}


// Process for annotating VCF file and converting to MAF
process VCF2MAF {

  tag "${meta.synapse_id}"

  container "sagebionetworks/vcf2maf:107.2"

  cpus   6
  memory { 32.GB * task.attempt }

  errorStrategy = 'retry'
  maxRetries 3

  afterScript "rm -f intermediate*"

  input:
  tuple val(meta), path(input_vcf)
  tuple path(reference_fasta), path(reference_fasta_fai)
  path vep_data

  output:
  tuple val(meta), path("*.maf")

  // TODO: Remove hard-coded VEP path
  // TODO: Handle VCF genotype columns per variant caller
  script:
  vep_path  = "/root/miniconda3/envs/vep/bin"
  vep_forks = 4
  basename  = input_vcf.name.replaceAll(/.gz$/, "").replaceAll(/.vcf$/, "")
  """
  if [[ ${input_vcf} == *.gz ]]; then
    zcat ${input_vcf} > intermediate.vcf
  else
    cat  ${input_vcf} > intermediate.vcf
  fi

  vcf2maf.pl \
    --input-vcf intermediate.vcf --output-maf intermediate.maf.raw \
    --ref-fasta ${reference_fasta} --vep-data ${vep_data}/ \
    --ncbi-build ${params.ncbi_build} --max-subpop-af ${params.max_subpop_af} \
    --vep-path ${vep_path} --maf-center ${params.maf_center} \
    --tumor-id '${meta.biospecimen_id}' --vep-forks ${vep_forks} \
    --species ${params.species}

  grep -v '^#' intermediate.maf.raw > '${basename}.maf'
  """

}


// Process for filtering MAF files for passed variants
process FILTER_MAF {

  tag "${meta.synapse_id}"

  container "python:3.10.4"

  input:
  tuple val(meta), path(input_maf)

  output:
  tuple val(meta), path("*.passed.maf")

  script:
  """
  filter_maf.py ${input_maf} '${input_maf.baseName}.passed.maf'
  """

}


// Process for merging study MAF files
process MERGE_MAFS {

  tag "${meta.study_id}-${meta.variant_class}-${meta.variant_caller}"

  container "python:3.10.4"

  memory { 16.GB * task.attempt }

  errorStrategy = 'retry'
  maxRetries 3

  input:
  tuple val(meta), path(input_mafs)

  output:
  tuple val(meta), path("*.merged.maf")

  script:
  prefix = "${meta.study_id}-${meta.variant_class}-${meta.variant_caller}"
  """
  merge_mafs.py -i '${input_mafs.join(',')}' -o '${prefix}.merged.maf'
  """
}


// Process for uploading files to Synapse
process SYNAPSE_STORE {

  tag "${parent_id}/${input.name}"

  container "sagebionetworks/synapsepythonclient:v2.6.0"

  secret "SYNAPSE_AUTH_TOKEN"

  input:
  tuple path(input), val(parent_id)

  script:
  """
  synapse store --parentId ${parent_id} ${input}
  """

}


// Workflow for generating sample-level MAF files
workflow SAMPLE_MAFS {

  take:
    sample_vcfs

  main:
    // Pair up FASTA and FAI reference files
    ref_fasta_pair = [params.reference_fasta, params.reference_fasta_fai]

    // Download VCF files from Synapse
    SYNAPSE_GET(sample_vcfs)

    // Decompress and extract VEP cache tarball
    EXTRACT_TAR_GZ(params.vep_tarball)

    // Run vcf2maf on each vcf file
    VCF2MAF(SYNAPSE_GET.out, ref_fasta_pair, EXTRACT_TAR_GZ.out)

    // Upload MAF files to Synapse
    sample_mafs_ch = VCF2MAF.out
      .map { meta, maf -> [ maf, meta.sample_parent_id ] }
    SYNAPSE_STORE(sample_mafs_ch)

  emit:
    VCF2MAF.out

}


// Workflow for generating study-level MAF files
workflow STUDY_MAFS {

  take:
    sample_mafs

  main:
    // Only consider releasable MAF files
    releasable_mafs = sample_mafs
      .filter { meta, maf -> meta.is_releasable }
    
    // Filter MAF files for passed variants
    FILTER_MAF(releasable_mafs)

    // Group MAF files by study and merge
    mafs_by_study_ch = FILTER_MAF.out
      .map { meta, maf -> subset_study_meta(meta, maf) }
      .groupTuple( by: 0 )
    MERGE_MAFS(mafs_by_study_ch)

    // Upload study MAF files to Synapse
    merged_mafs_ch = MERGE_MAFS.out
      .map { meta, maf -> [ maf, meta.merged_parent_id ] }
    SYNAPSE_STORE(merged_mafs_ch)

}


// Entrypoint workflow
workflow {

  // Parse input CSV sample sheet
  input_vcfs_ch = Channel
    .fromPath ( params.input )
    .splitCsv ( header:true, strip:true )
    .map { create_vcf_channel(it) }

  // Process individual sample VCF files
  SAMPLE_MAFS(input_vcfs_ch)

  // Filter and merge MAF files by study
  STUDY_MAFS(SAMPLE_MAFS.out)

}


// Function to get list of [ meta, vcf ]
def create_vcf_channel(LinkedHashMap row) {

  // Create metadata element
  def meta = [:]
  meta.synapse_id       = row.synapse_id
  meta.biospecimen_id   = row.biospecimen_id
  meta.sample_parent_id = row.sample_parent_id
  meta.merged_parent_id = row.merged_parent_id
  meta.study_id         = row.study_id
  meta.variant_class    = row.variant_class
  meta.variant_caller   = row.variant_caller
  meta.is_releasable    = row.is_releasable.toBoolean()

  // Combine with VCF file element
  def vcf_meta = [meta, row.synapse_id]

  return vcf_meta
}


// Function to get list of [ study_meta, maf ]
def subset_study_meta(vcf_meta, maf) {

  // Subset metadata element
  def study_meta = [:]
  study_meta.merged_parent_id = vcf_meta.merged_parent_id
  study_meta.study_id         = vcf_meta.study_id
  study_meta.variant_class    = vcf_meta.variant_class
  study_meta.variant_caller   = vcf_meta.variant_caller

  return [ study_meta, maf ]
}
