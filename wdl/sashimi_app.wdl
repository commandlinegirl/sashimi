task sashimi {

  # inputs for preprocessing
  File quantsf
  Array[String]? chromosomes
  Array[File]? blacklist
  # inputs to call_variants
  String smoothing_strategy = "medianfilter"
  Int smoothing_window = 3
  Array[Float]? del_hom_tpm_range = [0.0, 0.01]
  Array[Float]? del_het_tpm_range
  Array[Float]? dup_tpm_range
  # inputs for postprocessing
  Int? merge_distance = 0
  Int? min_variant_len = 500

  command {
    python /opt/sashimi.py ${quantsf} \
        --chromosomes $(echo ${sep=' ' chromosomes}) \
        --blacklist $(echo ${sep=' ' blacklist}) \
        --smoothing_strategy ${smoothing_strategy} \
        --smoothing_window ${smoothing_window} \
        --del_hom_tpm_range $(echo ${sep=' ' del_hom_tpm_range}) \
        --del_het_tpm_range $(echo ${sep=' ' del_het_tpm_range}) \
        --dup_tpm_range $(echo ${sep=' ' dup_tpm_range}) \
        --merge_distance ${merge_distance} \
        --min_variant_len ${min_variant_len}
  }

  output {
    # Raw data in BED, sorted
    File quant_sorted = "quant_sorted.bed" 

    # Files storing a caller output (1 or 0) for each region
    File marked_regions_del_ho = "marked_regions_del_ho.bed"
    File marked_regions_del_he = "marked_regions_del_he.bed"
    File marked_regions_dup = "marked_regions_dup.bed"

    # Final postprocessed BED with identified deletions and
    # duplications (merged)
    File merged_del_ho = "del_ho.bed"
    File merged_del_he = "del_he.bed"
    File merged_dup = "dup.bed"
  }

  runtime {
    docker: "commandlinegirl/sashimi:0.3"
  }
}

