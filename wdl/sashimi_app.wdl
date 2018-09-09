task sashimi_app {

  # inputs for preprocessing
  File quantsf
  Array[String]? chromosomes
  Array[File]? blacklist
  # inputs to call_variants
  String smoothing_strategy = "medianfilter"
  Int smoothing_window = 3
  # inputs for postprocessing
  Int? merge_distance = 0
  Int? min_variant_len = 500

  command {
    python /opt/sashimi.py ${quantsf} \
        --chromosomes $(echo ${sep=' ' chromosomes}) \
        --blacklist $(echo ${sep=' ' blacklist}) \
        --smoothing_strategy ${smoothing_strategy} \
        --smoothing_window ${smoothing_window} \
        --merge_distance ${merge_distance} \
        --min_variant_len ${min_variant_len}
  }

  output {
    # Raw data in BED, sorted
    File quant_sorted = "quant_sorted.bed" 

    # Files storing a caller output (1 or 0) for each region
    File all_events_del_ho = "all_events_del_ho.bed"
    File all_events_del_he = "all_events_del_he.bed"

    # Final postprocessed BED with identified deletions
    File merged_events_del_ho = "merged_events_del_ho.bed"
    File merged_events_del_he = "merged_events_del_he.bed"
  }

  runtime {
    docker: "commandlinegirl/sashimi:0.3"
  }
}

