
#########################################
# Preprocess quant.sf - the output of Salmon
#########################################

task sashimi_unified_task {
  File quantsf
  File? blacklist
  Int? merge_distance
  Int? min_variant_len
  Float? min_score
  Array[String]? chromosomes
  Array[String] smoothing_strategies = ["medianfilter"]
  Array[Int] smoothing_windows = [3]

  command {
    python /opt/sashimi_unified.py ${quantsf} \
        --merge_distance ${merge_distance} \
        --min_variant_len ${min_variant_len} \
        --min_score ${min_score} \
        --chromosomes $(echo ${sep=' ' chromosomes}) \
        --smoothing_strategies $(echo ${sep=' ' smoothing_strategies}) \
        --smoothing_windows $(echo ${sep=' ' smoothing_windows})
  }

  output {
    # Raw data in BED, sorted
    File salmon_all_sorted = "salmon_all_sorted.bed" 

    # Files storing the output of each caller in adjacent columns
    File all_classifier_outputs_ho = "all_classifier_outputs_ho.bed"
    File all_classifier_outputs_he = "all_classifier_outputs_he.bed"

    # Files storing a caller output (1 or 0) for each region
    # and a score/confidence attached to each call
    File integrated_output_ho = "integrated_output_ho.bed"
    File integrated_output_he = "integrated_output_he.bed"

    # Final postprocessed BED with identified deletions
    File result_dels_ho = "result_dels_ho.bed"
    File result_dels_he = "result_dels_he.bed"
  }

  runtime {
    docker: "commandlinegirl/sashimi"
  }
}

#########################################
# Assess the accuracy of the caller
#########################################

task evaluate_output {
  File sample_events
  File? truth_events
  Float? overlap = 0.5

  command {
    python /opt/evaluate.py \
      --sample_events ${sample_events} \
      --truth_events ${truth_events} \
      --overlap ${overlap}
  }

  output {
    File true_positives = "true_positives.bed"
    File false_positives = "false_positives.bed"
    File false_negatives = "false_negatives.bed"
    Float precision = read_float("score_precision.txt")
    Float recall = read_float("score_recall.txt")
    Float f_score = read_float("score_f_score.txt")
  }

  runtime {
    docker: "commandlinegirl/sashimi"
  }
}

#########################################
# Sashimi workflow
#########################################

workflow sashimi_simple {

  # general inputs
  Boolean analyse_hom = true
  Boolean analyse_het = false
  # inputs to read_salmon_output
  File quantsf
  Array[String]? chromosomes
  File? blacklist
  # inputs to call_variants scatter & gather
  Array[String] smoothing_strategies = ["medianfilter"]
  Array[Int] smoothing_windows = [3]
  # inputs to integrate
  Int? merge_distance = 500
  Int? min_variant_len = 1000
  Float? min_score = 1.0
  # inputs to evaluate
  Boolean evaluate = false
  File? truth_events_ho
  File? truth_events_he
  Float? overlap


  # read in data, preprocess, and normalize
  call sashimi_unified_task {
    input: 
      quantsf = quantsf,
      chromosomes = chromosomes,
      blacklist = blacklist,
      smoothing_strategies = smoothing_strategies,
      smoothing_windows = smoothing_windows,
      merge_distance = merge_distance,
      min_variant_len = min_variant_len,
      min_score = min_score
  }

  if (evaluate && analyse_hom) {
    call evaluate_output as evaluate_hom {
      input:
        sample_events = sashimi_unified_task.result_dels_ho,
        truth_events = truth_events_ho,
        overlap = overlap
    }
  }

  if (evaluate && analyse_het) {
    call evaluate_output as evaluate_het {
      input:
        sample_events = sashimi_unified_task.result_dels_he,
        truth_events = truth_events_he,
        overlap = overlap
    }
  }

}

