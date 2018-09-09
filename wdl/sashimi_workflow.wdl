
#########################################
# Preprocess quant.sf - the output of Salmon
#########################################

task read_salmon_output {
  File quantsf
  Array[String]? chromosomes

  command {
    python /opt/read_salmon_output.py ${quantsf} \
        --chromosomes $(echo ${sep=' ' chromosomes})
  }

  output {
    File salmon_all_sorted = "salmon_all_sorted.bed" 
  }

  runtime {
    docker: "commandlinegirl/sashimi:0.1"
  }
}

#########################################
# Call variants
#########################################

task call_variants {
  File salmon_all_sorted
  String? smoothing_strategy
  Int? smoothing_window

  command {
    python /opt/call_variants.py ${salmon_all_sorted} \
        --smoothing_strategy ${smoothing_strategy} \
        --smoothing_window ${smoothing_window}
  }

  output {
    File salmon_all_events_ho = "salmon_all_events_ho.bed"
    File salmon_all_events_he = "salmon_all_events_he.bed"
  }

  runtime {
    docker: "commandlinegirl/sashimi:0.1"
  }
}

#########################################
# Integrate callers' outputs
#########################################

task integrate_outputs {
  File salmon_all_sorted
  Array[File] call_outputs_ho
  Array[File] call_outputs_he
  Array[Int] smoothing_windows
  Int? merge_distance
  Int? min_variant_len
  Float? min_score
  Array[File]? blacklist = []

  command {
    python /opt/integrate_outputs.py \
        --salmon_all_sorted ${salmon_all_sorted} \
        --call_outputs_ho $(echo ${sep=' ' call_outputs_ho}) \
        --call_outputs_he $(echo ${sep=' ' call_outputs_he}) \
        --smoothing_windows $(echo ${sep=' ' smoothing_windows}) \
        --merge_distance ${merge_distance} \
        --min_variant_len ${min_variant_len} \
        --min_score ${min_score} \
        --blacklist $(echo ${sep=' ' blacklist})

  }

  output {
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
    docker: "commandlinegirl/sashimi:0.1"
  }
}

#########################################
# Assess the accuracy of the caller
#########################################

task evaluate_output {
  File sample_events
  File? truth_events
  Float overlap = 0.5

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
    docker: "commandlinegirl/sashimi:0.1"
  }
}

#########################################
# Sashimi workflow
#########################################

workflow sashimi {

  # general inputs
  Boolean analyse_hom = true
  Boolean analyse_het = false
  # inputs to read_salmon_output
  File quantsf
  Array[String]? chromosomes
  # inputs to call_variants scatter & gather
  Array[String] smoothing_strategies = ["medianfilter"]
  Array[Int] smoothing_windows = [3]
  # inputs to integrate
  Int? merge_distance = 500
  Int? min_variant_len = 500
  Float? min_score = 1.0
  Array[File]? blacklist
  # inputs to evaluate
  Boolean evaluate = false
  File? truth_events_ho
  File? truth_events_he
  Float? overlap

  Array[Pair[String, Int]] smoothing_params = zip(smoothing_strategies, smoothing_windows)

  # read in data, preprocess, and normalize
  call read_salmon_output {
    input: 
      quantsf = quantsf,
      chromosomes = chromosomes
  }

  # scatter (run the call_variants with different parameters to call dels)
  scatter(pair in smoothing_params) {
    call call_variants {
      input: 
        salmon_all_sorted = read_salmon_output.salmon_all_sorted,
        smoothing_strategy = pair.left,
        smoothing_window = pair.right
    }
  }

  # gather and integrate outputs
  call integrate_outputs {
    input:
      salmon_all_sorted = read_salmon_output.salmon_all_sorted,
      call_outputs_ho = call_variants.salmon_all_events_ho,
      call_outputs_he = call_variants.salmon_all_events_he,
      smoothing_windows = smoothing_windows,
      merge_distance = merge_distance,
      min_variant_len = min_variant_len,
      min_score = min_score
  }

  if (evaluate && analyse_hom) {
    call evaluate_output as evaluate_hom {
      input:
        sample_events = integrate_outputs.result_dels_ho,
        truth_events = truth_events_ho,
        overlap = overlap
    }
  }

  if (evaluate && analyse_het) {
    call evaluate_output as evaluate_het {
      input:
        sample_events = integrate_outputs.result_dels_he,
        truth_events = truth_events_he,
        overlap = overlap
    }
  }

}

