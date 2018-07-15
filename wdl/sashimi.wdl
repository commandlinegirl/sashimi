
#########################################
# Preprocess salmon output
#########################################

task read_salmon_output {
  File quantsf
  File? blacklist
  Array[String]? chromosomes

  command {
    python /opt/read_salmon_output.py ${quantsf} \
        --chromosomes $(echo ${sep=' ' chromosomes})
  }

  output {
    File salmon_all_sorted = "salmon_all_sorted.bed" 
  }

  runtime {
    docker: "commandlinegirl/sashimi"
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
    docker: "commandlinegirl/sashimi"
  }
}

#########################################
# Integrate callers' outputs
#########################################

task integrate_outputs {
  File salmon_all_sorted
  Array[File] event_outputs
  Array[Int] smoothing_windows
  Int? merge_distance

  command {
    python /opt/integrate_outputs.py \
        --salmon_all_sorted ${salmon_all_sorted} \
        --event_outputs $(echo ${sep=' ' event_outputs}) \
        --smoothing_windows $(echo ${sep=' ' smoothing_windows}) \
        --merge_distance ${merge_distance}
  }

  output {
    File all_classifier_outputs = "all_classifier_outputs.bed"
    File integrated_output = "integrated_output.bed"
    File merged_calls_ho = "merged_calls_ho.bed"
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
  Float overlap

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

workflow sashimi {

  # inputs to read_salmon_output
  File quantsf
  Array[String]? chromosomes
  File? blacklist
  # inputs to call_variants scatter & gather
  Array[String] smoothing_strategies = ["medianfilter"]
  Array[Int] smoothing_windows = [3]
  # inputs to integrate
  Int? merge_distance = 500
  # inputs to evaluate
  Boolean evaluate = false
  File? truth_events
  Float overlap = 0.5

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
      event_outputs = call_variants.salmon_all_events_ho,
      smoothing_windows = smoothing_windows,
      merge_distance = merge_distance
  }

  if (evaluate) {
    call evaluate_output {
      input:
        sample_events = integrate_outputs.merged_calls_ho,
        truth_events = truth_events,
        overlap = overlap
    }
  }

}

