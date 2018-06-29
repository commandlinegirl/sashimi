
#########################################
# Preprocess and normalize salmon output
#########################################

task read_salmon_output {
  File quantsf
  File?  blacklist
  Array[String]? chromosomes
  Boolean smooth_raw_output = false
  Int? smoothing_window_len
  String? smoothing_strategy

  command {
    python /opt/read_salmon_output.py ${quantsf}
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

task classifier {
  File salmon_all_sorted
  Int neighbor_size = 4

  command {
    python /opt/classifier.py ${salmon_all_sorted} \
        --neighbor_size ${neighbor_size} 
  }

  output {
    File salmon_dels_only_ho = "salmon_dels_only_ho.bed"
    File salmon_dels_only_he = "salmon_dels_only_he.bed"
    File salmon_all_events = "salmon_all_events.bed"
  }

  runtime {
    docker: "commandlinegirl/sashimi"
  }
}

#########################################
# Merge call outputs from different calls
#########################################

task merge_classifier_outputs {
  File salmon_all_sorted
  Array[File] event_outputs

  command {
    python /opt/merge_classifier_outputs2.py \
        --salmon_all_sorted ${salmon_all_sorted} \
        --event_outputs $(echo ${sep=' ' event_outputs})
  }

  output {
    File merged_events = "merged_events.bed" 
  }

  runtime {
    docker: "commandlinegirl/sashimi"
  }
}

#########################################
# Integrate outputs
#########################################

task integrate_outputs {
  File merged_events 
  Array[Int] neighbor_sizes

  command {
    python /opt/integrate_outputs.py \
        --merged_events ${merged_events} \
        --neighbor_sizes $(echo ${sep=' ' neighbor_sizes})
  }

  output {
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
    docker: "commandlinegirl/sashimi"
  }
}

#########################################
# Sashimi workflow
#########################################

workflow sashimi {
  File quantsf
  File? truth_events
  Array[Int] neighbor_sizes = [6, 10]
  Boolean evaluate = false

  # read in data, preprocess, and normalize
  call read_salmon_output {
    input: 
      quantsf = quantsf
  }

  # scatter (run the classifier with different parameters to call dels)
  scatter(ns in neighbor_sizes) {
    call classifier { 
      input: 
        neighbor_size = ns,
        salmon_all_sorted = read_salmon_output.salmon_all_sorted
    }
  }

  # gather
  call merge_classifier_outputs as merge_outputs {
    input:
      salmon_all_sorted = read_salmon_output.salmon_all_sorted,
      event_outputs = classifier.salmon_all_events
  }

  # integrate outputs
  call integrate_outputs {
    input:
        merged_events = merge_outputs.merged_events,
        neighbor_sizes = neighbor_sizes
  }

 if (evaluate) {
   call evaluate_output {
     input:
       sample_events = integrate_outputs.merged_calls_ho,
       truth_events = truth_events
   }
 }

}

