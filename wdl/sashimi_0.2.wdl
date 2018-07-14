
#########################################
# Preprocess salmon output
#########################################

task read_salmon_output {
  File quantsf
  File? blacklist
  Array[String]? chromosomes

  command {
    python /opt/read_salmon_output.py ${quantsf} \
        --blacklist ${blacklist} \
        --chromosomes ${chromosomes}
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
  Array[String]? smoothing_strategies
  Array[Int]? smoothing_windows
  Int? merge_distance

  command {
    python /opt/classifier.py ${salmon_all_sorted} \
        --smoothing_strategies $(echo ${sep=' ' smoothing_strategies}) \
        --smoothing_windows $(echo ${sep=' ' smoothing_windows}) \
        --merge_distance ${merge_distance}
  }

  output {
    File salmon_dels_only_ho = "salmon_dels_only_ho.bed"
    File salmon_dels_only_he = "salmon_dels_only_he.bed"
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
  Float overlap = 0.75

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

  # inputs to the classifier
  Array[String]? smoothing_strategies
  Array[Int]? smoothing_windows
  Int? merge_distance

  # inputs to evaluate
  File? truth_events
  Boolean evaluate = false

  # read in data and preprocess
  call read_salmon_output {
    input: 
      quantsf = quantsf,
      chromosomes = chromosomes
  }

  # run the classifier
  call classifier { 
    input: 
      salmon_all_sorted = read_salmon_output.salmon_all_sorted
      smoothing_strategies = smoothing_strategies,
      smoothing_windows = smoothing_windows,
      merge_distance = merge_distance
  }

 # if evaluate is true and HOM deletions were generated, run evaluations
 if (evaluate) {
   call evaluate_output {
     input:
       sample_events = integrate_outputs.salmon_dels_only_ho,
       truth_events = truth_events
   }
 }

}

