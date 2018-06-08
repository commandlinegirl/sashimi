
task postprocess {
  File quantsf

  command {
    python /opt/sashimi.py ${quantsf}
  }

  output {
    File salmon_all_sorted = "quant.bed" 
  }

  runtime {
    docker: "commandlinegirl/sashimi"
  }

}

workflow sashimi {
  File quantsf = "quant.sf"

  call postprocess {
    input :
        quantsf = quantsf
  }

  output {
    File salmon_all_sorted = postprocess.salmon_all_sorted
  }
}
