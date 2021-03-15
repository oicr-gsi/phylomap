version 1.0

workflow sarsCoV2heatmap {

  input {
  	String file_list
    Array[File] fasta_files = read_lines( file_list )
    File covariate_file
    File reference_fasta
    Float threshold
    File rename
  }

  parameter_meta {
    file_list: ""
    fasta_files: "A list of fasta files with consensus sequences of the reference, generated from different libraries"
    covariate_file: "A file with covariate information, in tabular format"
    reference_fasta: "Fasta file with the reference sequence"
    threshold: "Threshold for minimum coverage of reference genome"
    rename: ""
  }

  meta {
    author: "Lawrence Heisler"
    email: "lawrence.heisler@oicr.on.ca"
    description: "Generate variant calls from consensus sequence, a tree of distances, and a plot showing relationships with covariate data"
    dependencies: [
      {
        name: "ncov-tools/1.4",
        url: "https://github.com/jts/ncov-tools/archive/v1.4.tar.gz"
      },
    ]
  }

  call preprocess_consensus {
    input: threshold = threshold, fasta_files = fasta_files, rename = rename
  }
  
  call augur {
    input: fasta = preprocess_consensus.consensus_fa, ref = reference_fasta
  }

  call call_variants {
    input: ref = reference_fasta, aligned_fa = augur.aligned_fa
  }
  
  call plot {
    input: tree = augur.tree, variants=call_variants.alleles
  }

}

task preprocess_consensus {
  input {
    Float threshold
    Array[File] fasta_files
	  File rename
    String modules = "phylomap-tools/0"
  }

  parameter_meta {
      threshold: ""
      fasta_files: ""
      rename: ""
  }

  command <<<
    set -euo pipefail
    preprocess_fasta  \
    --threshold ~{threshold} \
    --fasta ~{sep=' --fasta ' fasta_files} > consensus.fasta \
    --stats "preprocess.stats.txt" \
    --rename_tbl ~{rename}
  >>>
 
   output {
     File consensus_fa = "consensus.fasta"
   }
}

task augur {
  input {
    File fasta
    File ref
    String modules = "ncov-tools/1 phylomap-tools/0"
  }

  command <<<
    set -euo pipefail
    augur align --sequences ~{fasta} --reference-sequence ~{ref} --output aligned.fasta --fill-gaps
    augur tree --alignment aligned.fasta --output tree_raw.nwk
    refid=$( head -n 1 ~{ref} | sed 's/>//' | sed 's/ .*//' )
    nw_reroot tree_raw.nwk $refid > tree.nwk
  >>>
  
  output {
    File aligned_fa = "aligned.fasta"
    File tree_raw = "tree_raw.nwk"
    File tree = "tree.nwk"
  }
}

task call_variants {
  input {
    File ref
    File aligned_fa
    String modules = "ncov-tools/1"
  }
  
  command <<<
    set -euo pipefail
    refid=$( head -n 1 ~{ref} | sed 's/>//' | sed 's/ .*//' )
    python ${NCOV_TOOLS_ROOT}/tree/align2alleles.py --reference-name $refid ~{aligned_fa} > alleles.tsv
  >>>
  
  output {
    File alleles = "alleles.tsv"
  }
  
}

task plot {
  input {
    File tree
    File variants
    String modules = "bis-rlibs ncov-tools/1 phylomap-tools/0"
  }
  
  command <<<
    set -euo pipefail
    script_R="/.mounts/labs/gsiprojects/gsi/lheisler/WDL/dev_cov2/scripts/plot_variant_tree.R"
    Rscript --vanilla $script_R --tree ~{tree}  --variants ~{variants} --out variant_tree.pdf
  >>>
  
  
  output {
    File plot1 = "variant_tree.pdf"
  }
}











