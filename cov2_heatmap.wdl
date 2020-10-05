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
    consensus_fasta: "a list of fasta files with consensus sequences of the reference,generated from different libraries"
    reference_fasta: "fasta file with the reference sequence"
    covariates: "a file with covariate information, in tabular format"
  }

  meta {
    author: "Lawrence Heisler"
    email: "lawrence.heisler@oicr.on.ca"
    description: "Generate variant calls from consensus sequence, a tree of distances, and a plot showing relationships with covariate data"
    dependencies: [
      {
        name: "ncov-tools/1",
        url: ""
      },
    ]
  }

  call preprocess_consensus {
    input: threshold = threshold,fasta_files = fasta_files,rename = rename
  }
  
  call augur {
    input: fasta = preprocess_consensus.consensus_fa, ref = reference_fasta
  }

  call call_variants {
    input: ref = reference_fasta, aligned_fa = augur.aligned_fa
  }
  
  call lineage {
    input: fasta = preprocess_consensus.consensus_fa
  }
  
  call plot {
    input: tree = augur.tree,variants=call_variants.alleles
  }

}


task plot{
  input{
    File tree
    File variants
    String modules = "bis-rlibs ncov-tools/1"
  }
  
  command <<<
    set -euo pipefail
    module purge
    module load ~{modules}
    script_R="/.mounts/labs/gsiprojects/gsi/lheisler/WDL/dev_cov2/scripts/plot_variant_tree.R"
    Rscript --vanilla $script_R --tree ~{tree}  --variants ~{variants} --out variant_tree.pdf
  >>>
  
  
  output{
    File plot1 = "variant_tree.pdf"
  }
}



task lineage{

  input{
    File fasta
    String modules = "ncov-tools/1"
  }

  command <<<
    set -euo pipefail
    module purge
    module load ~{modules}  
    pangolin -t 8 ~{fasta} -o .
    cat lineage_report.csv | sed 's/taxon,/name,/' | cut -d "," -f 1,2 | sed 's/,/\t/' > lineage_covariate.tsv
  >>>
  
  output{
    File lineage_report = "lineage_report.csv"
    File lineage_covariates = "lineage_covariate.tsv"

  }


}

task call_variants{
  input{
    File ref
    File aligned_fa
    String modules = "ncov-tools/1"
  }
  
  command <<<
    set -euo pipefail
    module purge
    module load ~{modules}  
    refid=$( head -n 1 ~{ref} | sed 's/>//' | sed 's/ .*//' )
    python ${NCOV_TOOLS_ROOT}/tree/align2alleles.py --reference-name $refid ~{aligned_fa} > alleles.tsv
  >>>
  
  output{
    File alleles = "alleles.tsv"
  
  }
  
}


task augur {
  input{
    File fasta
    File ref
    String modules = "ncov-tools/1"
  }

  command <<<
    set -euo pipefail
    module purge
    module load ~{modules}
    augur align --sequences ~{fasta} --reference-sequence ~{ref} --output aligned.fasta --fill-gaps
    augur tree --alignment aligned.fasta --output tree_raw.nwk
    refid=$( head -n 1 ~{ref} | sed 's/>//' | sed 's/ .*//' )
    nw_reroot tree_raw.nwk $refid > tree.nwk
  >>>
  
  output{
    File aligned_fa = "aligned.fasta"
    File tree_raw = "tree_raw.nwk"
    File tree = "tree.nwk"
  }
}




task preprocess_consensus {
  input{
    Float threshold
    Array[File] fasta_files
	File rename
  }

  command <<<
    set -euo pipefail
    /.mounts/labs/gsiprojects/gsi/lheisler/WDL/dev_cov2/scripts/preprocess_fasta.pl  \
      --threshold ~{threshold} \
      --fasta ~{sep=' --fasta ' fasta_files} > consensus.fasta \
      --stats "preprocess.stats.txt" \
      --rename_tbl ~{rename}
    
  >>>
 
   output{
     File consensus_fa = "consensus.fasta"
   }
}


