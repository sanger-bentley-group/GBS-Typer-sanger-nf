
def printHelp() {
  log.info"""
  Usage:
    nextflow run mic.nf --output [output prefix]

  Description:
    MIC predictions for Group B Strep.
    Serotyping, resistance and PBP typing pipelines must have been successfully run prior to running this pipeline.

  Options:

    --results_dir                 Results directory. (Default: './results')
    --output                      Output prefix that will be used to identify input files in the results directory.

  """.stripIndent()
}
