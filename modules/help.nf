
def printHelp() {
  log.info"""
  Usage:
    nextflow run main.nf --reads [paired reads] --output [output prefix]

  Description:
    Serotyping and resistance typing of Group B Strep.

  Additional options:

    Inputs:
      --contigs                     Path of file containing FASTA contigs. Only use when
                                    --run_pbptyper is specified. (Use wildcard '*' to specify multiple files,
                                    e.g. 'data/*.fa')
      --db_version                  Database version. (Default: 0.1.2)
      --other_res_dbs               Paths to other resistance reference database(s). Must be FASTA format.
                                    Specify 'none' to omit using other resistance databases.
                                    (Default: 'db/0.1.0/ARGannot-DB/ARG-ANNOT.fasta'
                                    from ARGannot_r3 of the SRST2 software, which includes non-redundant
                                    ResFinder and CARD database genes).

    Outputs:
      --results_dir                  Results directory for output files. (Default: './results')
      --tmp_dir                      Temporary directory for intermediate files. (Default: './tmp')

    Parameters:
      --gbs_res_min_coverage         Minimum coverage for mapping to the GBS resistance database. (Default: 99.9)
      --gbs_res_max_divergence       Maximum divergence for mapping to the GBS resistance database. (Default: 5,
                                     i.e. report only hits with <5% divergence)
      --other_res_min_coverage       Minimum coverage for mapping to other resistance reference database(s).
                                     Number of values must equal the number of resistance reference database files
                                     and must correspond to the order specified in --other_res_dbs. (Default: 70)
      --other_res_max_divergence     Maximum divergence for mapping to other resistance reference database(s).
                                     Number of values must equal the number of resistance reference database files
                                     and must correspond to the order specified in --other_res_dbs. (Default: 30,
                                     i.e. report only hits with <30% divergence)
      --restyper_min_read_depth      Minimum read depth where mappings to antibiotic resistance genes with fewer
                                     reads are excluded. (Default: 30)
      --serotyper_min_read_depth     Minimum read depth where mappings to serotyping genes with fewer reads are excluded.
                                     (Default: 30)

    Other Pipeline Options:
      --run_sero_res                Run the main serotyping and resistance pipelines. (Default: true)
                                    Use '--run_sero_res false' to override the default.
      --run_mlst                    Run the MLST pipeline to query new MLST alleles. (Default: false)
      --run_pbptyper                Run the PBP (Penicillin binding protein) allele typer pipeline. Must also specify
                                    --contigs input. (Default: false)
      --run_surfacetyper            Run the surface protein typing pipeline. (Default: false)


    Other Pipeline Parameters:
      --mlst_min_read_depth         Minimum read depth where mappings to alleles in MLST with fewer reads are excluded.
                                    Only operational with --run_mlst. (Default: 30)
      --pbp_frac_align_threshold    Minimum fraction of sequence alignment length of PBP gene. Only operational with
                                    --run_pbptyper. (Default: 0.5)
      --pbp_frac_identity_threshold Minimum fraction of alignment identity between PBP genes and assemblies. Only
                                    operational with --run_pbptyper. (Default: 0.5)
      --surfacetyper_min_coverage   Minimum coverage for mapping to the GBS surface protein database. Only
                                    operational with --run_surfacetyper. (Default: 70)
      --surfacetyper_max_divergence Maximum divergence for mapping to the GBS surface protein database. Only
                                    operational with --run_surfacetyper. (Default: 8, i.e. report only hits with <8%
                                    divergence)
      --surfacetyper_min_read_depth Minimum read depth for surface protein typing pipeline. Only operational with
                                    --run_surfacetyper. (Default: 30)
  """.stripIndent()
}
