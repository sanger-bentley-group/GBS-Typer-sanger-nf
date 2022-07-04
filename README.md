# GBS-Typer-sanger-nf
The GBS Typer is for characterising Group B Strep by serotyping, resistance typing, MLST, surface protein typing and penicillin-binding protein typing. It has been adapted from [Ben Metcalf's GBS Typer pipeline](https://github.com/BenJamesMetcalf/GBS_Scripts_Reference) in Nextflow for portability and reproducibility.

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/sanger-pathogens/GBS-Typer-sanger-nf/blob/main/LICENSE)   
![build](https://github.com/sanger-pathogens/GBS-Typer-sanger-nf/workflows/build/badge.svg)
[![codecov](https://codecov.io/gh/sanger-pathogens/GBS-Typer-sanger-nf/branch/main/graph/badge.svg)](https://codecov.io/gh/sanger-pathogens/GBS-Typer-sanger-nf)

## Contents
- [ Running the serotyping and resistance typing pipeline on a local machine ](#local)
    - [ Installation ](#localinstall)
    - [ Usage ](#localusage)
- [ Running the serotyping and resistance typing pipeline on the farm at Sanger](#farm)
- [ Outputs ](#outputs)
    - [ Main Report ](#main)
    - [ Serotyping and Resistance Typing Output ](#serores)
    - [ MLST Output ](#mlst)
    - [ Surface Protein Typing Output ](#surfacetyper)
    - [ Pencillin-binding protein Typing Output ](#pbp)
- [ Other examples of running pipelines ](#examples)
- [ Additional options ](#additional)
- [ Troubleshooting for errors](#errors)
- [ Clean Up ](#cleanup)
- [ Other information ](#info)
    - [ Software dependencies ](#dependencies)
    - [ For developers ](#developers)
- [ Acknowledge us ](#acknowledge)
- [ Reporting issues ](#issues)

<a name="local"></a>
## Running the pipeline on a local machine

<a name="localinstall"></a>
### Installation
GBS Typer relies on Nextflow and Docker.
Download:
1. [Docker](https://www.docker.com/).
2. [Nextflow](https://www.nextflow.io/).
3. Clone repository:
```
git clone https://github.com/sanger-pathogens/GBS-Typer-sanger-nf.git
cd GBS-Typer-sanger-nf
```

Please note that if you need reproducibility for a research paper, you should use a specific version
(usually the latest). This can be done by adding a **-b <git tag version>** to the git clone command above.
Alternatively, you can use:
```
nextflow clone -r <git tag version> sanger-pathogens/GBS-Typer-sanger-nf
```
If running on an LSF head node you will need to bsub the **nextflow clone** command.

<a name="localusage"></a>
### Usage
Note: Running the pipeline requires an internet connection to allow the pipeline to automatically download its [dependencies image](https://hub.docker.com/repository/docker/sangerpathogens/gbs-typer-sanger-nf).

- To run on one sample (replacing `sampleID`)
```
nextflow run main.nf --run_sero_res --run_surfacetyper --run_mlst --reads 'data/sampleID_{1,2}.fastq.gz' --output 'sampleID'
```

- To run on multiple samples in a directory (replacing `output_file_prefix`)
```
nextflow run main.nf --run_sero_res --run_surfacetyper --run_mlst --reads 'data/*_{1,2}.fastq.gz' --output 'output_file_prefix'
```

<a name="farm"></a>
## Running the pipeline on the farm at Sanger:
Please refer to the internal Sanger wiki 'Pathogen Informatics Nextflow Pipelines' and then follow the steps below.

Note: Running the pipeline requires an internet connection and should be done in lustre storage.

**1. Clone repository**
```
git clone https://github.com/sanger-pathogens/GBS-Typer-sanger-nf.git
cd GBS-Typer-sanger-nf
```

**2. Load modules**
```
module load ISG/singularity
module load nextflow
```

**3. If running on farm5, you will need to set the http/https proxy**
```
export http_proxy=http://wwwcache.sanger.ac.uk:3128
export https_proxy=http://wwwcache.sanger.ac.uk:3128
```

**4. Run using bsub**

-  For a single sample (replacing `sampleID`)
```
bsub -G <your_team> -J <job_name> -o %J.out -e %J.err -R "select[mem>1000] rusage[mem=1000]" -M1000 "nextflow run main.nf --reads 'data/sampleID_{1,2}.fastq.gz' --run_sero_res --run_surfacetyper --run_mlst --output 'sampleID' -profile sanger,lsf"
```

- For multiple samples (replacing `output_file_prefix`)
```
bsub -G <your_team> -J <job_name> -o %J.out -e %J.err -R "select[mem>1000] rusage[mem=1000]" -M1000 "nextflow run main.nf --reads 'data/*_{1,2}.fastq.gz' --run_sero_res --run_surfacetyper --run_mlst --output 'output_file_prefix' -profile sanger,lsf"
```
Specifying `-profile sanger,lsf` will instruct Nextflow to run tasks as separate LSF jobs in parallel and will instruct the pipeline to build a local Singularity image from the [quay.io Docker image](https://quay.io/repository/sangerpathogens/gbs-typer-sanger-nf). The default is to run up to 20 jobs at a time. The default settings can be tuned to your requirements by editing the **lsf** profile within the nextflow.config file.

Add a **-N 'my-email-address'** to the end of the command line if you wish to be sent a report by email upon completion of the pipeline.

If you are processing many samples and would like to speed up the pipeline, you can increase the `queue_size` (default 100 i.e. a maximum 100 jobs running concurrently). For example, to increase the number of jobs running concurrently to 500, set `--queue_size 500`. (Note, the number of jobs may already be limited by your HPC's available resources.)

<a name="outputs"></a>
## Outputs

<a name="main"></a>
### Main Report
Running the command with `--run_sero_res`, `--run_surfacetyper` and `--run_mlst` will generate the main report `<sampleID>_gbs_typer_report.txt` for one sample or `<output_file_prefix>_gbs_typer_report.txt` for multiple samples. This will include the serotype, MLST type, allelic frequencies from MLST, resistance gene incidence, surface protein types and GBS-specific resistance variants.

<a name="serores"></a>
### Serotyping and Resistance Typing Output
**When specifying `--run_sero_res`**

This will create three tab-delimited files in a 'results' directory within the current directory. Specifying multiple samples will produce a row per sample:
1. **\<output_file_prefix\>_serotype_res_incidence.txt** - Gives the serotype and presence/absence (i.e. +/-) of antibiotic resistance genes (GBS-specific alleles and ResFinder/ARG-ANNOT genes)
e.g. Isolate Strep B sample 25292_2#105 has serotype II and have genes: 23S1, 23S3, GYRA, lsaC and tetM

lane_id | Serotype | 23S1 | 23S3 | CAT | ermB | ermT | FOSA | GYRA | lnuB | lsaC | mefA | MPHC | MSRA | msrD | PARC | RPOBGBS-1 | RPOBGBS-2 | RPOBGBS-3 | RPOBGBS-4 | SUL2 | tetB | tetL | tetM | tetO | tetS
:---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---:
25292_2#105 | II | pos | pos | neg | neg | neg | neg | pos | neg | pos | neg | neg | neg | neg | neg | neg | neg | neg | neg | neg | neg | neg | pos | neg | neg

2. **\<output_file_prefix\>_gbs_res_variants.txt** - Gives the SNP variants for GBS-specific resistance genes
e.g. Isolate Strep B sample 25292_2#105 have common variants 23S1, 23S3 and GYRA, but replacement of amino acid S by Q in position 17 of the PARC protein sequence

lane_id | 23S1 | 23S3 | GYRA | PARC | RPOBGBS-1 | RPOBGBS-2 | RPOBGBS-3 | RPOBGBS-4
:---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---:
25292_2#105 | * | * | * | Q17S | | | | |

3. **\<output_file_prefix\>_drug_cat_alleles_variants.txt** - Gives the GBS-specific variants and other resistance genes and alleles for drug categories: EC (macrolides, lincosamides, streptogramins or oxazolidinones), FQ (fluoroquinolones), OTHER (other antibiotics) and TET (tetracyclines)
e.g. Isolate Step B sample 25292_2#105 have GBS-specific variants: erythromycin-resistant 23S1 and 23S3, fluoroquinolone-resistant PARC and GYRA, and other resistance allele tetracycline-resistant tet(M)_1 of gene tet(M) (as specified by gene[allele])

lane_id | EC | FQ | OTHER | TET
:---: | :---: | :---: | :---: | :---:
25292_2#105 | 23S1:23S3 | PARC-Q17S:GYRA | neg | tet(M)[tet(M)_1]

<a name="mlst"></a>
### MLST Output
**When specifying `--run_mlst`**

This will produce a new MLST log file `<output_file_prefix>_new_mlst_alleles.log` indicating whether new MLST alleles have been found for each sample (where there are mismatches with sufficient read depth at least the value specified --mlst_min_read_depth [Default: 30]). If it includes "<lane_id>: New MLST alleles found." then a FASTA file for the corresponding sample `<sampleID>_new_mlst_alleles.fasta` and a pileup file `<sampleID>_new_mlst_pileup.txt` are generated.

For other samples that have no new MLST alleles and only have existing sequence types, these existing types are generated in `<output_file_prefix>_existing_sequence_types.txt`. If "None found" for a sample then no sequence types were found with sufficient read depth at least the value specified by --mlst_min_read_depth [Default: 30]).

<a name="surfacetyper"></a>
### Surface Protein Typing Output
**When specifying `--run_surfacetyper`**

This will create two tab-delimited files in the 'results' directory
1. **\<output_file_prefix\>_surface_protein_incidence.txt**
This shows the incidence of different surface protein alleles in the Strep B sample(s), e.g.

lane_id | alp1 | alp2/3 | alpha | hvgA | PI1 | PI2A1 | PI2A2 | PI2B | rib | srr1 | srr2
:---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---:
26189_8#338 | neg | pos | neg | neg | pos | pos | neg | neg | neg | pos | neg

2. **\<output_file_prefix\>_surface_protein_variants.txt**
This shows all the surface proteins in the Strep B sample(s), e.g.

lane_id | ALPH | hvgA | PILI | SRR
:---: | :---: | :---: | :---: | :---:
26189_8#338 | alp2/3 | neg | PI1:PI2A1 | srr1

<a name="pbp"></a>
### PBP (Penicillin-binding protein) Typing Output
To enable the PBP typing pipeline provide the **--run_pbptyper** command line argument and specify the contig FASTA files using and **--contigs**:
```
nextflow run main.nf --reads 'data/*_{1,2}.fastq.gz' --output 'output_file_prefix' --run_pbptyper --contigs 'data/*.fa'
```

If existing PBP alleles are found, a tab-delimited file is created in the 'results' directory. The file contains the sample IDs (that are determined from the contig FASTA file names e.g. 25292_2#85 from data/25292_2#85.fa), the contig identifiers with start, end and forward(+)/reverse(-) positions, and the PBP allele identifier.

lane_id | Contig | PBP_allele
:---: | :---: | :---:
25292_2#85 | .25292_2_85.9:39495-40455(+) | 2\|\|GBS_1A
26077_6#118 | .26077_6_118.11:39458-40418(+) | 1\|\|GBS_1A

If a new PBP allele is found in a sample, a FASTA file of amino acids is created in the 'results' directory. For example, if contig .25292_2_85.9 from sample 25292_2#85 contained a new PBP-1A allele, then .25292_2_85.9_GBS1A-1_PBP_new_allele.faa is generated with contents:

    >.26077_6_118.11:39458-40418(+)
    DIYNSDTYIAYPNNELQIASTIMDATNGKVIAQLGGRHQNENISFGTNQSVLTDRDWGST
    MKPISAYAPAIDSGVYNSTGQSLNDSVYYWPGTSTQLYDWDRQYMGWMSMQTAIQQSRNV
    PAVRALEAAGLDEAKSFLEKLGIYYPEMNYSNAISSNNSSSDAKYGASSEKMAAAYSAFA
    NGGTYYKPQYVNKIEFSDGTNDTYAASGSRAMKETTAYMMTDMLKTVLTFGTGTKAAIPG
    VAQAGKTGTSNYTEDELAKIEATTGIYNSAVGTMAPDENFVGYTSKYTMAIWTGYKNRLT
    PLYGSQLDIATEVYRAMMSY

<a name="examples"></a>
## Other examples of running pipelines
It is recommended you use the default parameters for specifying other resistance databases. However, to use different or multiple resistance databases with the GBS-specific resistance database, e.g. ARG-ANNOT and ResFinder in the db/0.0.2 directory, both with a minimum coverage of 70 and maximum divergence of 30:
```
nextflow run main.nf --reads 'data/*_{1,2}.fastq.gz' --output 'output_file_prefix' --other_res_dbs 'db/0.0.2/ARGannot-DB/ARG-ANNOT.fasta db/0.0.2/ResFinder-DB/ResFinder.fasta' --other_res_min_coverage '70 70' --other_res_max_divergence '30 30'
```
To run **only** the surface protein typing pipeline, use:
```
nextflow run main.nf --reads 'data/*_{1,2}.fastq.gz' --output 'output_file_prefix' --run_sero_res false --run_surfacetyper
```
To run **only** the MLST pipeline, use:
```
nextflow run main.nf --reads 'data/*_{1,2}.fastq.gz' --output 'output_file_prefix' --run_sero_res false --run_mlst
```
The default configuration will run serotyping and resistance typing only.

The default configuration will run serotyping and resistance typing only.
To run **only** the surface protein typing pipeline, use:
```
nextflow run main.nf --reads 'data/*_{1,2}.fastq.gz' --output 'output_file_prefix' --run_sero_res false --run_surfacetyper
```
To run **only** the MLST pipeline, use:
```
nextflow run main.nf --reads 'data/*_{1,2}.fastq.gz' --output 'output_file_prefix' --run_sero_res false --run_mlst
```
To run **only** the PBP typing pipeline, use:
```
nextflow run main.nf --output 'output_file_prefix' --run_sero_res false --run_pbptyper --contigs 'data/*.fa'
```
The **--reads** parameter is not needed for the PBP typing pipeline.


<a name="additional"></a>
## Additional options
### Inputs
    --contigs                       Path of file containing FASTA contigs. Only use when --run_pbptyper is specified. (Use wildcard '*' to specify multiple files, e.g. 'data/*.fa')
    --db_version                    Database version. (Default: 0.1.2)
    --other_res_dbs                 Path to other resistance reference database. Must be FASTA format. (Default: 'db/0.1.2/ARGannot-DB/ARG-ANNOT.fasta' from ARGannot_r3 of the SRST2 software, which includes non-redundant ResFinder and CARD database genes).

### Outputs
    --results_dir                   Results directory for output files. (Default: './results')

### Parameters
    --gbs_res_min_coverage          Minimum coverage for mapping to the GBS resistance database. (Default: 99.9)
    --gbs_res_max_divergence        Maximum divergence for mapping to the GBS resistance database. (Default: 5, i.e. report only hits with <5% divergence)
    --other_res_min_coverage        Minimum coverage for mapping to other resistance reference database(s). (Default: 70)
    --other_res_max_divergence      Maximum divergence for mapping to other resistance reference database. (Default: 30, i.e. report only hits with <30% divergence)
    --restyper_min_read_depth       Minimum read depth where mappings to antibiotic resistance genes with fewer reads are excluded. (Default: 30)
    --serotyper_min_read_depth      Minimum read depth where mappings to serotyping genes with fewer reads are excluded. (Default: 30)

### Other Pipeline Options
    --run_sero_res                  Run the main serotyping and resistance pipelines. (Default: true)
                                    Use '--run_sero_res false' to override the default.
    --run_mlst                      Run the MLST pipeline to query existing sequence types and new MLST alleles. (Default: false)
    --run_pbptyper                  Run the PBP (Penicillin binding protein) allele typer pipeline. Must also specify --contigs input. (Default: false)
    --run_surfacetyper              Run the surface protein typing pipeline. (Default: false)

### Other Pipeline Parameters
    --mlst_min_read_depth           Minimum read depth where mappings to alleles in MLST with fewer reads are excluded. Only operational with --run_mlst. (Default: 30)
    --pbp_frac_align_threshold      Minimum fraction of sequence alignment length of PBP gene. Only operational with --run_pbptyper. (Default: 0.5)
    --pbp_frac_identity_threshold   Minimum fraction of alignment identity between PBP genes and assemblies. Only operational with --run_pbptyper. (Default: 0.5)
    --surfacetyper_min_coverage     Minimum coverage for mapping to the GBS surface protein database. Only operational with --run_surfacetyper. (Default: 70)
    --surfacetyper_max_divergence   Maximum divergence for mapping to the GBS surface protein database. Only operational with --run_surfacetyper. (Default: 8, i.e. report only hits with <8% divergence)
    --surfacetyper_min_read_depth   Minimum read depth for surface protein typing pipeline. Only operational with --run_surfacetyper. (Default: 30)

Options can also be changed by editing the ```nextflow.config``` file.


<a name="errors"></a>
## Troubleshooting for errors
It is possible that the pipeline may not complete successfully due to issues with input files and/or individual steps of the pipeline. To troubleshoot potential issues, you can use the `-with-trace` parameter e.g. `nextflow run main.nf --reads 'data/*_{1,2}.fastq.gz' --output 'output_file_prefix' -with-trace` to create a trace file in the current directory while the pipeline is running. The trace file will provide the status of the current directory that provides the status and the location of the log files for each sample and step. For example, if the step failed in `[00/8803ea]`, the command and the error from this step can be viewed by `cat work/00/8803ea01f8bc0fb43f68335d82831f/.command.log` (Hint: while typing `work/00/8803ea`, the complete file path can be completed with the TAB button.)

Alternatively, if `-with-trace` was not provided in the original command, then you can also find the location of these directories in the bsub output e.g. `grep 'Error' pipeline.o | sort | uniq`.


<a name="cleanup"></a>
### Clean up

The `work` directory keeps the intermediate files of the pipeline. You can use `bsub.py 8 nextflow clean [run_name|session_id]` to clean up or `rm -rf work` only when no other pipelines are running (more dangerous). The run names and session ids can be found by using `nextflow log -q`


<a name="info"></a>
## Other information
<a name="dependencies"></a>
## Software dependencies
All pipeline dependencies are built into the [quay.io dependencies image](https://quay.io/repository/sangerpathogens/gbs-typer-sanger-nf), used by the pipeline.
The current program versions in this image are as follows:

Program | Version
:---: | :---:
bedtools | 2.29.2
biopython | 1.78
bowtie | 2.2.9
freebayes | 1.3.3+
prodigal | 1:2.6.3
python 2 | 2.7
python 3 | 3.8
samtools | 0.1.18
srst2 | 0.2.0

<a name="developers"></a>
## For developers
### Run unit tests
```
pytest tests
```

<a name="acknowledge"></a>
## Acknowledge us
If you use this software, please acknowledge the contributors as appropriate.

<a name="issues"></a>
## Reporting issues
If any questions or problems, please post them under [Issues](https://github.com/sanger-pathogens/GBS-Typer-sanger-nf/issues).
