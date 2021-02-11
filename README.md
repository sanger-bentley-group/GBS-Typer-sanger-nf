# GBS-Typer-sanger-nf
An updated NextFlow version of [Ben Metcalf's GBS Typer pipeline](https://github.com/BenJamesMetcalf/GBS_Scripts_Reference).

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/sanger-pathogens/GBS-Typer-sanger-nf/blob/main/LICENSE)   
![build](https://github.com/sanger-pathogens/GBS-Typer-sanger-nf/workflows/build/badge.svg)  
![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/sangerpathogens/gbs-typer-sanger-nf)   
[![codecov](https://codecov.io/gh/sanger-pathogens/GBS-Typer-sanger-nf/branch/main/graph/badge.svg)](https://codecov.io/gh/sanger-pathogens/GBS-Typer-sanger-nf)   


## Running the pipeline on a local machine

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

### Usage
Note: Running the pipeline requires an internet connection to allow the pipeline to automatically download its [dependencies image](https://hub.docker.com/repository/docker/sangerpathogens/gbs-typer-sanger-nf).

#### To run on one sample:
```
nextflow run main.nf --reads 'data/sampleID_{1,2}.fastq.gz' --output 'sampleID'
```

#### To run on multiple samples in a directory:
```
nextflow run main.nf --reads 'data/*_{1,2}.fastq.gz' --output 'output_file_prefix'
```


## Running the pipeline at Sanger:
Please refer to the internal Sanger wiki 'Pathogen Informatics Nextflow Pipelines' and then follow the steps below.

Note: Running the pipeline requires an internet connection and should be done in lustre storage.

1. Clone repository
```
git clone https://github.com/sanger-pathogens/GBS-Typer-sanger-nf.git
cd GBS-Typer-sanger-nf
```

2. Load modules
```
module load ISG/singularity
module load nextflow
```

3. For a single sample, run using bsub and add '-profile singularity' as an option, e.g.
```
bsub -G <your_team> -J <job_name> -o %J.out -e %J.err -R "select[mem>1000] rusage[mem=1000]" -M1000 "nextflow run main.nf --reads 'data/sampleID_{1,2}.fastq.gz' --output 'sampleID' -profile singularity"
```

4. For multiple samples, also run using bsub and add '-profile singularity,lsf', e.g.
```
bsub -G <your_team> -J <job_name> -o %J.out -e %J.err -R "select[mem>1000] rusage[mem=1000]" -M1000 "nextflow run main.nf --reads 'data/*_{1,2}.fastq.gz' --output 'output_file_prefix' -profile singularity,lsf"
```
This will instruct Nextflow to run tasks as separate LSF jobs in parallel and can be significantly faster. The default is to run up to 20 jobs at a time. The default settings can be tuned to your requirements by editing the **lsf** profile within the nextflow.config file.

Specifying the **singularity** profile will instruct the pipeline to build a local singularity image from the [docker hub dependencies image](https://hub.docker.com/repository/docker/sangerpathogens/gbs-typer-sanger-nf).

Add a **-N 'my-email-address'** to the end of the command line if you wish to be sent a report by email upon completion of the pipeline.


## Output

### One sample
Using command:
```
nextflow run main.nf --reads 'data/sampleID_{1,2}.fastq.gz' --output 'sampleID'
```

This will create three tab-delimited files in a 'results' directory within the current directory:
1. **sampleID_serotype_res_incidence.txt** - Gives the serotype and presence/absence (i.e. +/-) of antibiotic resistance genes (GBS-specific alleles and ResFinder/ARG-ANNOT genes)
e.g. Isolate Strep B sample 25292_2#105 has serotype II and have genes: 23S1, 23S3, GYRA, LSAC and TETM

ID | Serotype | 23S1 | 23S3 | CAT | ERMB | ERMT | FOSA | GYRA | LNUB | LSAC | MEFA | MPHC | MSRA | MSRD | PARC | RPOBGBS-1 | RPOBGBS-2 | RPOBGBS-3 | RPOBGBS-4 | SUL2 | TETB | TETL | TETM | TETO | TETS
:---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---:
25292_2#105 | II | + | + | - | - | - | - | + | - | + | - | - | - | - | - | - | - | - | - | - | - | - | + | - | -

2. **sampleID_gbs_res_variants.txt** - Gives the SNP variants for GBS-specific resistance genes
e.g. Isolate Strep B sample 25292_2#105 have common variants 23S1, 23S3 and GYRA, but replacement of amino acid S by Q in position 17 of the PARC protein sequence

ID | 23S1 | 23S3 | GYRA | PARC | RPOBGBS-1 | RPOBGBS-2 | RPOBGBS-3 | RPOBGBS-4
:---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---:
25292_2#105 | 23S1 | 23S3 | GYRA | PARC-Q17S | | | | |

3. **sampleID_drug_cat_alleles_variants.txt** - Gives the GBS-specific variants and other resistance genes and alleles for drug categories: EC (macrolides, lincosamides, streptogramins or oxazolidinones), FQ (fluoroquinolones), OTHER (other antibiotics) and TET (tetracyclines)
e.g. Isolate Step B sample 25292_2#105 have GBS-specific variants: erythromycin-resistant 23S1 and 23S3, fluoroquinolone-resistant PARC and GYRA, and other resistance allele tetracycline-resistant tet(M)_1 of gene tet(M) (as specified by gene[allele])

ID | EC | FQ | OTHER | TET
:---: | :---: | :---: | :---: | :---:
25292_2#105 | 23S1:23S3 | PARC-Q17S:GYRA | neg | tet(M)[tet(M)_1]

### Multiple samples
Using command:
```
nextflow run main.nf --reads 'data/*_{1,2}.fastq.gz' --output 'output_file_prefix'
```
This will produce combined tables of output_file_prefix_serotype_res_incidence.txt, output_file_prefix_gbs_res_variants.txt and output_file_prefix_drug_cat_alleles_variants.txt in the 'results' directory that can be identified by sample ID (i.e. the name of the file before _1.fastq.gz or _2.fastq.gz).

## Additional options
### Inputs
    --contigs                       Path of file containing FASTA contigs. Only use when --run_pbptyper is specified. (Use wildcard '*' to specify multiple files, e.g. 'data/*.fa')
    --db_version                    Database version. (Default: 0.1.2)
    --other_res_dbs                 Paths to other resistance reference database(s). Must be FASTA format. Specify 'none' to omit using other resistance databases. (Default: 'db/0.1.1/ARGannot-DB/ARG-ANNOT.fasta' from ARGannot_r3 of the SRST2 software, which includes non-redundant ResFinder and CARD database genes).

### Outputs
    --results_dir                   Results directory for output files. (Default: './results')
    --tmp_dir                       Temporary directory for intermediate files. (Default: './tmp')

### Parameters
    --gbs_res_min_coverage          Minimum coverage for mapping to the GBS resistance database. (Default: 99.9)
    --gbs_res_max_divergence        Maximum divergence for mapping to the GBS resistance database. (Default: 5, i.e. report only hits with <5% divergence)
    --other_res_min_coverage        Minimum coverage for mapping to other resistance reference database(s). Number of values must equal the number of resistance reference database files and must correspond to the order specified in --other_res_dbs. (Default: 70)
    --other_res_max_divergence      Maximum divergence for mapping to other resistance reference database(s). Number of values must equal the number of resistance reference database files and must correspond to the order specified in --other_res_dbs. (Default: 30, i.e. report only hits with <30% divergence)
    --restyper_min_read_depth       Minimum read depth where mappings to antibiotic resistance genes with fewer reads are excluded. (Default: 30)
    --serotyper_min_read_depth      Minimum read depth where mappings to serotyping genes with fewer reads are excluded. (Default: 30)

### Other Pipeline Options
    --run_sero_res                  Run the main serotyping and resistance pipelines. (Default: true)
                                    Use '--run_sero_res false' to override the default.
    --run_mlst                      Run the MLST pipeline to query new MLST alleles. (Default: false)
    --run_pbptyper                  Run the PBP (Penicillin binding protein) allele typer pipeline. Must also specify --contigs input. (Default: false)
    --run_surfacetyper              Run the surface protein typing pipeline. (Default: false)

### Other Pipeline Parameters
    --mlst_min_read_depth           Minimum read depth where mappings to alleles in MLST with fewer reads are excluded. Only operational with --run_mlst. (Default: 30)
    --pbp_frac_align_threshold      Minimum fraction of sequence alignment length of PBP gene. Only operational with --run_pbptyper. (Default: 0.5)
    --pbp_frac_identity_threshold   Minimum fraction of alignment identity between PBP genes and assemblies. Only operational with --run_pbptyper. (Default: 0.5)
    --surfacetyper_min_coverage     Minimum coverage for mapping to the GBS surface protein database. Only operational with --run_surfacetyper. (Default: 70)
    --surfacetyper_max_divergence   Maximum divergence for mapping to the GBS surface protein database. Only operational with --run_surfacetyper. (Default: 8, i.e. report only hits with <8% divergence)
    --surfacetyper_min_read_depth   Minimum read depth for surface protein typing pipeline. Only operational with --run_surfacetyper. (Default: 30)

All default options can be changed by editing the ```nextflow.config``` file.

## Other Pipelines
### MLST Pipeline Usage
To include MLST pipeline to query new MLST alleles
```
nextflow run main.nf --reads 'data/*_{1,2}.fastq.gz' --output 'output_file_prefix' --new_mlst
```

This will produce a text file including sequences and pileups for each allele with sufficient read depth greater than or equal to the --mlst_min_read_depth for each sample called <sampleID>_new_mlst_alleles.txt in the 'results' directory.

    New MLST Allele Consensus:
    >adhP_1
    CCAGGACGCATTTTGGGTCACGAAGGCATTGGTATAGTAGAAGAAATTGGAGAAGGCGTA
    ACGTCTTTGAGGGTTGGTGATCGTGTCTCTATTGCATGGTTCTTTGAAGGGTGCGGTCAT
    TGCGAATACTGTACTACAGGACGTGAGACACTTTGTCGTAGTGTTAAAAATGCTGGATAC
    AGTGTTGATGGTGGTATGAGTGAATACGCTATTGTTACCGCGGACTATGCGGTTAAGGTT
    CCTGAGGGATTAGACCCAGCTCAAGCATCATCAATCACTTGTGCTGGAGTAACAACATAC
    AAGGCTATCAAAGAAGCTGGAGCTGCTCCTGGTCAGTGGATTGCAGTGTATGGTGCAGGT
    GGTCTTGGAAACTTAGCAGTCCAATATGCAAAAAAAGTATTCAATGCTCATGTTGTAGCT
    GTTGATATTAACGCAGATAAACTTCAATTAGCTAAAGAGGTTGGAGCAGATTTGACAGTT
    AATGGCAAAGAAATAAAA

    New MLST Allele Pileup:
    adhP_1	1	C	84	^#,^#.^".^".^".^".^".^".^".^".^".^".^".^".^".^".^".^".^".^".^".^".^".^".^".^".^".^,.^".^".^".^".^".^".^".^".^".^".^".^".^".^".^".^,.^".^".^".^".^#.^".^".^".^".^".^".^".^".^".^,.^".^".^".^".^".^".^,.^,.^".^".^".^,.^".^".^".^".^".^".^".^,.^".^".^".^,.^".	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

### Surface Protein Typing Pipeline Usage
To enable the surface typing pipeline provide the **--run_surfacetyper** command line argument:
```
nextflow run main.nf --reads 'data/*_{1,2}.fastq.gz' --output 'output_file_prefix' --run_surfacetyper
```

This will create two tab-delimited files in the 'results' directory
1. **<output_file_prefix>_surface_protein_incidence.txt**
This shows the incidence of different surface protein alleles in the Strep B sample(s), e.g.

ID | ALP1 | ALP23 | ALPHA | HVGA | PI1 | PI2A1 | PI2A2 | PI2B | RIB | SRR1 | SRR2
:---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---:
26189_8#338 | - | + | - | - | + | + | - | - | - | + | -

2. **<output_file_prefix>_surface_protein_variants.txt**
This shows all the surface proteins in the Strep B sample(s), e.g.

ID | ALPH | HVGA | PILI | SRR
:---: | :---: | :---: | :---: | :---:
26189_8#338 | ALP23 | neg | PI1:PI2A1 | SRR1

### PBP (Penicillin-binding protein) Typing Pipeline Usage
To enable the PBP typing pipeline provide the **--run_pbptyper** command line argument and specify the contig FASTA files using and **--contigs**:
```
nextflow run main.nf --reads 'data/*_{1,2}.fastq.gz' --output 'output_file_prefix' --run_pbptyper --contigs 'data/*.fa'
```

The contig fasta files must be named <sample_id>.fa.
If existing PBP alleles are found, a tab-delimited file is created in the 'results' directory. The file contains the sample IDs (that are determined from the contig FASTA file names e.g. 25292_2#85 from data/25292_2#85.fa), the contig identifiers with start, end and forward(+)/reverse(-) positions, and the PBP allele identifier.

ID | Contig | PBP_allele
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


### MIC (Minimum Inhibitory Concentration) Prediction
The MIC prediction pipeline can be run standalone after serotyping, resistance and PBP typing have been run successfully.
To run MIC prediction, use the following command:
```
nextflow run mic.nf --output 'output_file_prefix'
```
The output prefix must match that used for the prior typing pipeline run.
This will create one tab-delimited **<output_file_prefix>_mic_predictions.txt** file in the 'results' directory. 
e.g.

ID | PBP | TET | EC | FQ | OTHER
:---: | :---: | :---: | :---: | :---: | :---:
26189_8#338 | Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,NA,NA,NA,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag,Flag | neg,<=,1,S,<=,1,S,<=,1,U,<=,4,S,<=,0.5,U | ErmB_MLS[ErmB_841]:23S1:23S3,>=,1,R,>=,1,R,<=,2.0,S,<=,1.0,S,pos | PARC:GYRA,NA,NA,NA,Flag,Flag,Flag | TetM_Tet[TetM_831],>=,8,R
        

### Examples
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


## Dependencies
All pipeline dependencies are built into the [docker hub dependencies image](https://hub.docker.com/repository/docker/sangerpathogens/gbs-typer-sanger-nf), used by the pipeline.
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

## For developers
### Run unit tests
```
python3 -m pytest
```

### Making a release
Once all changes have been pushed to the main branch, confirm that the CI has run successfully and then execute:
```
release.sh <version number (without v)>
```
This will tag the main branch with the specified version number.
