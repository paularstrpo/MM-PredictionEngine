# predictionEngine

Prediction engine for MM Precision Medicine platform

## Install Dependencies

```
libs <- c('optparse', 'jsonlite', 'dplyr', 'RJSONIO', 'igraph', 'BiocManager')
install.packages(libs)
BiocManager::install('edgeR')
```

## General Usage

```
* --- [DAPHNI v2.0] Precision Medicine Drug Prediction Engine --- *

Options:
        --outdir=OUTDIR
                output directory path. [REQUIRED]

        --sampleID=SAMPLEID
                sample ID name for patient to run prediction engine on. [REQUIRED]

        --daphniDrugTable=DAPHNIDRUGTABLE
                Table of drugs and their related annotations to consider for predictions. [REQUIRED]

        --daphniMutationCensus=DAPHNIMUTATIONCENSUS
                Table of known activating / deactivating mutations for drugs considered in this tool. [REQUIRED]

        --daphniBiomarkerTable=DAPHNIBIOMARKERTABLE
                Table of manually curated drug-gene assocations for tier 1 drug buckets. [REQUIRED]

        --daphniPrognosticMarkers=DAPHNIPROGNOSTICMARKERS
                Table of manually curated prognostic biomarker information. [REQUIRED]

        --cytobandCoordinates=CYTOBANDCOORDINATES
                Table with per-arm cytoband coordinates [REQUIRED]

        --zScoreTable=ZSCORETABLE
                Table with expression zscores for all patients, including the patient of interest.

        --treeFile=TREEFILE
                json file with clonal tree structure.

        --somMutFileVICC=SOMMUTFILEVICC
                VICC annotations for SNV results, with clone annotations.

        --somMutFileCIVIC=SOMMUTFILECIVIC
                CIVIC annotations for SNV results, with clone annotations.

        --somMutFile=SOMMUTFILE
                annotated consensus vcf for somatic mutations

        --cnaFileVICC=CNAFILEVICC
                VICC annotations for CNV results, with clone annotations.

        --cnaFileCIVIC=CNAFILECIVIC
                CIVIC annotations for CNV results, with clone annotations.

        --cnaFile=CNAFILE
                gene-level CNA file from facets

        --geneExprFileVICC=GENEEXPRFILEVICC
                VICC annotations for expressed genes.

        --geneExprFileCIVIC=GENEEXPRFILECIVIC
                CIVIC annotations for expressed genes.

        --facetsFile=FACETSFILE
                facets cncf file

        --translocFile=TRANSLOCFILE
                predicted translocations

        --scarFile=SCARFILE
                Table containing scar scores at cohort level, including the patient of interest.

        --selineScoresFile=SELINESCORESFILE
                Table containing selinexor signature data at cohort level, including the patient of interest.

        -h, --help
                Show this help message and exit
```

## Example Command

``` 
docker build -t predengine .
docker run -v ${PWD}/example_outputs/:/data/ -v ${PWD}/example_inputs/:/example_inputs/ -v ${PWD}/refData/:/refData/ predengine \
Rscript predEngine.R \
--outdir '/data' \
--sampleID 'ISMMS01' \
--daphniDrugTable '/refData/daphniDrugTableSourcesAnnotated.tsv' \
--daphniPrognosticMarkers '/refData/prognosticMarkersTable.tsv' \
--daphniMutationCensus '/refData/daphni_mutation_census.tsv' \
--daphniBiomarkerTable '/refData/daphni_tier1_biomarkers.tsv' \
--cytobandCoordinates '/refData/hg38_cytoband_coordinates_perarm.tsv' \
--selineScoresFile "/example_inputs/selinescores_latest.tsv" \
--zScoreTable '/example_inputs/zscores_latest.csv' \
--scarFile "/example_inputs/scarScores_latest.tsv" \
--geneExprFileVICC "/example_inputs/vicc_expression_data.csv" \
--geneExprFileCIVIC "/example_inputs/civic_expression_data.csv" \
--treeFile "/example_inputs/sample.summ.json.gz" \
--somMutFile "/example_inputs/annotated_variants.consensus.vcf" \
--somMutFileCIVIC "/example_inputs/civic_mutation_clone_data.txt" \
--somMutFileVICC "/example_inputs/vicc_mutation_clone_data.txt" \
--cnaFile "/example_inputs/output_CNV_file.txt" \
--cnaFileCIVIC "/example_inputs/civic_cnv_clone_data.txt" \
--cnaFileVICC "/example_inputs/vicc_cnv_clone_data.txt" \
--facetsFile "/example_inputs/tumor.facets_cncf.txt" \
--translocFile "/example_inputs/predicted_translocations.csv"
```

## Input Files

Reference files for the `[daphniDrugTable]`, `[daphniMutationCensus]`, and `[daphniBiomarkerTable]`, arguments are provided under the `refData/` folder in this repository.

## Output Files

Example output files are shown in the `example_outputs/` folder in this directory