# # OPTION PARSER (INPUT FILE SPECS)
option_list <- list(
  optparse::make_option("--outdir", type='character', help='output directory path. [REQUIRED]', default='.'),
  optparse::make_option("--sampleID", type='character', help='sample ID name for patient to run prediction engine on. [REQUIRED]', default=NA),
  optparse::make_option("--daphniDrugTable", type='character', help='Table of drugs and their related annotations to consider for predictions. [REQUIRED]', default=NA),
  optparse::make_option("--daphniMutationCensus", type='character', help='Table of known activating / deactivating mutations for drugs considered in this tool. [REQUIRED]', default=NA),
  optparse::make_option("--daphniBiomarkerTable", type='character', help='Table of manually curated drug-gene assocations for tier 1 drug buckets. [REQUIRED]', default=NA),
  optparse::make_option("--zScoreTable", type='character', help='Table with expression zscores for all patients, including the patient of interest.', default=NA),
  optparse::make_option("--treeFile", type='character', help='json file with clonal tree structure.', default=NA),
  optparse::make_option("--somMutFileVICC", type='character', help='VICC annotations for SNV results, with clone annotations.', default=NA),
  optparse::make_option("--somMutFileCIVIC", type='character',help='CIVIC annotations for SNV results, with clone annotations.', default=NA),
  optparse::make_option("--somMutFile", type='character',help='annotated consensus vcf for somatic mutations', default=NA),
  optparse::make_option("--cnaFileVICC", type='character', help='VICC annotations for CNV results, with clone annotations.', default=NA),
  optparse::make_option("--cnaFileCIVIC",type='character', help='CIVIC annotations for CNV results, with clone annotations.', default=NA),
  optparse::make_option("--cnaFile", type='character', help='gene-level CNA file from facets', default=NA),
  optparse::make_option("--geneExprFileVICC", type='character', help='VICC annotations for expressed genes.', default=NA),
  optparse::make_option("--geneExprFileCIVIC", type='character', help='CIVIC annotations for expressed genes.',  default=NA),
  optparse::make_option("--mmPSNFile", type='character', help='MM-PSN Classifier results.',  default=NA),
  optparse::make_option("--mmHallMarks", type='character', help = 'mm-hallmarks results with cytogenetic markers', default=NA),
  optparse::make_option("--scarFile", type='character',help='Table containing scar scores at cohort level, including the patient of interest.', default=NA),
  optparse::make_option("--selineScoresFile", type='character',help='Table containing selinexor signature data at cohort level, including the patient of interest.', default=NA)
)

# # get command line options, if help option encountered print help and exit.
opt <- optparse::parse_args(
  optparse::OptionParser(option_list = option_list, description = "* --- [DAPHNI v2.0] Precision Medicine Drug Prediction Engine --- *")
  )

# load libraries
options(stringsAsFactors = FALSE, warn=-1)
libs <- c('jsonlite', 'tidyverse', 'RJSONIO', 'igraph', 'edgeR')
invisible(lapply(libs, library, character.only=TRUE))

# --- --- --- #
#  FUNCTIONS  #
# --- --- --- #

# check if the file is empty (returns TRUE if file is empty)
file.empty <- function(inputFile){file.size(inputFile) <= 1}

# Read table in only if the file exists and is not empty (otherwise return NA)
loadIfExists <- function(inputFile, sep='\n', ...){
  if(!is.na(inputFile)){
    if(file.exists(inputFile) && !file.empty(inputFile)){
      return(read.delim(inputFile, sep = sep, header = TRUE, ...))
      }else{
        return(NA)
        }
   } else{
   return(NA)
   }
  }

tree_mat <- function(z){
  b <- unlist(z)
  nm <- names(z)
  nx <- sapply(z, length)
  names(b) <- rep(nm, nx)
  c <- cbind(row.names(as.matrix(b)), as.matrix(b))
  row.names(c) <- NULL
  c <- apply(c, 2, as.integer)
  return(c)
}

get_clones <- function(populations){
  clone_df <- data.frame()
  for (clone in names(populations)){
    pop <- populations[[clone]]
    row <- data.frame(clone=clone, cellular_prevalence=pop$cellular_prevalence, num_cnvs=pop$num_cnvs, num_ssms=pop$num_ssms)
    clone_df <- rbind(clone_df, row)
    rm(pop, row)
  }
  return(clone_df)
}

annotate_disease_relevance <- function(association.disease_labels_truncated){
  ifelse(
    grepl('myeloma', association.disease_labels_truncated, ignore.case=TRUE), 2, # mutliple myeloma first
    ifelse(grepl('leuk|lymph|hemat|myel|b-cell|t-cell|Macroglobulinemia|sezary|polycythemia|macroglobulinemia|langerhan', association.disease_labels_truncated, ignore.case=TRUE), 1.5, # other hematological malignancies
    ifelse(grepl('onco|tumor|cancer|sarcoma|carcinoma', association.disease_labels_truncated, ignore.case=TRUE), 1, 0.75) # then cancer; things that are note even cancer related get slightly penalized.
    )
  )
}

cat("* Loading input data ...\n")
# ----

silent_muts <- c('Silent', 'Intron', "__UNKNOWN__", "COULD_NOT_DETERMINE", "3'UTR" , "5'UTR","IGR",  "5'Flank" , "RNA")

# Evidence and resonse based multiplier tables for VICC and CIVIC.
evidence2score <- c("A"=1, "B"=0.5, "C"=0.3)
evidence_score_df <- data.frame(association.evidence_level=names(evidence2score), evidence_score=evidence2score)
response_labels <- tolower(c("resistance","resistant", "sensitive", "Responsive", "Non Responsive", "Sensitivity/Response", "Sensitivity", "No Responsive", "no benefit", 'poor outcome'))
response_labels_sup <- tolower(c("sensitive", "Responsive", "Sensitivity/Response", "Sensitivity"))
response_label_df <- data.frame(association.response_type=response_labels, response.score_multiplier=ifelse(response_labels %in% response_labels_sup, 2, -1))
# --- 

# - Drug Tables - #
outdir <- opt$outdir
# All data that is processed and scored is appended to the pre-curated list of drugs we evaluate. (via a left_join)
# This pre-curated list contains annotations and evidence statements for various markers
drugMarkerList  <- read.delim(opt$daphniBiomarkerTable)

fullDrugList  <- read.delim(opt$daphniDrugTable)
fullDrugList$DrugNameSimple <- trimws(tolower(fullDrugList$DrugNameSimple)) # ensure cleanup

# This mutation census, based on oncokb primarily, gives us known gain of function and loss of function mutations in genes related to tier 1 biomarkers. We will use this to identify which druggable mutations are activating or deactivating.
mutationCensus <- read.delim(opt$daphniMutationCensus)

# Somatic Mutations
viccSOMTable <- loadIfExists(opt$somMutFileVICC, sep='\t', fill = TRUE, row.names=NULL)
civicSOMTable <- loadIfExists(opt$somMutFileCIVIC, sep='\t', fill = TRUE, row.names=NULL)
somTable <- loadIfExists(opt$somMutFile, sep='\t', fill=TRUE, comment.char='#', row.names=NULL)

# CNVs
viccCNATable <- loadIfExists(opt$cnaFileVICC, sep='\t', fill = TRUE, row.names=NULL)
civicCNATable <- loadIfExists(opt$cnaFileCIVIC, sep='\t', fill = TRUE, row.names=NULL)
cnaTable <- loadIfExists(opt$cnaFile, sep='\t', fill=TRUE, row.names=NULL)

# Expression
zscores <- loadIfExists(opt$zScoreTable, sep='\t')
selinescores <- loadIfExists(opt$selineScoresFile, sep='\t')

viccExpTable <- loadIfExists(opt$geneExprFileVICC, sep=',', fill = TRUE, row.names=NULL)
civicExpTable <- loadIfExists(opt$geneExprFileCIVIC, sep=',', fill = TRUE, row.names=NULL)

# genomic markers
scarData <- loadIfExists(opt$scarFile, sep="\t")

# mm PSN
mmPSNTable <- loadIfExists(opt$mmPSNFile, sep=',')

# mm hallmarks
mmHallMarksTable <- loadIfExists(opt$mmHallMarks, sep=',')

# sample ID
sampleID <- opt$sampleID

rnaSampleID <- gsub('-', '_', gsub('DNA', 'RNA', sampleID), fixed=TRUE)
dnaSampleID <- gsub('RNA', 'DNA', sampleID)

# set up list object for storing results
resultTableList <- list()

# --- --- --- --- --- --- --- #
cat("* --- PROCESSING DATA SOURCES --- *\n")

# --- --- --- --- --- --- ---  #
# General / Genomic Biomarkers #
# --- --- --- --- --- --- ---  #

# HRD / SCAR SCORE
if(!is.na(scarData) && length(scarData)>0 && dnaSampleID %in% scarData$sampleID){
  cat("* HRD-LOH data....\n")

  # find which drugs are predictable by SCAR score from our annotation table 
  scarDrugs <- fullDrugList$HRD.flag
  scarDrugs[is.na(scarDrugs)] <- 0
  scarDrugs <- fullDrugList$Rx.Bucket[scarDrugs > 0]

  # put scar scores into categories based on the quantiles and append score accordingly.
  scarData$riskCategory <- factor(cut(scarData$HRD.sum, breaks = quantile(scarData$HRD.sum), labels = c('low', 'moderate-low', 'moderate-high', 'high'), include.lowest=TRUE, ordered_result=TRUE))
  scarPatient <- scarData[scarData$sampleID == dnaSampleID, 'HRD.sum']
  scarRisk <- scarData[scarData$sampleID == dnaSampleID, 'riskCategory']

  scar_results <- data.frame(
    Rx.Bucket = scarDrugs,
    variant_type = 'GSS',
    variant_name = 'GSS',
    source = 'MSSM daphniDB',
    Variant_Classification = as.character(scarRisk),
    clone = NA,
    alteration_tier = 1,
    variant_statement = paste0("The patient's tumor has a scar score of ", as.numeric(scarPatient), ' which falls into the ', as.character(scarRisk), ' (', names(quantile(scarData$HRD.sum))[as.numeric(scarRisk)], ' quantile).'),
    variant_score = ifelse(scarRisk=="high", 1, 0) # make the score zero-indexed so that low scar score contributes nothing.
    )

  if(scar_results$variant_score != 0){
    resultTableList[['scarHRD']] <- scar_results %>% left_join(fullDrugList) %>% left_join(drugMarkerList) %>% unique()
  }

}


# SELINESCORES
if(!is.na(selinescores) && length(selinescores) > 0 && sum(grepl(rnaSampleID, rownames(selinescores))) > 0){
  cat("* Selinexor signature....\n")

  selinescores$category <- cut(selinescores$score, breaks=quantile(selinescores$score), labels = c('low', 'moderate-low', 'moderate-high', 'high'), include.lowest=TRUE, ordered_result=TRUE)
  selineScorePatient <- selinescores[grepl(rnaSampleID, rownames(selinescores)), 'score']
  selineCatPatient <- selinescores[grepl(rnaSampleID, rownames(selinescores)), 'category']
  
  seli_results <- data.frame(
    Rx.Bucket = 'XPO1 Inhibitor',
    DrugNameSimple = 'selinexor',
    source = 'MSSM daphniDB',
    Variant_Classification = as.character(selineCatPatient),
    variant_type = 'selinescore signature',
    variant_name = 'selinescore signature',
    clone = NA,
    alteration_tier = 1,
    variant_statement = paste("The patient's tumor has a selinexor score of", round(as.numeric(selineScorePatient), 2), 'which falls into the ',  as.character(selineCatPatient), '(', names(quantile(selinescores$score))[as.numeric(selineCatPatient)], 'quantile ).'),
    variant_score = ifelse(as.numeric(selineCatPatient) > 3, as.numeric(selineCatPatient), 0) # make the score zero-indexed so that selinescores below 50% quantile contribute nothing.
    )

  if(seli_results$variant_score != 0){
    resultTableList[['selinescore']] <- seli_results %>% left_join(fullDrugList) %>% left_join(drugMarkerList) %>% unique()
  }
}

# PSN Based Buckets: Venetoclax and Seli + Venetoclax:
if(!is.na(mmPSNTable) && length(mmPSNTable)>0){
  cat("* MM-PSN....\n")
  group <- mmPSNTable$Subgroup[1]
  mmPSNPredictedClass <- tolower(group)

  psn_results <- data.frame(
    source = 'MSSM daphniDB',
    Variant_Classification = tolower(as.character(mmPSNPredictedClass)),
    variant =  tolower(as.character(mmPSNPredictedClass)),
    variant_type = 'mm-psn',
    variant_name = paste0('MM-PSN Group ', tolower(as.character(mmPSNPredictedClass))),
    clone = NA,
    alteration_tier = 1,
    variant_statement = paste0("The predicted MM-PSN class for this patient's tumor is ", mmPSNPredictedClass),
    variant_score = 4
    )
  
  resultTableList[['mmpsn']] <- psn_results %>% left_join(drugMarkerList) %>% inner_join(fullDrugList)
}

# mm Hallmarks Based Buckets: Venetoclax and Seli + Venetoclax; FGFR, and other translocations...
if(!is.na(mmHallMarksTable) && length(mmHallMarksTable)>0 && nrow(mmHallMarksTable) > 0){
    cat("* mmHallMarks....\n")

    hallmark_results <- data.frame(
      source = 'MSSM daphniDB',
      Variant_Classification = mmHallMarksTable$variant_type,
      variant =  mmHallMarksTable$variant_name,
      variant_type = mmHallMarksTable$variant_type,
      variant_name = mmHallMarksTable$variant_name,
      clone = NA,
      alteration_tier = 1,
      variant_statement = paste(mmHallMarksTable$prevalence_statement, mmHallMarksTable$prognostic_statement, paste0('(PMID(s): ', mmHallMarksTable$pmid, ')')),
      variant_score = 4
      ) %>%
      left_join(drugMarkerList[is.na(drugMarkerList$Hugo_Symbol)|nchar(drugMarkerList$Hugo_Symbol) == 0,]) %>%
      inner_join(fullDrugList)
    
    resultTableList[['hallmarks']] <- hallmark_results
}


# --- --- --- --- --- --- --- --- #
#      RNA GENE EXPRESSION:       #
# --- --- --- --- --- --- --- --- #


# --- RNA Z-Scores --- --- --- ---

# INFO: Pre-process RNA data & annotate pathway table with expressed genes or drug targets
# load in the zscore table & use that; Assume this is already filtered for correct coverage
# +1 if the drug target is over or under expressed, and +2 if the up/down-regulated gene target is also a driver gene.
# Append related drug-alteration information to the drugPathwayList, and scale the tallied alterations downstream alongside the DNA pathway Alterations

# Additionally, grab expressed in-house custom markers.

resultTableList[['expression']]  <- data.frame()

if(!is.na(zscores) && length(zscores) > 0 && sum(grepl(rnaSampleID, colnames(zscores))) > 0){

  cat("* RNA Expression: Drug Targets & Pathways\n")
  
  zscoresPt <- as.numeric(zscores[, grepl(rnaSampleID, colnames(zscores))])
  names(zscoresPt) <- zscores$Gene
  expressedGenes <- zscoresPt[zscoresPt > 1.5 | zscoresPt < (-1.5)]

  if(length(expressedGenes) > 0){ 
    exprs_df <- data.frame(
        Zscore = expressedGenes,
        Hugo_Symbol = names(expressedGenes),
        Variant_Classification = ifelse(expressedGenes>0, 'over-expressed', 'under-expressed'),
        variant_effect = ifelse(expressedGenes>0, 'activating', 'deactivating'),
        variant = ifelse(expressedGenes>0, 'over-expressed', 'under-expressed'),
        variant_type = 'expression',
        disease_category = as.character('myeloma'),
        association.disease_labels_truncated = 'multiple myeloma',
        association.evidence_level = 'A',
        clone = NA,
        alteration_tier = 2,
        source = 'MSSM daphniDB'
        ) %>%
      inner_join(evidence_score_df) %>%
      mutate(
        disease.score_multiplier = annotate_disease_relevance(association.disease_labels_truncated),
        variant_score = evidence_score * disease.score_multiplier,
        variant_statement = paste0('The ', Hugo_Symbol, ' gene is ', Variant_Classification, ' in this tumor with a z-score of ', round(Zscore, 2), '.'),
        variant_name = paste(Hugo_Symbol, Variant_Classification)
        ) %>%
      filter(!is.na(Hugo_Symbol) & Hugo_Symbol != '')

    if(nrow(exprs_df) > 0){
      cat("* RNA Expression: Custom Markers\n")
      res <- drugMarkerList %>% inner_join(fullDrugList) %>% inner_join(exprs_df)

      if(nrow(res) > 0){
        
        res <- res %>%
          unique() %>% 
          inner_join(response_label_df) %>%
          mutate(
            variant_score = variant_score * response.score_multiplier
            ) %>%
          dplyr::select(
            Hugo_Symbol, Zscore, source, DrugNameSimple, Variant_Classification, Rx.Bucket,
            association.evidence_level, variant_statement, evidence_statement, variant_effect,
            association.response_type, disease_category, variant_name, variant_type,
            association.disease_labels_truncated, variant_score, clone, variant, Tier, alteration_tier
            ) %>%
          left_join(drugMarkerList)

        resultTableList[['expression']] <- tryCatch(
          {
            rbind(res, resultTableList[['expression']])
          },
          error = function(e){resultTableList[['expression']]}
        )
      }
    }
  }
}


# --- --- --- --- --- --- --- --- --- --- --- #
#     NOTE ON VICC & CIVIC SCORE CALCULATION  #
# --- --- --- --- --- --- --- --- --- --- --- #

# INFO: Each VICC or CIVIC entry is scored as a function of the entry's evidence level, the relevance of the associated disease, and whether the association is for response or resistance.
# Evidence levels for each entry are the starting score weight. This is further adjusted by response type and related disease relevance.
# Entries supportive of the drug are given a x2 multiplier and entries indicative of resistance are given a -1 multiplier, so that resistance counts *against* the final score.
# Additionally, entries whose association is in multiple myeloma are given another x2 multiplier, other hematological malignancies are given a x1.5 multiplier, and everything else is simply given x1.
# Therefore, each entry is given a score defined as: evidence_score * response_type_multiplier * disease_relevance_multiplier
# These scores are then summed for all entries associated with each given drug, and the resulting scores are appended to the fullDrugList table.
# In the case of SNVs and CNVs where clone information is available; scores are calculated on a per-drug, per-clone basis.

# --- VICC --- --- --- ---
if (!is.na(viccExpTable) && length(viccExpTable) > 0) {
  if(nrow(viccExpTable) > 0){
    cat("* RNA Expression: VICC ...\n")
    viccExpTable$Zscore <- viccExpTable[, grepl(sampleID, colnames(viccExpTable))]
    
    vicc_res <- viccExpTable %>%
      mutate(
        source = paste('VICC', source, sep='/'),
        association.disease_labels_truncated = tolower(association.disease_labels_truncated),
        association.response_type = tolower(association.response_type),
        DrugNameSimple = tolower(association.drug_labels),
        Variant_Classification = ifelse(Zscore>0, 'over-expressed', 'under-expressed'),
        variant_effect = ifelse(Zscore>0, 'activating', 'deactivating'),
        variant = ifelse(Zscore>0, 'over-expressed', 'under-expressed'),
        variant_type = 'expression',
        evidence_statement = association.description,
        Hugo_Symbol = Gene,
        clone = NA,
        alteration_tier = 2,
        variant_name = paste(Hugo_Symbol, Variant_Classification),
        variant_statement = paste0('The ', Hugo_Symbol, ' gene is ', Variant_Classification, ' in this tumor with a z-score of ', round(Zscore, 2), '.')
        ) %>% 
      inner_join(evidence_score_df) %>%
      inner_join(response_label_df) %>%
      mutate(
        association.response_type = as.character(ifelse(response.score_multiplier > 0, 'sensitivity', 'resistance')), 
        disease.score_multiplier = annotate_disease_relevance(association.disease_labels_truncated),
        variant_score = evidence_score * response.score_multiplier * disease.score_multiplier,
        disease_category = as.character(ifelse(disease.score_multiplier == 2, 'myeloma', ifelse(disease.score_multiplier == 1.5, 'hematological', ifelse(1, 'cancer', 'other'))))
        ) %>% unique() %>% 
      inner_join(fullDrugList) %>%
      dplyr::select(
        Hugo_Symbol, Zscore, source, DrugNameSimple, Variant_Classification, Rx.Bucket,
        association.evidence_level, variant_statement, evidence_statement, variant_effect,
        association.response_type, disease_category, variant_name,variant_type,
        association.disease_labels_truncated, variant_score, clone, variant, Tier, alteration_tier
        ) %>%
      left_join(drugMarkerList) 

    resultTableList[['expression']] <- tryCatch(
      {
      vicc_res %>% rbind(resultTableList[['expression']])
      },
      error = function(e){
      resultTableList[['expression']]
      })

    rm(vicc_res)
  }
} 

# --- CIVIC --- --- --- ---
if (!is.na(civicExpTable) && length(civicExpTable) >= length(cols)) {
  if(nrow(civicExpTable) > 0){
    cat("* RNA Expression: CIVIC ...\n")
    civicExpTable$Zscore <- civicExpTable[, grepl(sampleID, colnames(civicExpTable))]
    
    civic_res <- civicExpTable %>% 
      mutate(
        DrugNameSimple = tolower(drugs),
        source = 'CIVIC',
        association.response_type = tolower(civicExpTable$clinical_significance),
        association.disease_labels_truncated = tolower(civicExpTable$disease),
        association.evidence_level = civicExpTable$evidence_level,
        Hugo_Symbol = Gene,
        clone = NA,
        Variant_Classification = ifelse(Zscore>0, 'over-expressed', 'under-expressed'),
        variant_effect = ifelse(Zscore>0, 'activating', 'deactivating'),
        variant_type = 'expression',
        variant = ifelse(Zscore>0, 'over-expressed', 'under-expressed'),
        variant_name = paste(Hugo_Symbol, Variant_Classification),
        alteration_tier = 2,
        variant_statement = paste0('The ', Hugo_Symbol, ' gene is ', Variant_Classification, ' in this tumor with a z-score of ', round(Zscore, 2), '.')
      ) %>% 
      inner_join(evidence_score_df) %>% 
      inner_join(response_label_df) %>%
      mutate(
        association.response_type = as.character(ifelse(response.score_multiplier > 0, 'sensitivity', 'resistance')),
        disease.score_multiplier = annotate_disease_relevance(association.disease_labels_truncated),
        disease_category = as.character(ifelse(disease.score_multiplier == 2, 'myeloma', ifelse(disease.score_multiplier == 1.5, 'hematological', ifelse(1, 'cancer', 'other')))),
        variant_score = evidence_score * response.score_multiplier * disease.score_multiplier) %>% 
        unique() %>%
      inner_join(fullDrugList) %>% 
      dplyr::select(
        Hugo_Symbol, Zscore, source, DrugNameSimple, Variant_Classification, Rx.Bucket,
        association.evidence_level, variant_statement, evidence_statement, variant_type,
        association.response_type, disease_category, variant_name, variant_effect,
        association.disease_labels_truncated, variant_score, clone, variant, Tier, alteration_tier
        ) %>%
      left_join(drugMarkerList)

    resultTableList[['expression']] <- tryCatch({
        civic_res %>% rbind(resultTableList[['expression']])
      },
      error = function(e){
        resultTableList[['expression']]
      })

    rm(civic_res)
  }
} 

# --- --- --- --- --- #
#  Somatic Mutations: #
# --- --- --- --- --- #
resultTableList[['somatic_mutation']]  <- data.frame()

# Custom SOM markers
if(!is.na(somTable) && length(somTable)>0){
  if (nrow(somTable) > 0) {
  cat("* Somatic Mutations: Custom Markers ...\n")

  somTableFilt <- somTable %>%
    mutate(
      Start_Position = as.character(Start_Position),
      End_Position = as.character(End_Position),
      depth = as.numeric(t_alt_count)+as.numeric(t_ref_count),
      tumor_f = as.numeric(t_alt_count)/depth,
      vaf = tumor_f,
      protein_change = gsub('p.','', ifelse(Protein_Change!='', Protein_Change, 'unknown/not-applicable'), fixed=TRUE),
      cDNA_Change = gsub('c.', '', cDNA_Change, fixed=TRUE),
      source = 'MSSM daphniDB',
      variant_name = paste(Hugo_Symbol, ifelse(!is.na(protein_change), protein_change, cDNA_Change)),
      variant_type = 'somatic mutation',
      association.disease_labels_truncated = as.character('multiple myeloma'),
      association.evidence_level = 'A',
      clone = NA,
      alteration_tier = 1
      ) %>%
    inner_join(evidence_score_df) %>%
    mutate(
        disease.score_multiplier = annotate_disease_relevance(association.disease_labels_truncated),
        disease_category = as.character(ifelse(disease.score_multiplier == 2, 'myeloma', ifelse(disease.score_multiplier == 1.5, 'hematological', ifelse(1, 'cancer', 'other'))))
      ) %>%
    unique() %>%
    inner_join(mutationCensus)

  if(nrow(somTableFilt)>0){

    somTableFilt <- drugMarkerList %>%
      inner_join(fullDrugList) %>%
      right_join(somTableFilt) %>%
      unique() %>%
      inner_join(response_label_df) %>%
      mutate(
        variant_score = evidence_score * response.score_multiplier,
        variant_statement = paste0('This somatic ', tolower(gsub('_', ' ', Variant_Classification)) ,' was found in the ', Hugo_Symbol, ' gene with a variant allele frequency (VAF) of ', round(tumor_f*100, 2), '%, supported by ', t_alt_count, ' alternate reads in the tumor DNA. It is an ', variant_effect, ' alteration, with a functional protein change of ', protein_change, ', caused by an underlying cDNA change of ', cDNA_Change, '.')
        ) %>%
      dplyr::select(
          Hugo_Symbol, Chromosome, Start_Position, End_Position, variant_name, Rx.Bucket,
          Variant_Classification, cDNA_Change, protein_change, clone, variant_type,
          tumor_f, t_alt_count, t_ref_count, source, DrugNameSimple,  variant_effect,
          association.evidence_level,variant_statement, evidence_statement, association.response_type, 
          disease_category, association.disease_labels_truncated, variant_score, Tier, alteration_tier
          ) %>% 
      left_join(drugMarkerList)
      
    resultTableList[['somatic_mutation']] <-  tryCatch({
        somTableFilt %>% rbind(resultTableList[['somatic_mutation']])
        },
      error = function(e){resultTableList[['somatic_mutation']]})
    rm(somTableFilt)
    }
  }
}
# --- VICC --- --- --- ---
if (!is.na(viccSOMTable) && length(viccSOMTable)>0) {
  if(nrow(viccSOMTable) > 0){
    cat("* Somatic Mutations: VICC ...\n")

    vicc_res <- tryCatch({viccSOMTable %>%
      mutate(
        source = paste('VICC', source, sep='/'),
        Start_Position = as.character(Start_Position),
        End_Position = as.character(End_Position),
        tumor_f = as.numeric(tumor_f),
        t_alt_count = as.numeric(t_alt_count),
        t_ref_count = as.numeric(t_ref_count),
        association.disease_labels_truncated = tolower(association.disease_labels_truncated),
        association.response_type = tolower(association.response_type),
        DrugNameSimple = tolower(association.drug_labels),
        protein_change = gsub('p.','', ifelse(Protein_Change_x!='', Protein_Change_x, 'unknown/not-applicable'), fixed=TRUE),
        cDNA_Change = gsub('c.', '', cDNA_Change, fixed=TRUE),
        evidence_statement = association.description,
        variant_type = 'somatic mutation',
        variant_statement = paste0('This somatic ', tolower(gsub('_', ' ', Variant_Classification)) ,' was found in the ', Hugo_Symbol, ' gene with a variant allele frequency (VAF) of ', round(tumor_f*100, 2), '%, supported by ', t_alt_count, ' alternate reads in the tumor DNA. Its predicted change to the protein is ', protein_change, ' and its predicted change to the cDNA is ', cDNA_Change, '.'),
        variant_name = paste(Hugo_Symbol, ifelse(!is.na(protein_change), protein_change, cDNA_Change)),
        clone = as.character(clone),
        alteration_tier = 2
      ) %>%
      inner_join(evidence_score_df) %>% 
      inner_join(response_label_df) %>%
      mutate(
        association.response_type = ifelse(response.score_multiplier > 0, 'sensitivity', 'resistance'),
        disease.score_multiplier = annotate_disease_relevance(association.disease_labels_truncated),
        disease_category = as.character(ifelse(disease.score_multiplier == 2, 'myeloma', ifelse(disease.score_multiplier == 1.5, 'hematological', ifelse(1, 'cancer', 'other')))),
        variant_score = evidence_score * response.score_multiplier * disease.score_multiplier
        ) %>% 
      unique() %>%
      inner_join(fullDrugList) %>%
      dplyr::select(
        Hugo_Symbol, Chromosome, Start_Position, End_Position, variant_name,
        Variant_Classification, cDNA_Change, protein_change, clone, variant_type,
        tumor_f, t_alt_count, t_ref_count, source, DrugNameSimple, Rx.Bucket,
        association.evidence_level,variant_statement, evidence_statement,
        association.response_type, disease_category, Tier,
        association.disease_labels_truncated, variant_score, alteration_tier
        ) %>%
      left_join(drugMarkerList)},
      error = function(e){data.frame()})

    resultTableList[['somatic_mutation']] <- tryCatch({
      vicc_res %>% rbind(resultTableList[['somatic_mutation']])
      }, 
      error = function(e){
        resultTableList[['somatic_mutation']]
        })
    rm(vicc_res)
  }
}

# --- CIVIC --- --- --- ---
if (!is.na(civicSOMTable) && length(civicSOMTable)>0) {
  if(nrow(civicSOMTable) > 0){
    
    cat("* Somatic Mutations: CIVIC ...\n")

    civic_res <- tryCatch({civicSOMTable %>%
      mutate(
        Start_Position = as.character(Start_Position),
        End_Position = as.character(End_Position),
        tumor_f = as.numeric(tumor_f),
        t_alt_count = as.numeric(t_alt_count),
        t_ref_count = as.numeric(t_ref_count),
        association.evidence_level = evidence_level,
        association.disease_labels_truncated = tolower(disease),
        association.response_type = tolower(clinical_significance),
        DrugNameSimple = tolower(drugs),
        protein_change = gsub('p.','', ifelse(Protein_Change_x!='', Protein_Change_x, 'unknown/not-applicable'), fixed=TRUE),
        cDNA_Change = gsub('c.', '', cDNA_Change, fixed=TRUE),
        source = 'CIVIC',
        alteration_tier = 2,
        variant_type = 'somatic mutation',
        variant_statement = paste0('This somatic ', tolower(gsub('_', ' ', Variant_Classification)) ,' was found in the ', Hugo_Symbol, ' gene with a variant allele frequency (VAF) of ', round(tumor_f*100, 2), '%, supported by ', t_alt_count, ' alternate reads in the tumor DNA. Its predicted change to the protein is ', protein_change, ' and its predicted change to the cDNA is ', cDNA_Change, '.'),
        variant_name = paste(Hugo_Symbol, ifelse(!is.na(protein_change), protein_change, cDNA_Change)),
        clone = as.character(clone)
      ) %>%
      inner_join(evidence_score_df) %>%
      inner_join(response_label_df) %>%
      mutate(
        association.response_type = ifelse(response.score_multiplier > 0, 'sensitivity', 'resistance'),
        disease.score_multiplier = annotate_disease_relevance(association.disease_labels_truncated),
        disease_category = as.character(ifelse(disease.score_multiplier == 2, 'myeloma', ifelse(disease.score_multiplier == 1.5, 'hematological', ifelse(1, 'cancer', 'other')))),
        variant_score = evidence_score * response.score_multiplier * disease.score_multiplier
        ) %>%
      unique() %>%
      inner_join(fullDrugList) %>%
      dplyr::select(
        Hugo_Symbol, Chromosome, Start_Position, End_Position, variant_name,
        Variant_Classification, cDNA_Change, protein_change, clone, variant_type,
        tumor_f, t_alt_count, t_ref_count, source, DrugNameSimple, Rx.Bucket,
        association.evidence_level, variant_statement, evidence_statement, association.response_type, 
        disease_category, association.disease_labels_truncated, variant_score, Tier, alteration_tier
        ) %>%
      left_join(drugMarkerList)},
      error = function(e){data.frame()})

    resultTableList[['somatic_mutation']] <- tryCatch({
      civic_res %>% rbind(resultTableList[['somatic_mutation']])
      }, 
      error = function(e){
        resultTableList[['somatic_mutation']]
        })
    rm(civic_res)
  }
}

# --- --- --- --- #
#   COPY NUMBER   #
# --- --- --- --- #

resultTableList[['cna']] <- data.frame()

# Custom CNV markers
if(!is.na(cnaTable) && length(cnaTable)>0){
  if(nrow(cnaTable) > 0){
    cat("* CNA: Custom Markers ...\n")
    # filter germline SNPs for known clinvar associations, adequate depth, and non-silent functional effect.
    cnaTableFilt <- cnaTable %>%
      mutate(
        source = 'MSSM daphniDB',
        Chromosome = ifelse(startsWith(as.character(chr), 'chr'),chr, paste0('chr', as.character(chr))),
        Start_Position = as.character(start),
        End_Position = as.character(stop),
        Hugo_Symbol = as.character(geneName),
        copynumber = as.numeric(copynumber),
        Variant_Classification = ifelse(copynumber > 2, 'amplification', 'deletion'),
        variant = ifelse(copynumber > 2, 'amplification', 'deletion'),
        variant_type = 'CNA',
        variant_effect = ifelse(copynumber > 2, 'activating', 'deactivating'),
        tumor_f = as.numeric(cell_prev),
        association.evidence_level="A",
        variant_name = paste(Hugo_Symbol, Variant_Classification),
        variant_statement = paste0('The ', Hugo_Symbol, ' gene has a copy number ', Variant_Classification, ' in this tumor with an estimated ', copynumber, ' copies and an estimated cell fraction of ', as.character(round(as.numeric(tumor_f)*100, 2)), '%.'),
        association.disease_labels_truncated='multiple myeloma',
        alteration_tier = 1,
        clone=NA
      ) %>%
      inner_join(evidence_score_df) %>%
      mutate(
        disease.score_multiplier = annotate_disease_relevance(association.disease_labels_truncated),
        disease_category = as.character(ifelse(disease.score_multiplier == 2, 'myeloma', ifelse(disease.score_multiplier == 1.5, 'hematological', ifelse(1, 'cancer', 'other')))),
        variant_score = evidence_score * disease.score_multiplier
        ) %>%
      filter(tumor_f >= 0.15 & copynumber != 2) %>%
      inner_join(drugMarkerList)

    if(nrow(cnaTableFilt)>0){
      cnaTableFilt <- cnaTableFilt %>%
        inner_join(fullDrugList) %>%
        inner_join(response_label_df) %>%
        mutate(variant_score = evidence_score * response.score_multiplier) %>%
        ungroup() %>% data.frame() %>% unique() %>%
        dplyr::select(
          Hugo_Symbol, Chromosome, Start_Position, End_Position, variant_name,
          Variant_Classification, band, clone, tumor_f, source, variant_statement,
          DrugNameSimple, association.evidence_level, evidence_statement, Rx.Bucket,
          association.response_type, disease_category, variant_effect, variant_type,
          association.disease_labels_truncated, variant_score, Tier, variant, alteration_tier
          ) %>% 
        left_join(drugMarkerList)
        
        resultTableList[['cna']] <-  tryCatch({
          cnaTableFilt %>% rbind(resultTableList[['cna']])
          }, 
          error = function(e){
            resultTableList[['cna']]
            })
        rm(cnaTableFilt)
      }
  }
}

# --- VICC --- --- --- ---
if(!is.na(viccCNATable) && length(viccCNATable)>0) {
  cat("* CNA: VICC ...\n")

  vicc_res <- tryCatch({viccCNATable %>%
    mutate(
      source = paste('VICC', source, sep='/'),
      association.disease_labels_truncated = tolower(association.disease_labels_truncated),
      DrugNameSimple = tolower(association.drug_labels),
      association.response_type = tolower(association.response_type),
      evidence_statement = association.description,
      clone = as.character(clones),
      Chromosome = ifelse(startsWith(as.character(chr), 'chr'),chr, paste0('chr', as.character(chr))),
      Start_Position = as.character(start),
      End_Position = as.character(stop),
      Hugo_Symbol = as.character(geneName),
      Variant_Classification = ifelse(as.numeric(copynumber) > 2, 'amplification', 'deletion'),
      variant_type = 'CNA',
      variant_effect = ifelse(copynumber > 2, 'activating', 'deactivating'),
      tumor_f = as.numeric(cell_prev),
      variant_name = paste(Hugo_Symbol, Variant_Classification),
      alteration_tier = 2,
      variant = ifelse(copynumber > 2, 'amplification', 'deletion'),
      variant_statement = paste0('The ', Hugo_Symbol, ' gene has a copy number ', Variant_Classification, ' in this tumor with an estimated ', as.character(copynumber), ' copies and an estimated cell fraction of ', as.character(round(as.numeric(tumor_f)*100, 2)), '%.')
    ) %>%
    inner_join(evidence_score_df) %>%
    inner_join(response_label_df) %>%
    mutate(
      association.response_type = ifelse(response.score_multiplier > 0, 'sensitivity', 'resistance'),
      disease.score_multiplier = annotate_disease_relevance(association.disease_labels_truncated),
      variant_score = evidence_score * response.score_multiplier * disease.score_multiplier,
      disease_category = as.character(ifelse(disease.score_multiplier == 2, 'myeloma', ifelse(disease.score_multiplier == 1.5, 'hematological', ifelse(1, 'cancer', 'other'))))
      ) %>%
    ungroup() %>% data.frame() %>% unique() %>%
    inner_join(fullDrugList) %>%
    dplyr::select(
      Hugo_Symbol, Chromosome, Start_Position, End_Position, variant_name,
      Variant_Classification, band, clone, tumor_f, source, variant_statement,
      DrugNameSimple, association.evidence_level, evidence_statement, Rx.Bucket,
      association.response_type, disease_category, variant_effect, variant_type,
      association.disease_labels_truncated, variant_score, variant, Tier, alteration_tier
     ) %>% 
     left_join(drugMarkerList)},
     error = function(e){data.frame()})

  resultTableList[['cna']] <- tryCatch({
    vicc_res %>% rbind(resultTableList[['cna']])
    }, 
    error = function(e){
      resultTableList[['cna']]
      })
  rm(vicc_res)
}

# --- CIVIC --- --- --- ---
if (!is.na(civicCNATable) && length(civicCNATable)>0) {
  if(nrow(civicCNATable) > 0){
    cat("* CNA: CIVIC ...\n")

    civic_res <- tryCatch({civicCNATable %>% 
      mutate(
        association.evidence_level = evidence_level,
        association.disease_labels_truncated = tolower(as.character(disease)),
        association.response_type = tolower(clinical_significance),
        DrugNameSimple = tolower(drugs),
        clone = as.character(clones),
        Chromosome = ifelse(startsWith(as.character(chr), 'chr'), chr, paste0('chr', as.character(chr))),
        Start_Position = as.character(start),
        End_Position = as.character(stop),
        Hugo_Symbol = geneName,
        Variant_Classification = ifelse(copynumber > 2, 'amplification', 'deletion'),
        variant_type = 'CNA',
        variant_effect = ifelse(copynumber > 2, 'activating', 'deactivating'),
        tumor_f = as.numeric(cell_prev),
        alteration_tier = 2,
        variant_name = paste(Hugo_Symbol, Variant_Classification),
        variant = ifelse(copynumber > 2, 'amplification', 'deletion'),
        variant_statement = paste0('The ', Hugo_Symbol, ' gene has a copy number ', Variant_Classification, ' in this tumor with an estimated ', copynumber, ' copies and an estimated cell fraction of ', as.character(round(as.numeric(tumor_f)*100, 2)), '%.'),
        source = 'CIVIC'
      ) %>%
      inner_join(evidence_score_df) %>%
      inner_join(response_label_df) %>%
      mutate(
        association.response_type = ifelse(response.score_multiplier > 0, 'sensitivity', 'resistance'),
        disease.score_multiplier = annotate_disease_relevance(association.disease_labels_truncated),
        variant_score = evidence_score * response.score_multiplier * disease.score_multiplier,
        disease_category = as.character(ifelse(disease.score_multiplier == 2, 'myeloma', ifelse(disease.score_multiplier == 1.5, 'hematological', ifelse(1, 'cancer', 'other'))))
        ) %>%
      ungroup() %>% data.frame() %>% unique() %>%
      inner_join(fullDrugList) %>%
      dplyr::select(
      Hugo_Symbol, Chromosome, Start_Position, End_Position, variant_statement,
      Variant_Classification, band, clone, tumor_f, source, variant_name, variant,
      DrugNameSimple, association.evidence_level, evidence_statement, Rx.Bucket,
      association.response_type, disease_category, variant_effect, variant_type,
      association.disease_labels_truncated, variant_score, Tier, alteration_tier
      ) %>% unique() %>%
      left_join(drugMarkerList)},
      error = function(e){data.frame()})

    resultTableList[['cna']] <- tryCatch({
      civic_res %>% rbind(resultTableList[['cna']])
      },
      error = function(e){
        resultTableList[['cna']]
        })
    rm(civic_res)
  }
}

# --- --- --- --- --- --- --- --- #
#     CONSOLIDATE ALL RESULTS     #
# --- --- --- --- --- --- --- --- #
cat("* --- CONSOLIDATING RESULTS --- *\n")

cat("* Collating variant-drug associations across data types\n")

# get a big table with all the variants supporting the ranked drugs.

# ensure all entries in the results list are dataframes
# and then remove empty dataframes from the list
resultTableList <- lapply(resultTableList, as.data.frame)
resultTableList <- resultTableList[sapply(resultTableList, function(x){nrow(x)>0})]

# then bind it all into one enormous dataframe
variantResultDetails <- unique(bind_rows(resultTableList))

# ERROR HANDLING: case where no variants are found!
# No columns in variantResultDetails means no actionable variants AT ALL were found.
# This is expected to be a rare scenario, but in this case, save empty result files and quit.
if (ncol(variantResultDetails) == 0){
  cat("STOPPING SCRIPT: there are no actionable variants in this sample!\n")

  variantTable <- data.frame()
  variantSummaryTable <- data.frame()
  rankedDrugsTable <- data.frame()

  cat(paste('* Saving supporting variant summary (All Tiers):', paste(outdir, 'actionable_variants.prediction_engine.results.tsv', sep='/'), '\n'))
  write_tsv(variantTable, file=paste(outdir, 'actionable_variants.prediction_engine.results.tsv', sep='/'))

  cat(paste('* Saving supporting variant-drug associations (All Tiers):', paste(outdir, 'variant_associations.prediction_engine.results.tsv', sep='/'), '\n'))
  write_tsv(variantSummaryTable, file=paste(outdir, 'variant_associations.prediction_engine.results.tsv', sep='/'))

  cat(paste('* Saving ranked drug details table (All Tiers):', paste(outdir, 'drug_recommendations.prediction_engine.results.tsv', sep='/'), '\n'))
  write_tsv(rankedDrugsTable, file=paste(outdir, 'drug_recommendations.prediction_engine.results.tsv', sep='/'))

  vars <- c("expression", "cna", "somatic_mutation")

  # save per-tool detailed results
  for (name in vars){
      filename <- paste(outdir, paste0(name, ".prediction_engine.results.tsv"), sep='/')
      cat(paste('* Saving:', filename, '\n'))
      write_tsv(data.frame(), filename)
      }

  quit(status=0, save='no')
}

variantTable <- unique(variantResultDetails[, c("variant_name", "clone")])
variantTable$clone <- as.character(variantTable$clone)

# Clonal Analysis 
if (!is.na(opt$treeFile) && file.exists(opt$treeFile) && !file.empty(opt$treeFile) && 'clone' %in% colnames(variantResultDetails)) {
  
  cat("* Annotating clones\n")
  # get top clonal tree from structure file
  # based on the llh value closest to zero
  structFile <- RJSONIO::fromJSON(opt$treeFile)
  llh <- sapply(structFile$trees, function(x){abs(x$llh)})
  top_model_name <- strsplit2(names(llh[llh == min(llh)]), split = "[.]")[1]
  cloneTree <- structFile$trees[[top_model_name]]
  
  # if this fails; return a matrix of just the founder clone...
  treeMatrix <- tryCatch({as.matrix(tree_mat(cloneTree$structure))}, error = function(e){matrix(c(0,0),nrow=1, ncol=2)})
  
  # Error Handling: for malformed tree matrices
  if(ncol(treeMatrix) < 2){
    # case where there is one founder clone and main clone
    if(nrow(treeMatrix) == 2 && ncol(treeMatrix) == 1){
      treeMatrix <- t(treeMatrix)
    }else{
      # if the matrix is somehow malformed; just show the founder clone...
      treeMatrix <- matrix(c(0,0), nrow=1, ncol=2)
    }
  }
  
  rm(structFile, llh, top_model_name) # clean-up

  # From the tree structure matrix, annotate a clone table
  # with tree information, and clone specific ccf, cnv burden, snv burden,
  # and each clone's downstream descendants
  cat("* Perfoming clone-aware analysis\n")

  cloneData <- tryCatch({
    data.frame(
      from = as.character(treeMatrix[,1]),
      to = as.character(treeMatrix[,2])
      ) %>%
    mutate(clone=from) %>% 
    group_by(clone) %>%
    summarise(descendant_clones = paste(unique(to), collapse='|')) %>%
    right_join(get_clones(cloneTree$populations)) %>%
    mutate(
      tree.linearity_index = cloneTree$linearity_index,
      tree.llh = cloneTree$llh,
      tree.clustering_index = cloneTree$clustering_index,
      tree.branching_index = cloneTree$branching_index
    )
    },
    error = function(e){
      data.frame(clone = NA, tree.llh = NA, tree.clustering_index = NA, tree.branching_index = NA, descendant_clones = NA)
      }
    )

  if(nrow(variantTable) > 0 && nrow(cloneData) > 0){
    cloneData <- cloneData %>% left_join(variantTable)
    cat("* Performing clone-aware drug score adjustments\n")
    
    # Apply clonal adjustment so that per-clone drug scores are adjusted for resistance in descendant clones
    variantResultDetails <- variantResultDetails %>%
      select(-clone) %>%
      left_join(cloneData)
  }
  rm(cloneTree)
}

cat('* Implementing Conditional Reccommendations\n')

# 1. BCL2 + XPO1
if(('chr1q gain' %in% variantResultDetails$variant_name) && (! "venetoclax" %in% variantResultDetails$DrugNameSimple) ){
variantResultDetails <- variantResultDetails %>%
  filter(Rx.Bucket != 'XPO1 + BCL2 Inhibitior Combination')
}

cat("* Summarizing drug scores across data types\n")
variantResultDetails$alteration_tier <- ifelse(variantResultDetails$source == "MSSM daphniDB", 1, 2)
variantTable <- unique(variantResultDetails[, c("variant_name", "variant_statement", "variant_type", "variant_effect", "Hugo_Symbol", "Variant_Classification", "clone", 'alteration_tier')])

if(! "descendant_clones" %in% colnames(variantResultDetails)){
  variantResultDetails$descendant_clones <- NA
}

# Negative scores indicate possible resistance; greater scores indicate possible benefit. 
# Values at or near zero indicate that there is either not enough evidence to support use or avoidance of the drug, 
# or that the evidence available is inconclusive with mixed response/resistance data.
rankedDrugsTable <- variantResultDetails %>%
  group_by(DrugNameSimple, Rx.Bucket, alteration_tier) %>%
  summarise(
    drug.summary_score = sum(na.omit(variant_score)), # adjust for combinations
    variant_names = paste(unique(variant_name), collapse='|')
    ) %>%
  ungroup() %>%
  arrange(-drug.summary_score) %>%
  inner_join(fullDrugList) 

rankedDrugsTable$drug.summary_score[is.na(rankedDrugsTable$drug.summary_score)] <- 0
rankedDrugsTable <- rankedDrugsTable[rankedDrugsTable$drug.summary_score != 0, ]
rankedDrugsTable$Tier[is.na(rankedDrugsTable$Tier)] <- '2'

cat("* Getting per-bucket rankings and results\n")
drugBucketSummaryResults <- rankedDrugsTable %>%
  inner_join(fullDrugList) 

drugBucketSummaryResults$drug.summary_score[is.na(drugBucketSummaryResults$drug.summary_score)] <- 0
drugBucketSummaryResults <- drugBucketSummaryResults %>%
  left_join(variantResultDetails) %>%
  group_by(Rx.Bucket, Tier, alteration_tier) %>%
  summarise(
    drug.count = n_distinct(DrugNameSimple),
    bucket.summary_score = sum(drug.summary_score)/drug.count, # overall summary of support for the bucket; summarised for evidence across all drugs in a given bucket
    affected_clones = gsub('|', ', ', paste(unique(na.omit(c(clone, descendant_clones))), collapse='|'), fixed=TRUE)
  ) %>% 
  arrange(-bucket.summary_score) %>% 
  ungroup() %>%
  dplyr::select(-drug.count)

rankedDrugsTable <- rankedDrugsTable %>% inner_join(drugBucketSummaryResults)

cat("* Updating variant-drug assocations with per-drug summary scores\n")
variantSummaryTable <- variantResultDetails %>%
  left_join(rankedDrugsTable) %>%
  group_by(Tier, Rx.Bucket, DrugNameSimple, variant_statement, variant_type, variant_effect, association.response_type, variant_name, Hugo_Symbol, Variant_Classification, affected_clones, alteration_tier) %>%
  summarise(
    evidence_statement.combined = paste(unique(evidence_statement[!is.na(evidence_statement) & evidence_statement!='']), collapse=' '),
    sources_statement.combined = paste('The assocation between this alteration and', DrugNameSimple, 'is supported in this report by data compiled from the', knitr::combine_words(unique(na.omit(source))), 'database(s).')
    ) %>%
  filter(!is.na(variant_statement)) %>% ungroup() %>% unique() %>%
  inner_join(variantResultDetails)

cat("* --- SAVING RESULTS --- *\n")
cat(paste('* Saving output results to output directory:', outdir, '\n'))

# save per-tool detailed results
for (name in names(resultTableList)){
  if(nrow(resultTableList[[name]]) > 0) {
    filename <- paste(outdir, paste0(name, ".prediction_engine.results.tsv"), sep='/')
    cat(paste('* Saving:', filename, '\n'))
    write_tsv(resultTableList[[name]], filename)
    rm(filename, name)
    }
}

cat(paste('* Saving supporting variant summary (All Tiers):', paste(outdir, 'actionable_variants.prediction_engine.results.tsv', sep='/'), '\n'))
write_tsv(variantTable, file=paste(outdir, 'actionable_variants.prediction_engine.results.tsv', sep='/'))

cat(paste('* Saving supporting variant-drug associations (All Tiers):', paste(outdir, 'variant_associations.prediction_engine.results.tsv', sep='/'), '\n'))
write_tsv(variantSummaryTable, file=paste(outdir, 'variant_associations.prediction_engine.results.tsv', sep='/'))

cat(paste('* Saving ranked drug details table (All Tiers):', paste(outdir, 'drug_recommendations.prediction_engine.results.tsv', sep='/'), '\n'))
write_tsv(rankedDrugsTable, file=paste(outdir, 'drug_recommendations.prediction_engine.results.tsv', sep='/'))

## prognostic_marker_results <- prognosticMarkerList[prognostic_marker_filt, ]
## if(nrow(prognostic_marker_results)>0){
##   cat(paste('* Saving prognostic biomarker findings:', paste(outdir, 'prognostic_markers.prediction_engine.results.tsv', sep='/'), '\n'))
##   write_tsv(prognostic_marker_results, file=paste(outdir, 'prognostic_markers.prediction_engine.results.tsv', sep='/'))
## }  

cat("* --- DONE! ---  *\n")