# temporary list for testing
opt <- list(
  outdir='data/outputs/MM_4289_BM_01',
  sampleID='MM_4289_BM_01',
  # reference tables
  daphniDrugTable='refData/daphniDrugTableSourcesAnnotated.tsv',
  daphniMutationCensus='refData/daphni_mutation_census.tsv',
  daphniBiomarkerTable='refData/daphni_tier1_biomarkers.tsv',
  # expression
  zScoreTable='data/merge_results/zscores_latest.csv',  
  geneExprFileVICC="data/MM_4289_BM_01/vicc_expression_data.csv",
  geneExprFileCIVIC="data/MM_4289_BM_01/civic_expression_data.csv",
  selineScoresFile='data/MM_4289_BM_01/selinescores_latest.tsv',
  # clonality
  treeFile="data/MM_4289_BM_01/sample.summ.json.gz",
  # psn
  mmPSNFile='data/MM_4289_BM_01/Predicted_class.csv',
  # snv
  somMutFileCIVIC="data/MM_4289_BM_01/civic_mutation_clone_data.txt",
  somMutFileVICC="data/MM_4289_BM_01/vicc_mutation_clone_data.txt",
  somMutFile='data/MM_4289_BM_01/annotated_variants.consensus.vcf',
  # cnv
  cnaFileCIVIC="data/MM_4289_BM_01/civic_cnv_clone_data.txt",
  cnaFileVICC="data/MM_4289_BM_01/vicc_cnv_clone_data.txt",
  cnaFile='data/MM_4289_BM_01/output_CNV_file.txt',
  # scar / genomic instability
  scarFile="data/merge_results/scarScore_latest.tsv",
  # hallmarks
  mmHallMarksFile='data/MM_4289_BM_01/mm_hallmarks_results.csv'
)