SAMPLEID="MM_4289_BM_01"
docker build -t sinaiiidgst/predengine:test .
docker run -v ${PWD}/data/:/data/ -v ${PWD}/refData/:/refData/ sinaiiidgst/predengine:test \
Rscript /bin/predEngine.R \
--outdir /data/$SAMPLEID \
--daphniBiomarkerTable /refData/daphni_tier1_biomarkers.tsv \
--daphniDrugTable /refData/daphniDrugTableSourcesAnnotated.tsv \
--daphniMutationCensus /refData/daphni_mutation_census.tsv \
--sampleID $SAMPLEID \
--geneExprFileCIVIC /data/$SAMPLEID/civic_expression_data.csv \
--geneExprFileVICC /data/$SAMPLEID/vicc_expression_data.csv \
--cnaFile /data/$SAMPLEID/output_CNV_file.txt \
--cnaFileVICC /data/$SAMPLEID/vicc_cnv_clone_data.txt \
--selineScoresFile /data/$SAMPLEID/selinescores_latest.tsv \
--somMutFileCIVIC /data/$SAMPLEID/civic_mutation_clone_data.txt \
--somMutFile /data/$SAMPLEID/annotated_variants.consensus.vcf \
--somMutFileVICC /data/$SAMPLEID/vicc_mutation_clone_data.txt \
--scarFile /data/merge_results/scarScore_latest.tsv \
--treeFile /data/$SAMPLEID/sample.summ.json.gz \
--mmPSNFile /data/$SAMPLEID/Predicted_class.csv \
--mmHallMarks /data/$SAMPLEID/mm_hallmarks_results.csv \
--cnaFileCIVIC /data/$SAMPLEID/civic_cnv_clone_data.txt \
--zScoreTable /data/merge_results/zscores_latest.csv
