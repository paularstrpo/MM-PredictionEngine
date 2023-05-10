#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
  DockerRequirement:
    dockerPull: "sinaiiidgst/predengine:954ca51"

inputs:
  outdir:
    type: string
    default: $(runtime.outdir)
    inputBinding:
      valueFrom: $(runtime.outdir)
      prefix: "--outdir"

  sampleID:
    type: string
    inputBinding:
      prefix: "--sampleID"

  daphniDrugTable:
    type: File
    inputBinding:
      prefix: "--daphniDrugTable"

  daphniMutationCensus:
    type: File
    inputBinding:
      prefix: "--daphniMutationCensus"

  daphniBiomarkerTable:
    type: File
    inputBinding:
      prefix: "--daphniBiomarkerTable"

  somMutFileVICC:
    type: File?
    inputBinding:
      prefix: "--somMutFileVICC"

  somMutCIVIC:
    type: File?
    inputBinding:
      prefix: "--somMutFileCIVIC"

  somMutFile:
    type: File?
    inputBinding:
      prefix: "--somMutFile"

  cnaFileVICC:
    type: File?
    inputBinding:
      prefix: "--cnaFileVICC"

  cnaFileCIVIC:
    type: File?
    inputBinding:
      prefix: "--cnaFileCIVIC"

  cnaFile:
    type: File?
    inputBinding:
      prefix: "--cnaFile"

  geneExprFileVICC:
    type: File?
    inputBinding:
      prefix: "--geneExprFileVICC"

  geneExprFileCIVIC:
    type: File?
    inputBinding:
      prefix: "--geneExprFileCIVIC"

  scarFile:
    type: File?
    inputBinding:
      prefix: "--scarFile"

  treeFile:
    type: File?
    inputBinding:
      prefix: "--treeFile"

  zScoreTable:
    type: File?
    inputBinding:
      prefix: "--zScoreTable"

  selineScoresFile:
    type: File?
    inputBinding:
      prefix: "--selineScoresFile"

  mmPSNFile:
    type: File?
    inputBinding:
      prefix: "--mmPSNFile"

  mmHallMarks:
    type: File?
    inputBinding:
      prefix: "--mmHallMarks"

baseCommand: [Rscript, /bin/predEngine.R]

outputs:
  output_preds:
    type:
        - "null"
        - File[]
    outputBinding:
      glob: "*.prediction_engine.results.tsv"

  all_drug_details_ranked:
    type: ["null", File]
    outputBinding:
      glob: "drug_recommendations.prediction_engine.results.tsv"

  variant_details:
    type: ["null", File]
    outputBinding:
      glob: "actionable_variants.prediction_engine.results.tsv"

  variant_summary_ranked:
    type: ["null", File]
    outputBinding:
      glob: "variant_associations.prediction_engine.results.tsv"

  somatic_mutation:
    type: ["null", File]
    outputBinding:
      glob: "somatic_mutation.prediction_engine.results.tsv"

  cna:
    type: ["null", File]
    outputBinding:
      glob: "cna.prediction_engine.results.tsv"

  expression:
    type: ["null", File]
    outputBinding:
      glob: "expression.prediction_engine.results.tsv"
