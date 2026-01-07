# Tutorial Workflow
See ALTER-code/5_tutorial-workflows/2_MARIA_Pipeline for information on how to use the accompanying code.

Multiple instances of MARIA can sometimes generate slightly differing (<Δ1) MARIA percentile scores. If comparing the presentation prediction of two constructs, it is recommended to use the combine_split_MARIA_input.py function for local runs, or combine the inputs to be processed at once for computing cluster runs. In these cases, it is advised to use the "Gene Symbol" column to denote which construct each analysed peptide originated from for downstream processing.
