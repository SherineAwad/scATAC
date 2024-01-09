with open(config['SAMPLES']) as fp:
    SAMPLES = fp.read().splitlines()

print(SAMPLES) 

rule all:
         input:
            expand("{sample}.rds", sample = SAMPLES),
            expand("{sample}_filtered.rds", sample = SAMPLES), 
            "merged.rds",  

rule analyse:
      input: 
           "{sample}_atac_fragments.tsv.gz", 
           "{sample}_filtered_feature_bc_matrix.h5" 
      params: 
          "{sample}"
      output: 
           "{sample}.rds",
 
      shell: 
           """
           Rscript preprocess.R {params}  
           """

rule filter: 
      input: 
          "{sample}.rds" 
      params: 
          "{sample}"
      output: 
          "{sample}_filtered.rds",  
      shell: 
          "Rscript filter.R {params}"
 
rule run_all: 
      input: 
          expand("{sample}_filtered.rds", sample =SAMPLES),
      params:
          "{sample}"  
      output: 
          "merged.rds" 
      shell: 
          """
          Rscript merge.R {input} {output} 
          Rscript geneExp.R {params} 
          Rscript callPeaks.R {params}
          """
 
