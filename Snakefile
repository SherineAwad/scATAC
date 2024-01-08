with open(config['SAMPLES']) as fp:
    SAMPLES = fp.read().splitlines()

print(SAMPLES) 

rule all:
         input:
            expand("{sample}.rds", sample = SAMPLES),
            expand("{sample}_filtered.rds", sample = SAMPLES), 
            "merged.rds"  

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
           Rscript scATAC.R {params}  
           """

rule filter: 
      input: 
          "{sample}.rds" 
      output: 
          "{sample}_filtered.rds",  
      shell: 
          "Rscript filter.R {input}"
 
rule merge: 
      input: 
          expand("{sample}_filtered.rds", sample =SAMPLES),
      params:
          "merged"  
      output: 
          "merged.rds" 
      shell: 
          """
          Rscript merge.R {input} {output} 
          """
 
