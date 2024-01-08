with open(config['SAMPLES']) as fp:
    SAMPLES = fp.read().splitlines()

print(SAMPLES) 

rule all:
         input:
            expand("{sample}_QC_vlnplot.pdf", sample = SAMPLES),
            "merged.rds"  

rule analyse:
      input: 
           "{sample}_atac_fragments.tsv.gz", 
           "{sample}_filtered_feature_bc_matrix.h5" 
      params: 
          "{sample}"
      output: 
           "{sample}_QC_vlnplot.pdf",
 
      shell: 
           """
           Rscript scATAC.R {params}  
           """

rule merge: 
      input: 
          expand("{sample}.rds", sample =SAMPLES),
      params:
          "merged"  
      output: 
          "merged.rds" 
      shell: 
          """
          Rscript merge.R {input} {output} 
          """
 
