with open(config['SAMPLES']) as fp:
    SAMPLES = fp.read().splitlines()

print(SAMPLES) 


rule all:
         input:
            expand("{sample}_preprocessed.rds", sample = SAMPLES),
            expand("{sample}_filtered.rds", sample = SAMPLES), 
            expand("{sample}_preprocessed.rds", sample = SAMPLES), 
            expand("{sample}_filtered.rds", sample = SAMPLES),
            expand("{myobject}.rds", myobject = config['OBJECT']),
            expand("{myobject}_analysed.rds", myobject = config['OBJECT']),
            expand("{myobject}_renamedClusters.rds", myobject = config['OBJECT']),
            expand("{myobject}_FinalClusters.rds", myobject = config['OBJECT']),
            expand("{myobject}_diffPeaks.rds", myobject = config['OBJECT']),
            expand("{myobject}_linkPeaks.rds", myobject = config['OBJECT']),
            expand("{myobject}_linkCoverage1.pdf",myobject = config['OBJECT']),
            expand("{myobject}_linkCoverage2.pdf", myobject = config['OBJECT']),
            expand("{myobject}_linkCoverage3.pdf", myobject = config['OBJECT']),
            expand("{myobject}_linkCoverage4.pdf", myobject = config['OBJECT']),
            #expand("{myobject}_peaks.rds", myobject = config['OBJECT']),

rule preprocess:
      input: 
           "{sample}_atac_fragments.tsv.gz", 
           "{sample}_filtered_feature_bc_matrix.h5" 
      params: 
          "{sample}"
      output: 
           "{sample}_preprocessed.rds",
 
      shell: 
           """
           Rscript preprocess.R {params}  
           """

rule filter: 
      input: 
          "{sample}_preprocessed.rds" 
      params: 
          "{sample}", 
          nATAC1 = config['nATAC1'],
	  nRNA1 = config['nRNA1'], 
          nATAC2= config['nATAC2'], 
          nRNA2= config['nRNA2'],
          features =config['features'], 
          nucl= config['nucl'], 
          tss= config['tss'],
          mt= config['mt']
      output: 
         "{sample}_filtered.rds",  
      shell: 
          "Rscript filter.R {params}"

rule merge: 
      input: 
         expand("{sample}_filtered.rds", sample =SAMPLES)
      params:
         lambda w: "merged ".join(expand("{sample}", sample =SAMPLES)), 
      output: 
        expand("{myobject}.rds", myobject=config['OBJECT']) 
      shell: 
        "Rscript merge.R {params}"
 
rule analyse: 
      input: 
          "{myobject}.rds"
      params: 
          "{myobject}"
      output: 
          "{myobject}_analysed.rds"
      shell: 
          "Rscript analyse.R {params}"

rule annotateClusters: 
      input: 
           "{myobject}_analysed.rds"  
      params: 
          "{myobject}"
      output: 
           "{myobject}_renamedClusters.rds" 
      shell: 
          "Rscript annotateClusters.R {params}" 
  
rule removeClusters:
      input:
           "{myobject}_renamedClusters.rds"
      params:
          "{myobject}"
      output:
           "{myobject}_FinalClusters.rds"
      shell:
          "Rscript removeClusters.R {params}"


rule diffPeaks:
      input:
           "{myobject}_FinalClusters.rds"
      params:
          "{myobject}"
      output:
           "{myobject}_diffPeaks.rds"
      shell:
          "Rscript diffPeaks.R {params}"


rule linkPeaks:
      input:
           "{myobject}_diffPeaks.rds"
      params:
          "{myobject}"
      output:
           "{myobject}_linkPeaks.rds"
      shell:
          "Rscript linkPeaks.R {params}"

rule plotPeaks:
      input:
           "{myobject}_linkPeaks.rds"
      params:
          "{myobject}"
      output:
           "{myobject}_linkCoverage1.pdf",
           "{myobject}_linkCoverage2.pdf",
           "{myobject}_linkCoverage3.pdf",
           "{myobject}_linkCoverage4.pdf",
      shell:
          "Rscript plotPeaks.R {params}" 
