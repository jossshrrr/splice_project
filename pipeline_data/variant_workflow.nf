#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//input parameters
params.id = 'NA'
params.input_vcf = 'NA'
params.reference_fasta = 'NA'

params.filter = false
params.gene_bed = 'NA' 

params.gnomad_toml = 'NA'
params.gnomad_count = 100

params.clinvar_toml = 'NA'

params.spliceai_distance = 50
params.pangolin_distance = 50

// check input VCF file
if (!params.input_vcf) {
    exit 1, "No input VCF file provided by user"}

if (!file(params.input_vcf).exists()) {
    exit 1, "Invalid input VCF file provided by user"}

log.info """\
    ==============================================
        Splice Variant Annotation Workflow    
    ==============================================
    cohort id      : ${params.id}
    input vcf      : ${params.input_vcf}
    gene list      : ${params.gene_bed}
    ==============================================
    """
    .stripIndent()

//make output dir and set prefix
outputdir = file(params.id)
outputdir.mkdirs()

//create a subdirectory within outputdir
subdir = outputdir + "/processed"
subdir.mkdirs()

prefix = (params.id)
bed = (params.gene_bed)
fasta = (params.reference_fasta)

gnomad_toml = (params.gnomad_toml)
gnomad_count = (params.gnomad_count)

clinvar_toml = (params.clinvar_toml)
vep_cache = (params.vep_cache)

spliceai_distance = (params.spliceai_distance)
pangolin_distance = (params.pangolin_distance)


workflow {

    input_vcf = Channel.fromPath(params.input_vcf).map { tuple("$it", "${it}.tbi")}

    outputs = Gene_Filter(input_vcf,bed)
    | Gnomad
    | Gnomad_Filter
    | Split_Vcf
    | flatten
    | SpliceAI
    | Pangolin
    | SQUIRLS
    | CADD
    | VEP
    | ClinVar

    annotated_vcf = Join_VCF(outputs.collect())

    Process_Splice_Scores(annotated_vcf)
    | Filter_Variants
    
}


process Gene_Filter {
    label 'C1M1T1'

    module 'bcftools'

    input:
        tuple path(vcf), path(vcf_index)
        path(bed) 

    output:
        path("*filtered.vcf")

    shell:
    '''
    if [ "!{params.filter}" = true ] ; then
        bcftools filter --regions-file !{bed} !{vcf} > gene_regions.vcf
        bcftools norm --fasta-ref !{fasta} --multiallelics '-' gene_regions.vcf | bcftools view -Ov -o filtered.vcf
    else
        bcftools norm --fasta-ref !{fasta} --multiallelics '-' !{vcf} | bcftools view -Ov -o filtered.vcf
    fi
    '''
}

process Gnomad {
    label 'C1M1T1'

    conda 'bioconda::vcfanno=0.3.3'

    input:
        path(vcf) 
        
    output:
        path("*.gnomad.vcf")
    
    shell:
    '''
    output=$(basename $PWD)
    vcfanno ${workflow.projectDir}/pipeline_data/vcfanno_files/gnomad_toml !{vcf} > $output.gnomad.vcf 
    '''
}

process Gnomad_Filter {
    label 'C1M1T1'

    module 'gatk'   

    input:
        path(vcf) 

    output:
        path("*.filter.vcf*")
    
    shell:
    '''
    output=$(basename $PWD)
    gatk --java-options "-Xmx7g" VariantFiltration \
    -R !{fasta} \
    -V !{vcf}  \
    -O $output.flag.vcf \
    --filter-name "High_Pop_Freq" \
    --filter-expression "AC_gnomad3.0 > !{gnomad_count}"
    
    grep -v High_Pop_Freq $output.flag.vcf > $output.filter.vcf          
    '''
}

process Split_Vcf {
    label 'C1M1T1'
	
    conda 'bioconda::snpsift'	

    input:
        path(vcf) 

    output:
        path("*.vcf*")

    shell:
    '''
    variant_number=`expr $(zgrep -v "^#" !{vcf} | wc -l)`
    
    if [ $variant_number -gt 1000 ]
    then
        task_number=500
    elif [ $variant_number -gt 10 ]
    then
        task_number=10
    else
        task_number=2
    fi
    
    split_number=`expr $variant_number / $task_number`

    split_number=`expr $variant_number / $task_number`
    if [ $split_number -lt 5 ]; then
        split_number=5
    fi

    echo "Variant number = $variant_number" >> error_check.txt
    echo "Number of tasks = $task_number" >> error_check.txt
    echo "Number of variants per split = $split_number" >> error_check.txt

    SnpSift split -l $split_number !{vcf} 
    '''
}

process SpliceAI {
    label 'C1M4T2'

    container 'joshrrr/spliceai'	

    input:
        path(vcf) 

    output:
        path("spliceai.vcf")

    script:
    """
    spliceai -I $vcf \
    -O spliceai.vcf \
    -R $fasta \
    -D $spliceai_distance \
    -A ${workflow.projectDir}/pipeline_data/spliceai/gencode.v38.annotation.txt  
    """
}

process Pangolin {
    label 'C1M4T1'
        
    container 'joshrrr/pangolin'

    input:
        path(vcf)

    output:
        path("pangolin_file*")

    script:
    """
    pangolin $vcf \
    $fasta \
    ${workflow.projectDir}/pipeline_data/pangolin/gencode.v38.annotation.db \
    -d $pangolin_distance  \
    pangolin_file
    """
}

process SQUIRLS {
    label 'C1M4T1'

    input:
        path(vcf) 

    output:
        path("squirls.vcf")

    script:
    """
	java -jar ${workflow.projectDir}/pipeline_data/squirls \
	annotate-vcf --all-transcripts=true -f=vcf -d ${workflow.projectDir}/pipeline_data/squirls/squirls_data \
	$vcf ./squirls 
    """
}

process CADD {
    label 'C1M4T1'

    module 'htslib'

    conda 'bioconda::vcfanno=0.3.3'

    input:
        path(vcf) 
        
    output:
        path("*.cadd_run.vcf")
    
    shell:
    '''
	output=$(basename $PWD)

	#CADD processes variants in a stripped VCF format (i.e. no other annotations or genotype info)
	#Preprocessing steps below convert the input VCF to this format (tmp.vcf)
	grep -v "^#" !{vcf} | cut -f 1-5 | sed 's/^chr//' > tmp.vcf

	#Apply CADD scores to variants in tmp.vcf (scores are output to a file tmp.tsv.gz)
	${workflow.projectDir}/pipeline_data/cadd/CADD.sh tmp.vcf
	tabix -f -b 2 -e 2 -s 1 tmp.tsv.gz

	#Create .toml file to annotate the original VCF file with CADD scores
	cadd_path=$(realpath tmp.tsv.gz)
	echo '[[annotation]]' > cadd.toml
	echo "file= '${cadd_path}'" >> cadd.toml
	echo 'names=["CADD_Score"]' >> cadd.toml
	echo 'ops=["mean"]' >> cadd.toml
	echo 'columns=[6]' >> cadd.toml

	vcfanno cadd.toml !{vcf} > $output.cadd_run.vcf 
	'''
}

process VEP {
    label 'C1M1T1'

    module 'ensembl-vep'
    
    input:
        path(vcf)

    output:
        path("*.vep.vcf")

    shell:
    '''
    output=$(basename $PWD)
    vep --cache --dir ${workflow.projectDir}/pipeline_data/vep_cache --cache_version 104 --assembly GRCh38 \
    -i !{vcf} -o $output.vep.vcf --format vcf --vcf --symbol --terms SO --tsl --hgvs \
    --fasta !{fasta} --offline --sift b --polyphen b --ccds --hgvs --hgvsg --symbol \
    --numbers --protein --af --af_1kg --af_gnomad --max_af --variant_class --pick_allele_gene
    '''
}

process ClinVar {
    label 'C1M1T1'

    conda 'bioconda::vcfanno=0.3.3'

    input:
        path(vcf)
                
    output:
        path("*.clinvar.vcf")

    shell:
    '''
    output=$(basename $PWD)
    vcfanno ${workflow.projectDir}/pipeline_data/vcfanno_files/clinvar_toml !{vcf} > $output.clinvar.vcf
    '''
}

process Join_VCF {
    label 'C1M1T1'

    module 'bcftools:htslib'

    conda 'bioconda::snpsift'

    publishDir "$outputdir", mode: "copy"

    input:
        path("*.clinvar.vcf") 

    output:
        path("*.gz")
        path("*.tbi")

    shell:
    '''
    SnpSift split -j *.clinvar.vcf | bcftools sort -Oz -o final.vcf.gz
    tabix final.vcf.gz
    '''
}

process Process_Splice_Scores {
    label 'C1M1T1'

	conda 'bioconda::vcfanno=0.3.3'	

	input:
		path(vcf)
        path(vcf_index)
				
	output:
		path("*processed.vcf")

	shell:
	'''
	vcfanno -p 1 -lua ${workflow.projectDir}/pipeline_data/vcfanno_files/spliceai.lua ${workflow.projectDir}/pipeline_data/vcfanno_files/spliceai_postprocess.toml !{vcf} > file1.vcf
	vcfanno -p 1 -lua ${workflow.projectDir}/pipeline_data/vcfanno_files/pangolin.lua ${workflow.projectDir}/pipeline_data/vcfanno_files/pangolin_postprocess.toml file1.vcf > file2.vcf
	vcfanno -p 1 -lua ${workflow.projectDir}/pipeline_data/vcfanno_files/squirls.lua ${workflow.projectDir}/pipeline_data/vcfanno_files/squirls_postprocess.toml file2.vcf > !{prefix}_processed.vcf
	'''
}

process Filter_Variants {
    label 'C1M1T1'

    module 'R/4.2.0'

    input:
        path(vcf) 

    output:
        path("*.csv")

    shell:
    '''
    Rscript ${workflow.projectDir}/pipeline_data/filter_variants.R !{vcf}  ${workflow.projectDir}/pipeline_data/EpilepsyGenes_v2022-09_Full.tsv
    '''
}

