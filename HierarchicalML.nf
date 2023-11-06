#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
    Given a directory of paired-end FASTQs, make a tab-delimited file listing
    each sample name and its R1 + R2 files separated by tabs on each line.
    All of the files created per sample are collected in the workflow via collectFile into one file listing all of the samples
*/
process MAKE_TAB_FILE {
    debug true

    input:
    tuple val(sample_id), path(reads)
    val readdir

    output:
    path 'sample_reads.txt'

    script:
    """
    echo "${sample_id}\t${readdir}${reads[0]}\t${readdir}${reads[1]}" > sample_reads.txt
    """
}

/*
    Performs read QC, trimming, downsampling and kmer counting/unitig generation
    for all samples
*/
process GENERATE_UNITIGS {
    debug true
    conda params.condaenv1
    publishDir "${params.outdir}", mode:'copy'

    input:
    val scriptdir
    path tabfile
    val cpus

    output:
    path("unitigs")

    script:
    """
    ${scriptdir}/unitig_pipeline.pl --reads ${tabfile} -o unitigs/ --cpus ${cpus}
    """
}

/*
    Call unitigs against query unitigs in the UKHSA dataset
*/
process CALL_UNITIGS_FROM_LIST {
    debug true
    conda params.condaenv1
    publishDir "${params.outdir}/unitigs", mode:'copy'

    input:
    val scriptdir
    path reads
    path tabfile
    val datadir
    val cpus

    output:
    path("processed_to_patterns")

    script:
    """
    ${scriptdir}/call_unitigs_from_list.pl -r ${reads}/ -o processed_to_patterns/ -l ${tabfile} -q ${datadir}/unitigs_unitigs.fasta -t ${cpus}
    """
}

/*
    Convert the variants to known patterns
*/
process CONVERT_TO_PRESENT_PATTERNS {
    debug true
    conda params.condaenv1
    publishDir "${params.outdir}/unitigs/processed_to_patterns", pattern: "*.tab", mode:'copy'

    input:
    val scriptdir
    path patterns_dir
    val datadir

    output:
    path("${patterns_dir}")

    script:
    """
    ${scriptdir}/variants_to_known_patterns -i ${patterns_dir}/counts.rtab -o ${patterns_dir}/patterns.tab -c ${datadir}/pattern_conversion.tab ${datadir}/patterns_in_model.txt
    """
}

/*
    Create sample_list and sample.location files necessary to run the
    classify_new_samples script
*/
process CREATE_CLASSIFICATION_FILES {
    debug true
    shell '/bin/bash'
    publishDir "${params.outdir}/unitigs/processed_to_patterns", mode:'copy'

    input:
    path tabfile
    path patterns_dir

    output:
    path("sample_list.txt"), emit: sample_list
    path("sample.location"), emit: sample_location

    script:
    """
    cut -f1 -d '\t' ${tabfile} > sample_list.txt
    cut -f1 -d '\t' ${tabfile} | awk 'BEGIN { FS = OFS = "\t" } { \$(NF+1)="Poland"; print \$0 }' > sample.location
    """
}

/*
    Copy in necessary files to working directory and classify the samples
*/
process CLASSIFY_SAMPLES {
    debug true
    conda params.condaenv2
    publishDir "${params.outdir}/unitigs/", pattern: '*/*.tsv', mode:'copy'
    publishDir "${params.outdir}/unitigs/", pattern: '*/*.tab', mode:'copy'
    publishDir "${params.outdir}/unitigs/", pattern: '*/*.rtab', mode:'copy'
    publishDir "${params.outdir}/unitigs/", pattern: '*/*.fasta', mode:'copy'
    publishDir "${params.outdir}/unitigs/", pattern: '*/*.gfa', mode:'copy'
    publishDir "${params.outdir}/unitigs/", pattern: '*/*.txt', mode:'copy'
    publishDir "${params.outdir}/unitigs/", pattern: '*/*.bfg_colors', mode:'copy'

    input:
    path patterns_dir
    path sample_list
    path sample_location
    path notebook

    output:
    path("${patterns_dir}/*")

    script:
    """
    cp ${notebook} ${patterns_dir}
    cp ${sample_list} ${patterns_dir}
    cp ${sample_location} ${patterns_dir}
    cp ${params.optimised_model} ${patterns_dir}
    cp ${params.hc_package} ${patterns_dir}
    cd ${patterns_dir}
    ipython classify_new_samples.ipynb
    """
}


workflow {
    reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    tabfile = MAKE_TAB_FILE(reads_ch, params.outdir)
    .collectFile(name: params.tabfile_name, storeDir: params.outdir)

    unitig_path = GENERATE_UNITIGS(params.scriptdir, tabfile, params.cpus)

    patterns_dir = CALL_UNITIGS_FROM_LIST(params.scriptdir, unitig_path, tabfile, params.datadir, params.cpus)

    patterns_dir2 = CONVERT_TO_PRESENT_PATTERNS(params.scriptdir, patterns_dir, params.datadir)

    CREATE_CLASSIFICATION_FILES(tabfile, patterns_dir)

    CLASSIFY_SAMPLES(patterns_dir2, CREATE_CLASSIFICATION_FILES.out.sample_list, CREATE_CLASSIFICATION_FILES.out.sample_location, params.notebook)
}
