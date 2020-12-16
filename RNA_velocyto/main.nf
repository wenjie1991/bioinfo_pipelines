nextflow.enable.dsl=2

work_dir = "/data/users/wsun/workdir/rna_velo/"
cellranger_ref="/bioinfo/local/build/Centos/cellranger/refdata/refdata-cellranger-GRCh38-3.0.0"
custom_path= "$PATH:/bioinfo/local/build/Centos/cellranger/cellranger-3.1.0"
repeat_msk="/data/users/wsun/annotation/GRCh38_rmsk.gtf"
ref_gtf="$cellranger_ref/genes/genes.gtf"

params.sample_table
params.output_dir = "$work_dir/data"


process download_data {

    cpus 3
    time '3h'
    memory '3 GB'
    executor 'pbs'


    input: val "patient_name"

    output: tuple val("$patient_name"), path("fastq")

    script:
    """
    mkdir -p fastq
    gsutil -m cp -r gs://fc-d0d37f95-5eba-42ea-8967-a2c4a4b676c4/${patient_name} fastq
    """
}

process run_cell_ranger {

    cpus 10
    time '100h'
    memory '50 GB'
    executor 'pbs'


    input: tuple val(patient_name), path(fastq), val(indices)

    output: tuple val("$patient_name"), path("${patient_name}_${indices}"), val(indices)

    script:
    job_name="${patient_name}_${indices}"
    """
    export PATH=$custom_path
    rm -rf $patient_name
    cellranger count --id=${job_name} \
        --fastqs=${fastq}/${patient_name} \
        --indices=${indices} --localcores=10 \
        --localmem=48 --transcriptome=${cellranger_ref}
    """
}

process run_rnavelocyto {

    cpus 5
    time '50h'
    memory '50 GB'
    executor 'pbs'

    publishDir "${params.output_dir}/", mode: 'symlink'

    input: tuple val(patient_name), path(cell_ranger_output), val(indices)

    output: 
    path "${patient_name}_${indices}.log"
    path cell_ranger_output

    script:
    """
    export PS1=
    source ~/.bashrc
    conda activate velocyto
    eval `modulecmd bash load samtools`
    rm -rf $cell_ranger_output/velocyto
    velocyto run10x -@ 4 --samtools-memory 1000 -m $repeat_msk $cell_ranger_output $ref_gtf > ${patient_name}_${indices}.log
    """
}

workflow {
    indices = Channel.fromPath(params.sample_table).splitCsv()
    patient_name = indices.map { it[0] }.unique { it }
    patient_name | download_data | combine(indices, by: 0) | map { it } | run_cell_ranger | run_rnavelocyto | mix | view
}
