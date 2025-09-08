process DATA_PREPROCESSING {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::pandas=1.5.2 bioconda::numpy=1.21.6 bioconda::scipy=1.9.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'quay.io/biocontainers/pandas:1.5.2' }"

    input:
    path expression_matrix
    path tissue_attributes
    path subject_attributes
    path annotation_gtf

    output:
    path "processed_data.h5"        , emit: processed_data
    path "preprocessing_stats.json" , emit: multiqc_files
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    data_preprocessing.py \\
        --expression-matrix ${expression_matrix} \\
        --tissue-attributes ${tissue_attributes} \\
        --subject-attributes ${subject_attributes} \\
        --annotation-gtf ${annotation_gtf} \\
        --output processed_data.h5 \\
        --stats-output preprocessing_stats.json \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch processed_data.h5
    touch preprocessing_summary.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
    END_VERSIONS
    """
}