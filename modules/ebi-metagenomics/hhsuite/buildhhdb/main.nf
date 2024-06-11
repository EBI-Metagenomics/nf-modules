process HHSUITE_BUILDHHDB {
    tag "$meta.id"
    label 'process_high'

    container 'microbiome-informatics/hh-suite-db-builder:1.0.0'

    input:
    tuple val(meta), path(a3m)

    output:
    tuple val(meta), path("${meta.id}_hh_db"), emit: hh_db
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_hh_db"

    """
    ffindex_build -s ${a3m}_a3m.ff{data,index} ${a3m}
    mpirun -np ${task.cpus} ffindex_apply_mpi ${a3m}_a3m.ff{data,index} -i ${a3m}_hhm.ffindex -d ${a3m}_hhm.ffdata -- hhmake -i stdin -o stdout -v 0
    mpirun -np ${task.cpus} cstranslate_mpi -f -x 0.3 -c 4 -I a3m -i ${a3m}_a3m -o ${a3m}_cs219
    sort -k3 -n -r ${a3m}_cs219.ffindex | cut -f1 > sorting.dat
    ffindex_order sorting.dat ${a3m}_hhm.ff{data,index} ${a3m}_hhm_ordered.ff{data,index}
    mv ${a3m}_hhm_ordered.ffindex ${a3m}_hhm.ffindex
    mv ${a3m}_hhm_ordered.ffdata ${a3m}_hhm.ffdata
    ffindex_order sorting.dat ${a3m}_a3m.ff{data,index} ${a3m}_a3m_ordered.ff{data,index}
    mv ${a3m}_a3m_ordered.ffindex ${a3m}_a3m.ffindex
    mv ${a3m}_a3m_ordered.ffdata ${a3m}_a3m.ffdata

    mkdir ${prefix}
    mv ${a3m}_* ${prefix}/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hh-suite: \$(hhblits -h | grep 'HHblits' | sed -n -e 's/.*\\([0-9]\\+\\.[0-9]\\+\\.[0-9]\\+\\).*/\\1/p')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_hh_db"

    """
    mkdir ${prefix}
    touch ${prefix}/test_a3m_a3m.ffdata
    touch ${prefix}/test_a3m_a3m.ffindex
    touch ${prefix}/test_a3m_cs219.ffdata
    touch ${prefix}/test_a3m_cs219.ffindex
    touch ${prefix}/test_a3m_cs219.log1
    touch ${prefix}/test_a3m_hmm.ffdata
    touch ${prefix}/test_a3m_hmm.ffindex

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hh-suite: \$(hhblits -h | grep 'HHblits' | sed -n -e 's/.*\\([0-9]\\+\\.[0-9]\\+\\.[0-9]\\+\\).*/\\1/p')
    END_VERSIONS
    """
}
