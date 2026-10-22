// Stage an existing Cell Ranger (count or arc) output into a Cell Ranger-count-like layout made of
// symlinks, so downstream steps never write into the original run directory.
process STAGE_CELLRANGER {
    tag "$meta.id"
    label 'process_med'

    input:
    tuple val(meta), path(cr_dir, stageAs: 'cr_in')

    output:
    tuple val(meta), path("${meta.id}")

    script:
    """
    src=\$(readlink -f cr_in)
    # accept either the sample directory (containing outs/) or the outs directory itself
    [ -d "\$src/outs" ] && src="\$src/outs"

    mkdir -p ${meta.id}/outs
    for f in raw_feature_bc_matrix.h5 filtered_feature_bc_matrix.h5 filtered_feature_bc_matrix; do
        [ -e "\$src/\$f" ] || { echo "Missing \$src/\$f" >&2; exit 1; }
        ln -s "\$src/\$f" ${meta.id}/outs/\$f
    done

    # count: possorted_genome_bam.bam ; arc: gex_possorted_bam.bam (both exposed under the count name)
    if   [ -e "\$src/possorted_genome_bam.bam" ]; then bam="\$src/possorted_genome_bam.bam"
    elif [ -e "\$src/gex_possorted_bam.bam" ];    then bam="\$src/gex_possorted_bam.bam"
    else echo "No possorted_genome_bam.bam or gex_possorted_bam.bam in \$src" >&2; exit 1; fi
    [ -e "\$bam.bai" ] || { echo "Missing BAM index \$bam.bai" >&2; exit 1; }
    ln -s "\$bam"     ${meta.id}/outs/possorted_genome_bam.bam
    ln -s "\$bam.bai" ${meta.id}/outs/possorted_genome_bam.bam.bai
    """
}
