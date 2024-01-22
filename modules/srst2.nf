// run srst2

process SRST2 {
    tag "${prefix}"
    errorStrategy { task.exitStatus in 148 ? 'ignore' : 'terminate' }

    input:
    tuple val(prefix), path(reads)
    val(procdir)
    path(srst2DB)

    output:
    publishDir "$procdir/${prefix}/srst2/${srst2DB.name}", mode: 'copy'
    path("${prefix}*"), emit: output
    val("$prefix"), emit: prefix
    stdout emit: stdout

    script:
    def outputdir = file("$procdir/${prefix}/srst2/${srst2DB.name}")
    if (outputdir.exists() && !params.overwrite) {
        println "$outputdir exists. Skipping $prefix ..."
        """
        exit 148
        """
    }else{
        """
        srst2 --input_pe $reads --output $prefix --log --gene_db $srst2DB --threads ${task.cpus}
        # Get reads
        if [[ -f "${prefix}__fullgenes__CARD_v3.0.8_SRST2__results.txt" ]]; then
            > temp.txt
            for i in \$(cut -f13 ${prefix}__fullgenes__CARD_v3.0.8_SRST2__results.txt); do
              counts=\$(samtools view *.bam | cut -f3 | grep \$i | wc -l)
              echo \$counts >> temp.txt
            done
            sed -i "1 s/0/reads/" temp.txt
            paste ${prefix}__fullgenes__CARD_v3.0.8_SRST2__results.txt temp.txt >> ${prefix}_results.csv
            rm temp.txt
            # Get total reads
            totalreads=\$(zcat $reads | wc -l)
            totalreads=\$((totalreads/2))
            echo $totalreads >> ${prefix}_totalreads.csv
            ######## RPKM, FPKM, TPM https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/
            # Get RPKM
            awk -v OFS='\\t' '{print \$0, '\$totalreads'}' ${prefix}_results.csv > ${prefix}_results_temp.csv
            awk -F '\\t' -v OFS='\\t' '{\$(NF+1)=(\$15/\$16*1000000/\$10)} {print \$0}' ${prefix}_results_temp.csv > ${prefix}_results_temp2.csv
            awk -F '\\t' 'NR==1{\$17="RPKM"} NR==1{\$16="totalcounts"} {print \$0}' ${prefix}_results_temp2.csv > ${prefix}_results_RPKM.csv
            # Get TPM

            # Clear files
            rm ${prefix}_results_temp.csv ${prefix}_results_temp2.csv
        fi
        """
    }
}

//
// process SRST2_RPKM {
//     tag "${prefix}"
//     errorStrategy { task.exitStatus in 148 ? 'ignore' : 'terminate' }
//
//     input:
//     tuple val(prefix), val(results)
//     tuple val(prefix2), val(totalreads)
//     val(procdir)
//
//     output:
//     publishDir "$procdir/${prefix}/srst2/${srst2DB.name}", mode: 'copy'
//     path("${prefix}*"), emit: output
//     val("${prefix}"), emit: prefix
//     stdout emit: stdout
//
//     script:
//     """
//
//     """
// }
