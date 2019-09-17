#!/usr/bin/env nextflow

version = "0.2"
version_date = "September 16th, 2019"

def helpMessage() {
    log.info"""
     megahit-nf: simple Megahit assembler Nextflow pipeline
     Homepage: https://github.com/maxibor/megahit-nf
     Author: Maxime Borry <borry@shh.mpg.de>
    =========================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run maxibor/megahit-nf --reads '/path/to/paired_end_reads_*.{1,2}.fastq.gz'
    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)

    Settings:
      --phred                       Specifies the fastq quality encoding (33 | 64). Defaults to ${params.phred}
      --pairedEnd                   Specifies if reads are paired-end (true | false). Default = ${params.pairedEnd}
      --silva_db                    Silva database for dada2. Default = ${params.silva_db}
      --silva_specie_db             Silva species db for dada2. Default = ${silva_specie_db}


    Options:
      --results                     The output directory where the results will be saved. Defaults to ${params.results}
      --help  --h                   Shows this help page
    """.stripIndent()
}

// Show help message
if (params.help){
    helpMessage()
    exit 0
}

//Logging parameters
log.info "========================================="
log.info " Pipeline dada_pipe ${version}"
log.info " Last updated on ${version_date}"
log.info "========================================="
def summary = [:]
summary['Reads'] = params.reads
summary['phred quality'] = params.phred
summary['Paired end'] = params.pairedEnd
summary['Silva DB'] = params.silva_db
summary['Silva species DB'] = params.silva_specie_db
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="



Channel
    .fromFilePairs( params.reads, size: params.pairedEnd ? 2 : 1 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\n" }
	.set {reads_to_trim}

process AdapterRemoval {
    tag "$name"

    label 'adaprem'

    input:
        set val(name), file(reads) from reads_to_trim

    output:
        set val(name), file('*.trimmed.fastq') into trimmed_reads
        file("*.settings") into adapter_removal_results

    script:
        out1 = name+".pair1.trimmed.fastq"
        out2 = name+".pair2.trimmed.fastq"
        se_out = name+".trimmed.fastq"
        non_col1 = name+".discarded.1.fq"
        non_col2 = name+".discarded.2.fq"
        settings = name+".settings"
        if (params.pairedEnd){
            """
            AdapterRemoval --basename $name --file1 ${reads[0]} --file2 ${reads[1]} --trimns --trimqualities --maxns 0 --minquality 20 --minlength 30 --output1 $out1 --output2 $out2 --threads ${task.cpus} --qualitybase ${params.phred} --settings $settings
            """
        } else {  
            """
            AdapterRemoval --basename $name --file1 ${reads[0]} --trimns --trimqualities --minquality 20 --maxns 0 --minlength 30 --output1 $se_out --threads ${task.cpus} --qualitybase ${params.phred} --settings $settings
            """
        }      
}

process dada2 {
    tag "$name"

    label 'dada'

    errorStrategy 'ignore'

    publishDir "${params.results}/dada", mode: 'copy'

    input:
        set val(name), file(fq) from trimmed_reads
    output:
        set val(name), file("*.dada2.csv") into dada_out
        set val(name), file("*.read_count.csv") into dada_read_count_table
    script:
        outname = name+".dada2.csv"
        read_count_name = name+".read_count.csv"
        if (params.pairedEnd) {
            """
            #!/usr/bin/env Rscript
            
            library(dada2)
            library(plyr)

            setwd(".")

            sample.name = "$name"

            fwd = "${fq[0]}"
            rev = "${fq[1]}"

            silva_db = "${params.silva_db}"
            silva_specie_db = "${params.silva_specie_db}"

            # Dereplication            
            derepFs <- derepFastq(fwd, verbose=TRUE)
            derepRs <- derepFastq(rev, verbose=TRUE)

            # Learning error rate
            err_forward_reads <- learnErrors(derepFs, multithread = ${task.cpus})
            err_reverse_reads <- learnErrors(derepRs, multithread=${task.cpus})

            # Sample Inference
            dadaFs <- dada(derepFs, err=err_forward_reads, multithread=${task.cpus})
            dadaRs <- dada(derepRs, err=err_reverse_reads, multithread=${task.cpus})

            # Merge paired reads 
            mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

            # Make sequence table
            seqtab <- makeSequenceTable(mergers)

            # Remove chimeras
            seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=${task.cpus}, verbose=TRUE)

            # Track read numbers
            getN <- function(x) sum(getUniques(x))
            track <- cbind(getN(dadaFs), getN(dadaRs), getN(mergers), rowSums(seqtab.nochim))
            colnames(track) <- c("denoisedF", "denoisedR", "merged", "nonchim")
            rownames(track) <- sample.name
            write.csv(track, "${read_count_name}")

            # Assign taxonomy 
            taxa <- assignTaxonomy(seqtab.nochim, silva_db, tryRC = TRUE, multithread=${task.cpus})

            # Assign species 
            taxa = addSpecies(taxtab = taxa, refFasta = silva_specie_db)
            
            # Write results to disk 
            print(taxa)

            spec_count = table(paste(as.data.frame(taxa)[,"Genus"], as.data.frame(taxa)[,"Species"]))

            print(spec_count)

            write.csv(spec_count, "$outname")
            """
        } else {
            """
            #!/usr/bin/env Rscript
            
            library(dada2)
            library(plyr)

            setwd(".")

            sample.name = "$name"

            fwd = "${fq[0]}"

            silva_db = "${params.silva_db}"
            silva_specie_db = "${params.silva_specie_db}"

            # Dereplication            
            derepFs <- derepFastq(fwd, verbose=TRUE)

            # Learning error rate
            err_forward_reads <- learnErrors(derepFs, multithread = ${task.cpus})

            # Sample Inference
            dadaFs <- dada(derepFs, err=err_forward_reads, multithread=${task.cpus})

            # Make sequence table
            seqtab <- makeSequenceTable(dadaFs)

            # Remove chimeras
            seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=${task.cpus}, verbose=TRUE)

             # Track read numbers
            getN <- function(x) sum(getUniques(x))
            track <- cbind(getN(dadaFs), rowSums(seqtab.nochim))
            colnames(track) <- c("denoisedF", "nonchim")
            rownames(track) <- sample.name
            write.csv(track, "${read_count_name}")

            # Assign taxonomy 
            taxa <- assignTaxonomy(seqtab.nochim, silva_db, tryRC = TRUE, multithread=${task.cpus})

            # Assign species 
            taxa = addSpecies(taxtab = taxa, refFasta = silva_specie_db)

            spec_count = table(paste(as.data.frame(taxa)[,"Genus"], as.data.frame(taxa)[,"Species"]))
            
            # Write results to disk 
            write.csv(spec_count, "$outname")
            """
        }
        
}



process dada2_to_taxo {
    tag "$name"

    label 'ristretto'

    errorStrategy 'ignore'

    publishDir "${params.results}/taxo", mode: 'copy'

    input:
        set val(name), file(dd) from dada_out
    output:
        set val(name), file("*.dadataxo.csv") into dada_taxo
    script:
        """
        dada2taxo.py -s $name $dd
        """
}