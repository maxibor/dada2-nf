#!/usr/bin/env nextflow

version = "0.2"
version_date = "September 16th, 2019"

def helpMessage() {
    log.info"""
     dada2: simple dada2 16s classifier pipeline
     Homepage: https://github.com/maxibor/dada2-nf
     Author: Maxime Borry <borry@shh.mpg.de>
    =========================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run maxibor/dada2-nf --reads '/path/to/paired_end_reads_*.{1,2}.fastq.gz'
    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)

    Settings:
      --phred                       Specifies the fastq quality encoding (33 | 64). Defaults to ${params.phred}
      --pairedEnd                   Specifies if reads are paired-end (true | false). Default = ${params.pairedEnd}
      --silva_db                    Silva database for dada2. Default = ${params.silva_db}
      --silva_species_db            Silva species db for dada2. Default = ${params.silva_species_db}
      --rank                        Taxonomic rank to retain (Genus | Species). Default = ${params.rank}


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
summary['Silva species DB'] = params.silva_species_db
summary['Rank retained'] = params.rank
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="



Channel
    .fromFilePairs( params.reads, size: params.pairedEnd ? 2 : 1 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\n" }
	.set {reads_to_trim}

if (params.silva_db == '' || params.silva_species_db == ''){
    println('Downloading SILVA Databases\n')
    process download_db {

        output:
            file('silva_nr_v132_train_set.fa') into silva_db
            file('silva_species_assignment_v132.fa') into silva_species_db
        script:
            """
            wget https://zenodo.org/record/1172783/files/silva_nr_v132_train_set.fa.gz?download=1 -O silva_nr_v132_train_set.fa.gz
            gunzip silva_nr_v132_train_set.fa.gz
            wget https://zenodo.org/record/1172783/files/silva_species_assignment_v132.fa.gz?download=1 -O silva_species_assignment_v132.fa.gz
            gunzip silva_species_assignment_v132.fa.gz
            """
    }
} else {
    Channel
        .fromPath(params.silva_db)
        .set {silva_db}
    Channel
        .fromPath(params.silva_species_db)
        .set {silva_species_db}
}




process AdapterRemoval {
    tag "$name"

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
        file(silva) from silva_db.first()
        file(silva_species) from silva_species_db.first()
    output:
        set val(name), file("*.species_dada2.csv") into dada_out
        set val(name), file("*.dada2.csv") into dada_classify
        set val(name), file("*.read_count.csv") into dada_read_count_table
    script:
        outname = name+".dada2.csv"
        species_out = name+".species_dada2.csv"
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

            silva_db = "${silva}"
            silva_species_db = "${silva_species}"

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
            taxa <- assignTaxonomy(seqtab.nochim, silva_db, tryRC = TRUE, taxLevels = c("Genus","Species"), multithread=${task.cpus})
            write.csv(taxa, "$outname")

            # Assign species 
            taxa2 = addSpecies(taxtab = taxa, refFasta = silva_species_db)

            rownames(mergers) = mergers[,'sequence']
            
            # Write results to disk 
            print(taxa2)

            spec = taxa2[,c(ncol(taxa2)-1,ncol(taxa2))]
            colnames(spec) = c('Genus','Species')
            spec_abund = merge(mergers, spec, by = 'row.names')[,c('Genus','Species','abundance')]


            write.csv(spec_abund, "$species_out")
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
            silva_species_db = "${params.silva_species_db}"

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
            taxa <- assignTaxonomy(seqtab.nochim, silva_db, tryRC = TRUE, taxLevels = c("Genus","Species"), multithread=${task.cpus})
            write.csv(taxa, "$outname")

            # Assign species 
            taxa2 = addSpecies(taxtab = taxa, refFasta = silva_species_db)

            rownames(mergers) = mergers[,'sequence']
            
            # Write results to disk 
            print(taxa2)

            spec = taxa2[,c(ncol(taxa2)-1,ncol(taxa2))]
            colnames(spec) = c('Genus','Species')
            spec_abund = merge(mergers, spec, by = 'row.names')[,c('Genus','Species','abundance')]


            write.csv(spec_abund, "$species_out")
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
        dada2taxo.py -s $name -r ${params.rank} $dd
        """
}

process dada_merge {

    label 'ristretto'

    publishDir "${params.results}/merged", mode: 'copy'

    input:
        file(csv_count) from dada_taxo.collect()

    output:
        file('dada2_otu_table.csv') into dada_merged

    script:
        out = "dada2_otu_table.csv"
        """
        dada_merge.py -o $out
        """    
}