version 1.0

import "structs.wdl"

task SplitBarcodes {
    input {
        File sample_sheet
        FastQSet fastq_set
    }

    Int disk_size = 10 + ceil(3*size(fastq_set.fastq_r1, "GiB")) + ceil(3*size(fastq_set.fastq_r1, "GiB"))
    String output_prefix = basename(fastq_set.fastq_r1, ".unmapped.1.fastq.gz")

    command <<<
        set -e

        mkdir outputs
        java -jar /usr/local/share/fgbio/fgbio.jar DemuxFastqs --inputs ~{fastq_set.fastq_r1} ~{fastq_set.fastq_r2} --read-structures 9B41T 50T -x ~{sample_sheet} -o outputs
        mv outputs/unmatched.bam unmatched.bam
    >>>

    runtime {
            docker: "quay.io/biocontainers/fgbio@sha256:5d61a29190c50cfc1d56e00c7390b4fadc59f975d5084e73c1324c4052d52fef"
            disks: "local-disk "+disk_size+" HDD"
            memory: "16 GB"
    }

    output {
        Array[File] sample_unmapped_bams = glob("outputs/*.bam")
        File no_sample_unmapped_bams = "unmatched.bam"
    }
}

task BuildCombinedName {
    input {
        Array[String]+ names
    }

    command <<<
        set -e

        echo ~{sep="_" names}.fna
    >>>

    runtime {
       docker: "ubuntu:19.10"
       disks: "local-disk 5 HDD"
       memory: "1 GB"
    }

    output {
        String combined_name = read_string(stdout())
    }
}

task BuildReferenceResources {
    input {
        Array[File]+ fastas
        Array[File]+ gffs
        String combined_name
    }

    Int disk_size = 10 + ceil(2.2*size(gffs, "GiB")) + ceil(2.2*size(fastas, "GiB"))

    command <<<
        set -e

        cat ~{sep=" " fastas} > ~{combined_name}.fna
        samtools faidx ~{combined_name}.fna
        samtools dict ~{combined_name}.fna > ~{combined_name}.dict
        /usr/gitc/bwa index ~{combined_name}.fna

        cat ~{sep=" " gffs} >  ~{combined_name}.gff
    >>>

    runtime {
       docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud@sha256:82dd1af86c9e6d4432170133382053525864d8f156a352e18ecf5947542e0b29"
       disks: "local-disk "+disk_size+" HDD"
       memory: "16 GB"
    }

    output {
        File combined_fasta = "~{combined_name}.fna"
        File combined_index = "~{combined_name}.fna.fai"
        File combined_dict = "~{combined_name}.dict"
        File combined_pac = "~{combined_name}.fna.pac"
        File combined_bwt = "~{combined_name}.fna.bwt"
        File combined_ann = "~{combined_name}.fna.ann"
        File combined_amb = "~{combined_name}.fna.amb"
        File combined_sa = "~{combined_name}.fna.sa"
        File combined_gff = "~{combined_name}.gff"
    }

}

task CreateAccensionFiles {
    input {
        Array[File] gffs
        Array[File] fnas
        Array[String] accs
        Int add3
        Int add5
    }

    Int disk_size = 10 + ceil(3.2*size(gffs, "GiB")) + ceil(3.2*size(fnas, "GiB"))

    command <<<
        set -e

        mkdir acc_links
        ln -s ~{sep=" " gffs} ~{sep=" " fnas} acc_links/
        mkdir out_dir
        mkdir temp_dir
        /root/tools/gff_parse3.sh acc_links out_dir temp_dir ~{add5} ~{add3} ~{sep=" " accs}
    >>>

    runtime {
        docker: "quay.io/ckachuli/microbial_rna_tag_seq:0.1"
        disks: "local-disk "+disk_size+" HDD"
        memory: "16 GB"
    }

    output {
        File all_gff = glob("out_dir/*_ALL.gff")[0]
        File genes_gff = glob("out_dir/*_GENES.gff")[0]
        File out_fna = glob("out_dir/*.fna")[0]
    }
}

task BuildBWAIndex {
    input {
        File fasta
    }

    Int disk_size = 10 + ceil(3*size(fasta, "GiB"))
    String algorithm = if ceil(size(fasta, "GiB")) >= 2 then "bwtsw" else "is"

    command <<<
        set -e

        mkdir output
        cd output
        /root/tools/bwa index -a ~{algorithm} ~{fasta}
    >>>

    runtime {
        docker: "quay.io/ckachuli/microbial_rna_tag_seq:0.1"
        disks: "local-disk "+disk_size+" HDD"
        memory: "16 GB"
    }

    output {
        Array[File] index_resources = glob("output/*")
    }
}

# Get version of BWA
task GetBwaVersion {
  command {
    # not setting set -o pipefail here because /bwa has a rc=1 and we dont want to allow rc=1 to succeed because
    # the sed may also fail with that error and that is something we actually want to fail on.
    /usr/gitc/bwa 2>&1 | \
    grep -e '^Version' | \
    sed 's/Version: //'
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud@sha256:82dd1af86c9e6d4432170133382053525864d8f156a352e18ecf5947542e0b29"
    memory: "1 GB"
  }
  output {
    String bwa_version = read_string(stdout())
  }
}

task AlignPairedReadsBWABacktrack {
    input {
        File unmapped_bam
        ReferenceFasta ref
        String bwa_version
    }

    Float ref_size = size(ref.ref_fasta, "GB") + size(ref.ref_fasta_index, "GB") + size(ref.ref_dict, "GB")
    Float bwa_size = size(ref.ref_amb, "GB") + size(ref.ref_ann, "GB") + size(ref.ref_bwt, "GB") + size(ref.ref_pac, "GB") + size(ref.ref_sa, "GB")
    Int disk_size = 10 + ceil(ref_size + bwa_size + 4.5*size(unmapped_bam, "GiB"))
    String name = basename(unmapped_bam, ".bam")

    command <<<
        set -e

        java -jar /usr/gitc/picard.jar \
              SamToFastq \
              INPUT=~{unmapped_bam} \
              F=f1.fastq \
              F2=f2.fastq \
              NON_PF=true

        /usr/gitc/bwa aln ~{ref.ref_fasta} f1.fastq > r1.sai
        /usr/gitc/bwa aln ~{ref.ref_fasta} f2.fastq > r2.sai
        /usr/gitc/bwa sampe ~{ref.ref_fasta} r1.sai r2.sai f1.fastq f2.fastq 2> >(tee bwa.stderr.log >&2) | \
        java -jar /usr/gitc/picard.jar \
            MergeBamAlignment \
            VALIDATION_STRINGENCY=SILENT \
            EXPECTED_ORIENTATIONS=FR \
            ALIGNED_BAM=/dev/stdin \
            UNMAPPED_BAM=~{unmapped_bam} \
            OUTPUT=~{name}.bam \
            REFERENCE_SEQUENCE=~{ref.ref_fasta} \
            CLIP_ADAPTERS=false \
            MAX_INSERTIONS_OR_DELETIONS=-1 \
            PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
            PROGRAM_RECORD_ID="bwasampe" \
            PROGRAM_GROUP_VERSION="~{bwa_version}" \
            PROGRAM_GROUP_COMMAND_LINE="bwa sampe" \
            PROGRAM_GROUP_NAME="bwasampe" \
            UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
            ALIGNER_PROPER_PAIR_FLAGS=true \
            UNMAP_CONTAMINANT_READS=true \
            ADD_PG_TAG_TO_READS=false \
            CREATE_INDEX=true
    >>>

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud@sha256:82dd1af86c9e6d4432170133382053525864d8f156a352e18ecf5947542e0b29"
        disks: "local-disk "+disk_size+" HDD"
        memory: "16 GB"
    }

    output {
        File aligned_bam = "~{name}.bam"
        File aligned_bam_index = "~{name}.bai"
        File bwa_sterr_log = "bwa.stderr.log"
    }
}

task AlignBWAMem {
    input {
        File unmapped_bam
        ReferenceFasta ref
        Array[File] fasta_resources
        String bwa_version
    }

    Float ref_size = size(ref.ref_fasta, "GB") + size(ref.ref_fasta_index, "GB") + size(ref.ref_dict, "GB")
    Float bwa_size = size(ref.ref_amb, "GB") + size(ref.ref_ann, "GB") + size(ref.ref_bwt, "GB") + size(ref.ref_pac, "GB") + size(ref.ref_sa, "GB")
    Int disk_size = 10 + ceil(ref_size + bwa_size + 2.5*size(unmapped_bam, "GiB"))
    String name = basename(unmapped_bam, ".bam")
    command <<<
        set -e

        java -jar /usr/gitc/picard.jar \
              SamToFastq \
              INPUT=~{unmapped_bam} \
              FASTQ=/dev/stdout \
              INTERLEAVE=true \
              NON_PF=true | \
        /usr/gitc/bwa mem /dev/stdin - 2> >(tee bwa.stderr.log >&2) | \
             java -jar /usr/gitc/picard.jar \
               MergeBamAlignment \
               VALIDATION_STRINGENCY=SILENT \
               EXPECTED_ORIENTATIONS=FR \
               ATTRIBUTES_TO_RETAIN=X0 \
               ATTRIBUTES_TO_REMOVE=NM \
               ATTRIBUTES_TO_REMOVE=MD \
               ALIGNED_BAM=/dev/stdin \
               UNMAPPED_BAM=~{unmapped_bam} \
               OUTPUT=~{name}.bam \
               REFERENCE_SEQUENCE=~{ref.ref_fasta} \
               PAIRED_RUN=true \
               SORT_ORDER="unsorted" \
               IS_BISULFITE_SEQUENCE=false \
               ALIGNED_READS_ONLY=false \
               CLIP_ADAPTERS=false \
               MAX_RECORDS_IN_RAM=2000000 \
               ADD_MATE_CIGAR=true \
               MAX_INSERTIONS_OR_DELETIONS=-1 \
               PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
               PROGRAM_RECORD_ID="bwamem" \
               PROGRAM_GROUP_VERSION="~{bwa_version}" \
               PROGRAM_GROUP_COMMAND_LINE="bwa mem" \
               PROGRAM_GROUP_NAME="bwamem" \
               UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
               ALIGNER_PROPER_PAIR_FLAGS=true \
               UNMAP_CONTAMINANT_READS=true \
               ADD_PG_TAG_TO_READS=false \
               CREATE_INDEX=true
    >>>

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud@sha256:82dd1af86c9e6d4432170133382053525864d8f156a352e18ecf5947542e0b29"
        disks: "local-disk "+disk_size+" HDD"
        memory: "16 GB"
    }

    output {
        File aligned_bam = "~{name}.bam"
        File aligned_bam_index = "~{name}.bai"
        File bwa_sterr_log = "bwa.stderr.log"
    }
}

task MergeBamsByBarcodes {
    input {
        Array[File] aligned_bams
        Array[File] aligned_bam_indecies
        File? picardJar
        String sample_name
    }

    Int disk_size = 10 + ceil(3*size(aligned_bams, "GiB"))

    command <<<
        set -e

        #check that all input bams are the same sample
        for bam in ~{sep=" " aligned_bams}; do
            sample=`samtools view -H $bam | \
                awk '$1=="@RG" {
                        for(i=2; i<=NF; i++) {
                            split($i,split_str,":")
                            if(split_str[1]=="SM") {
                                print split_str[2]
                            }
                        }
                    }' | sort | uniq`

            if [ $sample != ~{sample_name} ]; then
                exit 1
            fi
        done

        java -jar ~{default="/usr/gitc/picard.jar" picardJar} MergeSamFiles I=~{sep=" I=" aligned_bams} O=~{sample_name}.bam CREATE_INDEX=true
    >>>

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud@sha256:82dd1af86c9e6d4432170133382053525864d8f156a352e18ecf5947542e0b29"
        preemptible: 0
        disks: "local-disk "+disk_size+" HDD"
    }

    output {
        File merged_bam = "~{sample_name}.bam"
        File merged_bam_index = "~{sample_name}.bai"
    }
}

task GetSampleName {
    input {
        File bam
    }

    Int disk_size = 10 + ceil(size(bam, "GiB"))
    command <<<
        set -e

        samtools view -H ~{bam} | \
            awk '$1=="@RG" {
                    for(i=2; i<=NF; i++) {
                        split($i,split_str,":")
                        if(split_str[1]=="SM") {
                            print split_str[2]
                        }
                    }
                }' | sort | uniq
    >>>
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud@sha256:82dd1af86c9e6d4432170133382053525864d8f156a352e18ecf5947542e0b29"
    disks: "local-disk "+disk_size+" HDD"
    memory: "1 GB"
  }
  output {
    String sample_name = read_string(stdout())
  }
}

task CountFeaturesJL {
    input {
        File aligned_bam
        File aligned_bam_index
        File feature_file
    }

    Int disk_size = 10 + ceil(10*size(aligned_bam, "GiB"))
    String feature_name = basename(feature_file, "_ALL.gff")
    String sample_name = basename(aligned_bam, ".bam")

    command <<<
        set -e

        /root/tools/samtools sort -n ~{aligned_bam} -o aligned.sam
        mkdir temp
        /root/tools/SAM_to_counts2.sh aligned.sam ~{feature_file} /root/pipeline/shell_scripts temp/ ~{sample_name}_~{feature_name}.counts
    >>>

    runtime {
        docker: "quay.io/ckachuli/microbial_rna_tag_seq:0.1"
        disks: "local-disk "+disk_size+" HDD"
        memory: "16 GB"
    }

    output {
        File counts = "~{sample_name}_~{feature_name}.counts"
        File mets = "~{sample_name}_~{feature_name}.mets"
    }
}

task CountFeaturesFeatureCounts {
    input {
        File aligned_bam
        File aligned_bam_index
        File feature_file
    }

    Int disk_size = 10 + ceil(10*size(aligned_bam, "GiB"))
    String feature_name = basename(feature_file, ".gtf")
    String sample_name = basename(aligned_bam, ".bam")

    command <<<
        set -e

        /usr/local/bin/featureCounts -a ~{feature_file} -t feature -o ~{sample_name}_~{feature_name}.counts -O --fraction -p -s 1 --primary ~{aligned_bam}
        /usr/local/bin/featureCounts -a ~{feature_file} -t feature -o ~{sample_name}_~{feature_name}.reversely_stranded.counts -O --fraction -p -s 2 --primary ~{aligned_bam}
    >>>

    runtime {
            docker: "quay.io/biocontainers/subread:1.6.4--h84994c4_1"
            disks: "local-disk "+disk_size+" HDD"
            memory: "16 GB"
    }

    output {
        File counts = "~{sample_name}_~{feature_name}.counts"
        File counts_rev_stranded = "~{sample_name}_~{feature_name}.reversely_stranded.counts"
        File counts_summary = "~{sample_name}_~{feature_name}.counts.summary"
        File counts_rev_stranded_summary = "~{sample_name}_~{feature_name}.reversely_stranded.counts.summary"
    }
}

task CountFeaturesGATK {
    input {
        File aligned_bam
        File aligned_bam_index
        File feature_file
        File gatk_jar
        File sequence_dict
        String sample_name
    }

    Int disk_size = 10 + ceil(1.2*size(aligned_bam, "GiB"))
    String feature_file_basename = basename(feature_file, ".gff")
    String name = basename(aligned_bam, ".bam")

    command <<<
        set -e

        export GATK_LOCAL_JAR=~{gatk_jar}

        #sort and index feature_file
        /usr/gitc/gatk4/gatk-launch SortGff -F ~{feature_file} -O ~{feature_file_basename}_sorted.gff
        /usr/gitc/gatk4/gatk-launch IndexFeatureFile -F ~{feature_file_basename}_sorted.gff

        /usr/gitc/gatk4/gatk-launch CollectFragmentCounts \
                                    -I ~{aligned_bam} \
                                    -G ~{feature_file_basename}_sorted.gff \
                                    -O ~{name}.counts \
                                    --sequence-dictionary ~{sequence_dict} \
                                    --gene_id_key ID \
                                    --type CDS \
                                    --type ncRNA \
                                    --type rRNA \
                                    --type tRNA \
                                    --name ~{sample_name}
    >>>

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud@sha256:82dd1af86c9e6d4432170133382053525864d8f156a352e18ecf5947542e0b29"
        disks: "local-disk "+disk_size+" HDD"
        memory: "16 GB"
    }

    output {
        File counts = "~{name}.counts"
    }

}

task CountFeaturesFragments {
    input {
        File aligned_bam
        File aligned_bam_index
        File feature_file
        File picard_jar
    }

    Int disk_size = 10 + ceil(10*size(aligned_bam, "GiB"))
    String feature_name = basename(feature_file, ".gtf")
    String sample_name = basename(aligned_bam, ".bam")

    command <<<
        set -e

        java -jar ~{picard_jar} ReadsToFragments I=~{aligned_bam} O=~{sample_name}_fragments.bam
        /root/tools/featureCounts -a ~{feature_file} -t feature -o ~{sample_name}_~{feature_name}.counts -O --fraction -p -s 1 --primary --fracOverlap 0.1 -Q 10 ~{sample_name}_fragments.bam
        /root/tools/featureCounts -a ~{feature_file} -t feature -o ~{sample_name}_~{feature_name}.reversely_stranded.counts -O --fraction -p -s 2 --primary --fracOverlap 0.1 -Q 10 ~{sample_name}_fragments.bam
    >>>

    runtime {
            docker: "quay.io/ckachuli/microbial_rna_tag_seq:0.1"
            disks: "local-disk "+disk_size+" HDD"
            memory: "16 GB"
    }

    output {
        File counts = "~{sample_name}_~{feature_name}.counts"
        File counts_rev_stranded = "~{sample_name}_~{feature_name}.reversely_stranded.counts"
        File counts_summary = "~{sample_name}_~{feature_name}.counts.summary"
        File counts_rev_stranded_summary = "~{sample_name}_~{feature_name}.reversely_stranded.counts.summary"
    }
}

task MarkDuplicates {
    input {
        File aligned_bam
        File aligned_bam_index
        File? picardJar
        Boolean removeDups = false
    }

    Int disk_size = 10 + ceil(2.2*size(aligned_bam, "GiB"))
    String name = basename(aligned_bam, ".bam")

    command <<<
        set -e

        java -Xmx6G -jar ~{default="/usr/gitc/picard.jar" picardJar} MarkDuplicates I=~{aligned_bam} O=~{name}.duplicate_marked.bam \
            M=~{name}.duplicate_marking.metrics CREATE_INDEX=true ~{if(removeDups) then "REMOVE_DUPLICATES=true" else ""}
    >>>

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud@sha256:82dd1af86c9e6d4432170133382053525864d8f156a352e18ecf5947542e0b29"
        preemptible: 0
        disks: "local-disk "+disk_size+" HDD"
        memory: "16 GB"
    }

    output {
        File deduped_bam = "~{name}.duplicate_marked.bam"
        File deduped_bam_index = "~{name}.duplicate_marked.bai"
        File duplicate_marking_metrics = "~{name}.duplicate_marking.metrics"
    }
}

task CollectAlignmentSummaryMetrics {
    input {
        File bam
        File bam_index
        File ref_fasta
        File? picardJar
    }

    Int disk_size = 10 + ceil(2.2*size(bam, "GiB")) + ceil(size(ref_fasta, "GiB"))
    String name = basename(bam, ".bam")

    command <<<
        set -e
        java -jar ~{default="/usr/gitc/picard.jar" picardJar} CollectAlignmentSummaryMetrics R=~{ref_fasta} I=~{bam} O=~{name}_AlignmentSummaryMetrics.txt
    >>>

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud@sha256:82dd1af86c9e6d4432170133382053525864d8f156a352e18ecf5947542e0b29"
        preemptible: 0
        disks: "local-disk "+disk_size+" HDD"
        memory: "16 GB"
    }

    output {
        File alignment_summary_metrics = "~{name}_AlignmentSummaryMetrics.txt"
    }
}

task MeanQualityByCycle {
    input {
        File bam
        File bam_index
        File? picardJar
    }

    Int disk_size = 10 + ceil(2.2*size(bam, "GiB"))
    String name = basename(bam, ".bam")

    command <<<
        set -e
        java -jar ~{default="/usr/gitc/picard.jar" picardJar} MeanQualityByCycle I=~{bam} O=~{name}.mean_quality_by_cycle.txt \
            CHART=~{name}.mean_quality_by_cycle.pdf PF_READS_ONLY=false ALIGNED_READS_ONLY=false
    >>>

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud@sha256:82dd1af86c9e6d4432170133382053525864d8f156a352e18ecf5947542e0b29"
        preemptible: 0
        disks: "local-disk "+disk_size+" HDD"
    }

    output {
        File mean_quality_txt = "~{name}.mean_quality_by_cycle.txt"
        File mean_quality_pdf = "~{name}.mean_quality_by_cycle.pdf"
    }
}

task CollectInsertSizeMetrics {
    input {
        File bam
        File bam_index
        File? picardJar
    }

    Int disk_size = 10 + ceil(1.2*size(bam, "GiB"))
    String name = basename(bam, ".bam")

    command <<<
        set -e
        java -jar ~{default="/usr/gitc/picard.jar" picardJar} CollectInsertSizeMetrics I=~{bam} O=~{name}.insert_size.metrics \
            H=~{name}.insert_size_histogram.pdf M=0.5
    >>>

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud@sha256:82dd1af86c9e6d4432170133382053525864d8f156a352e18ecf5947542e0b29"
        preemptible: 0
        disks: "local-disk "+disk_size+" HDD"
    }

    output {
        File insert_size_metrics = "~{name}.insert_size.metrics"
        File insert_size_histogram = "~{name}.insert_size_histogram.pdf"
    }
}

task CollectRnaSeqMetrics {
    input {
        File bam
        File bam_index
        File gff
        File? picardJar
    }

    Int disk_size = 10 + ceil(1.2*size(bam, "GiB"))
    String name = basename(bam, ".bam")

    command <<<
        set -e
        java -jar ~{default="/usr/gitc/picard.jar" picardJar} CollectRnaSeqMetrics I=~{bam} ANNOTATION_FILE=~{gff} STRAND=SECOND_READ_TRANSCRIPTION_STRAND O=~{name}.rna_seq_metrics
    >>>

    runtime {
        docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud@sha256:82dd1af86c9e6d4432170133382053525864d8f156a352e18ecf5947542e0b29"
        preemptible: 0
        disks: "local-disk "+disk_size+" HDD"
    }

    output {
        File rna_seq_metrics = "~{name}.rna_seq_metrics"
    }
}

task BamToTDF {
    input {
        File bam
        File bam_index
    }

    Int disk_size = 10 + ceil(3.2*size(bam, "GiB"))
    String name = basename(bam)
    command <<<
        java -jar /root/tools/genometools-24.jar bam2tdf ~{bam}
        mv ~{bam}.tdf .
    >>>

    runtime {
        docker: "quay.io/ckachuli/microbial_rna_tag_seq:0.1"
        disks: "local-disk "+disk_size+" HDD"
        memory: "16 GB"
    }

    output {
        File tdf = "~{name}.tdf"
    }
}

task CombineGeneCounts {
    input {
        Array[File] count_files
        String accension
        Array[String] sample_names
    }

    Int disk_size = 10 + ceil(2.2*size(count_files, "GiB"))

    command <<<
        set -e

        echo "Geneid" > ~{accension}_counts.tsv
        cut ~{count_files[0]} -f 1 >> ~{accension}_counts.tsv

        for file in ~{sep=" " count_files}; do
            name=(`basename $file`)
            for sample in ~{sep=" " sample_names}; do
                if [[ $name == "$sample"* ]]; then
                    echo "$sample" > temp_counts.tsv
                    cut $file -f 3 >> temp_counts.tsv
                    paste -d"\t" ~{accension}_counts.tsv temp_counts.tsv > temp_counts2.tsv
                    cat temp_counts2.tsv > ~{accension}_counts.tsv
                    break
                fi
            done
        done
    >>>

    runtime {
        docker: "ubuntu:16.04"
        preemptible: 0
        disks: "local-disk "+disk_size+" HDD"
    }

    output {
        File combined_gene_counts = "~{accension}_counts.tsv"
    }
}

task CombineGeneCountsFeatureCounts {
    input {
        Array[File] count_files
        String outName
    }

    Int disk_size = 10 + ceil(2.2*size(count_files, "GiB"))

    command <<<
        set -e

        Rscript -<< "EOF"
        library(data.table)
        file_list <- c("~{sep='","' count_files}")
        d_list <- lapply(file_list, fread)
        d_list <- lapply(d_list, function(dt) {names(dt)[7]<-tools::file_path_sans_ext(basename(names(dt)[7]));dt})
        lapply(d_list,setkey,Geneid,Chr,Start,End,Strand,Length)
        merged <- Reduce(function(...) merge(..., all = TRUE), d_list)
        fwrite(merged, "~{outName}")
        EOF
    >>>

    runtime {
        docker: "rocker/tidyverse:3.6"
        disks: "local-disk "+disk_size+" HDD"
        memory: "16 GB"
    }

    output {
        File combined = "~{outName}"
    }
}


task CombineGeneCountsGATK {
    input {
        Array[File] count_files
        String outName
    }

    Int disk_size = 10 + ceil(2.2*size(count_files, "GiB"))

    command <<<
        set -e

        Rscript -<< "EOF"
        library(data.table)
        file_list <- c("~{sep='","' count_files}")
        d_list <- lapply(file_list, fread)
        lapply(d_list,setkey,gene_id,contig,start,stop,strand,sense_antisense)
        merged <- Reduce(function(...) merge(..., all = TRUE), d_list)
        fwrite(merged, "~{outName}")
        EOF
    >>>

    runtime {
        docker: "rocker/tidyverse:3.6"
        disks: "local-disk "+disk_size+" HDD"
        memory: "16 GB"
    }

    output {
        File combined = "~{outName}"
    }
}

task CountsToFPKMFeatureCounts {
    input {
        File counts
        File alignmentMetrics
    }

    Int disk_size = 10 + ceil(4.2*size(counts, "GiB"))
    String fpkm_name = sub(basename(counts), "counts", "fpkm")

    command <<<
        set -e

        Rscript -<< "EOF"
        library(data.table)
        d_counts <- fread("~{counts}")
        d_alignment <- fread("~{alignmentMetrics}")
        n_fragments <- d_alignment[CATEGORY=="FIRST_OF_PAIR", PF_READS_ALIGNED]
        d_fpkm <- data.table(d_counts)
        colCounts <- names(d_fpkm)[7]
        counts <- d_fpkm[,..colCounts]
        d_fpkm[,7:=counts/(Length*1e-3*n_fragments*1e-6)]
        fwrite(d_fpkm, "~{fpkm_name}")
        EOF
    >>>

    runtime {
        docker: "rocker/tidyverse:3.6"
        disks: "local-disk "+disk_size+" HDD"
        memory: "16 GB"
    }

    output {
        File fpkm = "~{fpkm_name}"
    }
}

task CountsToFPKMGATK {
    input {
        File counts
        File alignmentMetrics
    }

    Int disk_size = 10 + ceil(4.2*size(counts, "GiB"))
    String fpkm_name = sub(basename(counts), "counts", "fpkm")

    command <<<
        set -e

        Rscript -<< "EOF"
        library(data.table)
        d_counts <- fread("~{counts}", skip="gene_id")
        d_alignment <- fread("~{alignmentMetrics}")
        n_fragments <- d_alignment[CATEGORY=="FIRST_OF_PAIR", PF_READS_ALIGNED]
        d_fpkm <- data.table(d_counts)
        colCounts <- names(d_fpkm)[7]
        counts <- d_fpkm[,..colCounts]
        names(d_fpkm)[7] <- sub("_counts", "_fpkm", colCounts)
        d_fpkm[,7:=counts/((stop-start+1)*1e-3*n_fragments*1e-6)]
        fwrite(d_fpkm, "~{fpkm_name}")
        EOF
    >>>

    runtime {
        docker: "rocker/tidyverse:3.6"
        disks: "local-disk "+disk_size+" HDD"
        memory: "16 GB"
    }

    output {
        File fpkm = "~{fpkm_name}"
    }
}

task CountsToFPKM {
    input {
        File counts
    }

    Int disk_size = 10 + ceil(4.2*size(counts, "GiB"))
    String fpkm_name = sub(basename(counts), "counts", "fpkm")

    command <<<
        set -e

        mkdir temp
        /root/tools/counts_to_FPKMO.sh ~{counts} ~{fpkm_name} temp
    >>>

    runtime {
        docker: "quay.io/ckachuli/microbial_rna_tag_seq:0.1"
        disks: "local-disk "+disk_size+" HDD"
        memory: "16 GB"
    }

    output {
        File fpkm = "~{fpkm_name}"
    }
}

task SampleCorrelations {
    input {
        File counts_file
    }

    Int disk_size = 10 + ceil(2*size(counts_file, "GiB"))
    String corr_file_name = sub(basename(counts_file), "counts.tsv", "corr.txt")

    command <<<
        set -e

        Rscript -<< "EOF"
        FPKMO= read.table(file = "~{counts_file}", header=T, row.names = 1, sep = "\t", quote = "", comment.char = "")
        sorted_FPKMO <- FPKMO[ , order(names(FPKMO))]

        log_FPKMO <- log10(sorted_FPKMO+1)

        corr <- cor(log_FPKMO)

        write.table(corr, file = "~{corr_file_name}",  sep = "\t", quote = F, row.names = T)
        EOF
    >>>

    runtime {
        docker: "quay.io/ckachuli/microbial_rna_tag_seq:0.1"
        disks: "local-disk "+disk_size+" HDD"
        memory: "16 GB"
    }

    output {
        File corr_file= "~{corr_file_name}"
    }
}

task SampleCorrelationsFeatureCounts {
    input {
        File fpkm_sense
        File fpkm_antisense
    }

    Int disk_size = 10 + ceil(2*size(fpkm_sense, "GiB") + 2*size(fpkm_antisense, "GiB"))
    String corr_file_name = sub(basename(fpkm_sense), "sense.counts", "corr.txt")

    command <<<
        set -e

        Rscript -<< "EOF"
        library(data.table)
        fpkm_sense <- fread("~{fpkm_sense}")
        fpkm_antisense <- fread("~{fpkm_antisense}")

        fpkm <- rbind(fpkm_sense, fpkm_antisense)
        fpkm <- fpkm[,7:ncol(fpkm)]

        log_fpkm <- log10(fpkm+1)

        corr <- cor(log_fpkm)

        fwrite(corr, "~{corr_file_name}")
        EOF
    >>>

    runtime {
        docker: "rocker/tidyverse"
        disks: "local-disk "+disk_size+" HDD"
        memory: "16 GB"
    }

    output {
        File corr_file= "~{corr_file_name}"
    }
}

task SampleCorrelationsGATK {
    input {
        File fpkm
    }

    Int disk_size = 10 + ceil(2*size(fpkm, "GiB"))
    String corr_file_name = sub(basename(fpkm), ".fpkm", ".corr.txt")

    command <<<
        set -e

        Rscript -<< "EOF"
        library(data.table)
        fpkm<- fread("~{fpkm}")
        fpkm <- fpkm[,7:ncol(fpkm)]

        log_fpkm <- log10(fpkm+1)

        corr <- cor(log_fpkm)

        fwrite(corr, "~{corr_file_name}")
        EOF
    >>>

    runtime {
        docker: "rocker/tidyverse"
        disks: "local-disk "+disk_size+" HDD"
        memory: "16 GB"
    }

    output {
        File corr_file= "~{corr_file_name}"
    }
}


task Picard_metrics_parse {
    input {
        Array[File] metrics
        String project
    }

    command <<<
        set -e

        mkdir metrics_links
        ln -s ~{sep=" " metrics} metrics_links
        /root/tools/PICARD_metrics_parse.sh metrics_links ~{project}_AlignmentSummaryMetrics.txt
    >>>

    runtime {
        docker: "quay.io/ckachuli/microbial_rna_tag_seq:0.1"
        disks: "local-disk 25 HDD"
        memory: "16 GB"
    }

    output {
        File alignment_summary_metrics = "~{project}_AlignmentSummaryMetrics.txt"
    }
}

task RPG_metrics {
    input {
        Array[File] counts
        Array[File] mets
        String project_id
    }

    command <<<
        set -e

        mkdir -p temp/patho_result/~{project_id}
        ln -s ~{sep=" " counts} temp/patho_result/~{project_id}
        ln -s ~{sep=" " mets} temp/patho_result/~{project_id}
        /root/tools/RPG_metrics5.sh -p ~{project_id} -t temp -r results
    >>>

    runtime {
        docker: "quay.io/ckachuli/microbial_rna_tag_seq:0.1"
        disks: "local-disk 25 HDD"
        memory: "16 GB"
    }

    output {
        File metrics = "results/~{project_id}/~{project_id}_metrics.txt"
    }
}

