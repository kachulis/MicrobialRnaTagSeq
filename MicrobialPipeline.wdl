version 1.0

import "tasks.wdl" as tasks
import "structs.wdl"

workflow MicrobialPipeline {
    input {
        Array[FastQSet] fastq_sets
        File sample_sheet
        String project_id
        Array[File]+ gffs
        Array[File]+ fastas
        Array[String]+ microbeNames
    }

#    call tasks.CreateAccensionFiles {
#        input:
#            gffs = gffs,
#            fnas = fnas,
#            accs = accs,
#            add3 = add3,
#            add5 = add5
#    }
#    scatter (sample in samples) {
#        String input_sample_names = sample.name
#        String input_sample_barcodes = sample.barcode
#    }

#    call tasks.BuildBWAIndex {
#        input:
#            fasta = CreateAccensionFiles.out_fna
#    }

    call tasks.GetBwaVersion

    call tasks.BuildCombinedName {
        input:
            names = microbeNames
    }

    call tasks.BuildReferenceResources {
        input:
            fastas = fastas,
            gffs = gffs,
            combined_name = BuildCombinedName.combined_name

    }

    ReferenceFasta ref = object {
                        ref_dict : BuildReferenceResources.combined_dict,
                        ref_fasta : BuildReferenceResources.combined_fasta,
                        ref_fasta_index : BuildReferenceResources.combined_index,
                        ref_sa : BuildReferenceResources.combined_sa,
                        ref_amb : BuildReferenceResources.combined_amb,
                        ref_bwt : BuildReferenceResources.combined_bwt,
                        ref_ann : BuildReferenceResources.combined_ann,
                        ref_pac : BuildReferenceResources.combined_pac,
                                }

    scatter (fastq_set in fastq_sets) {
        call tasks.SplitBarcodes {
            input:
                fastq_set = fastq_set,
                sample_sheet = sample_sheet
        }

        scatter (unmapped_bam in SplitBarcodes.sample_unmapped_bams) {
            call tasks.AlignPairedReadsBWABacktrack {
                input:
                    unmapped_bam = unmapped_bam,
                    ref = ref,
                    bwa_version = GetBwaVersion.bwa_version
            }
        }
    }

    Array[Array[File]] bams_to_merge_arrays = transpose(AlignPairedReadsBWABacktrack.aligned_bam)
    Array[Array[File]] bais_to_merge_arrays = transpose(AlignPairedReadsBWABacktrack.aligned_bam_index)

    scatter (i in range(length(bams_to_merge_arrays))) {
        Array[File] bams_to_merge = bams_to_merge_arrays[i]
        Array[File] bais_to_merge = bais_to_merge_arrays[i]

        call tasks.GetSampleName {
            input:
                bam = bams_to_merge[0]
        }
        call tasks.MergeBamsByBarcodes {
            input:
                aligned_bams = bams_to_merge,
                aligned_bam_indecies = bais_to_merge,
                sample_name = GetSampleName.sample_name
        }

#        call tasks.CountFeaturesJL {
#            input:
#                aligned_bam = MergeBamsByBarcodes.merged_bam,
#                aligned_bam_index = MergeBamsByBarcodes.merged_bam_index,
#                feature_file = gff
#        }

#        call tasks.CountFeaturesFeatureCounts {
#            input:
#                aligned_bam = MergeBamsByBarcodes.merged_bam,
#                aligned_bam_index = MergeBamsByBarcodes.merged_bam_index,
#                feature_file = gff
#        }

#        call tasks.CountFeaturesFragments {
#            input:
#                aligned_bam = MergeBamsByBarcodes.merged_bam,
#                aligned_bam_index = MergeBamsByBarcodes.merged_bam_index,
#                feature_file = gtf,
#                picard_jar = "gs://broad-dsde-methods-ckachulis/MicrobialRNASeq/jars/picard_ReadsToFragments.jar"
#        }

        call tasks.CountFeaturesGATK {
            input:
                aligned_bam = MergeBamsByBarcodes.merged_bam,
                aligned_bam_index = MergeBamsByBarcodes.merged_bam_index,
                feature_file = BuildReferenceResources.combined_gff,
                gatk_jar = "gs://broad-dsde-methods-ckachulis/jars/gatk_collect_fragment_counts_microbial_v4.jar",
                sample_name = GetSampleName.sample_name,
                sequence_dict = ref.ref_dict
        }

        call tasks.CountsToFPKMGATK {
                    input:
                        counts = CountFeaturesGATK.counts,
                        alignmentMetrics = CollectAlignmentSummaryMetrics.alignment_summary_metrics
        }
#        call tasks.CountsToFPKMFeatureCounts as CountsToFPKMSense {
#            input:
#                counts = CountFeaturesFragments.counts,
#                alignmentMetrics = CollectAlignmentSummaryMetrics.alignment_summary_metrics
#        }

#        call tasks.CountsToFPKMFeatureCounts as CountsToFPKMAntiSense {
#            input:
#                counts = CountFeaturesFragments.counts_rev_stranded,
#                alignmentMetrics = CollectAlignmentSummaryMetrics.alignment_summary_metrics
#        }

        call tasks.BamToTDF {
            input:
                bam = MergeBamsByBarcodes.merged_bam,
                bam_index = MergeBamsByBarcodes.merged_bam_index
        }

        call tasks.CollectAlignmentSummaryMetrics {
            input:
                bam = MergeBamsByBarcodes.merged_bam,
                bam_index = MergeBamsByBarcodes.merged_bam_index,
                ref_fasta = ref.ref_fasta
        }

        call tasks.MeanQualityByCycle {
            input:
                bam = MergeBamsByBarcodes.merged_bam,
                bam_index = MergeBamsByBarcodes.merged_bam_index
        }

        call tasks.CollectInsertSizeMetrics {
            input:
                bam = MergeBamsByBarcodes.merged_bam,
                bam_index = MergeBamsByBarcodes.merged_bam_index
        }

        call tasks.CollectRnaSeqMetrics {
            input:
                bam = MergeBamsByBarcodes.merged_bam,
                bam_index = MergeBamsByBarcodes.merged_bam_index,
                gff = BuildReferenceResources.combined_gff,
                picardJar = "gs://broad-dsde-methods-ckachulis/jars/picard_CollectRNAMetrics.jar"
        }
    }

#    call tasks.CombineGeneCounts {
#        input:
#            count_files = CountFeaturesFeatureCounts.counts,
#            accension = accs[0],
#            sample_names = GetSampleName.sample_name
#    }

#    call tasks.CombineGeneCountsFeatureCounts as CombineGeneCountsSense {
#        input:
#            count_files = CountFeaturesFragments.counts,
#            outName = "Sense.counts"
#    }

    call tasks.CombineGeneCountsGATK as CombineCounts {
        input:
            count_files = CountFeaturesGATK.counts,
            outName = "All.counts"
    }

#    call tasks.CombineGeneCountsFeatureCounts as CombineGeneCountsAntiSense {
#         input:
#             count_files = CountFeaturesFragments.counts_rev_stranded,
#             outName = "AntiSense.counts"
#    }

#    call tasks.CombineGeneCountsFeatureCounts as CombineFPKMSense {
#        input:
#            count_files = CountsToFPKMSense.fpkm,
#            outName = "Sense.fpkm"
#    }

    call tasks.CombineGeneCountsGATK as CombineFPKM {
        input:
            count_files = CountsToFPKMGATK.fpkm,
            outName = "All.fpkm"
    }

#    call tasks.CombineGeneCountsFeatureCounts as CombineFPKMAntiSense {
#         input:
#             count_files = CountsToFPKMAntiSense.fpkm,
#             outName = "AntiSense.fpkm"
#    }

#    call tasks.CountsToFPKM {
#        input:
#            counts = CombineGeneCounts.combined_gene_counts
#    }

#    call tasks.SampleCorrelations {
#        input:
#            counts_file = CountsToFPKM.fpkm
#    }

#    call tasks.SampleCorrelationsFeatureCounts {
#        input:
#            fpkm_sense = CombineFPKMSense.combined,
#            fpkm_antisense = CombineFPKMAntiSense.combined
#    }

        call tasks.SampleCorrelationsGATK {
            input:
                fpkm = CombineFPKM.combined
        }

    call tasks.Picard_metrics_parse {
        input:
            metrics = CollectAlignmentSummaryMetrics.alignment_summary_metrics,
            project = project_id
    }

#    call tasks.RPG_metrics {
#        input:
#            counts = CountFeaturesJL.counts,
#            mets = CountFeaturesJL.mets,
#            project_id = project_id
#    }
}






