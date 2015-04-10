#
# Define configuration.
#
import os

configfile: "config.json"
SPECIES = config["species"]
CLONES = sorted(config["clone_paths"].keys())

#
# Define helper functions.
#

def _get_sequence_for_clone(wildcards):
    """
    Get the absolute path to the FASTA sequence for the given clone name.
    """
    return config["clone_paths"][wildcards.clone]

def _get_clones_for_species(wildcards):
    """
    Get a list of clone names for the given species.
    """
    return ["query_regions/%s.fasta" % clone for clone in config["clones_by_species"][wildcards.species]]

#
# Define rules.
#

rule all:
    input:
        expand("multiple_sequence_alignments_by_species/{species}.fasta", species=SPECIES),
        expand("pairwise_identity/{species}.pdf", species=SPECIES),
        "dotplots.pdf"
    params: sge_opts=""

rule merge_dotplots:
    input: expand("dotplots/{clone}.pdf", clone=CLONES)
    output: "dotplots.pdf"
    params: sge_opts=""
    shell: "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile={output} `for file in {input}; do if [[ -e $file ]]; then echo $file; fi; done`"

rule plot_pairwise_identity_by_species:
    input: "multiple_sequence_alignments_by_species/{species}.fasta"
    output: "pairwise_identity/{species}.pdf"
    params: sge_opts=""
    shell: "Rscript plot_pairwise_identity.R `dirname {input}` {output}"

# rule mask_repeats_in_alignment:
#     input: alignment="multiple_sequence_alignments_by_species/{species}.fasta", repeats=config["repeat_sequences"]
#     output: "masked_alignments_by_species/{species}.consensus.fa"
#     params: sge_opts="", window="1000", window_slide="100"
#     shell: "mkdir -p `dirname {output}`; cd `dirname {output}`; ln -s ../{input.alignment} .; mam `basename {input.alignment}` -sw={params.window_slide} -ww={params.window} -program=crossmatch -fasta=on -exonfile={input.repeats} -slider=on -pc=c -identity=on -consensus=on -alnstats=on"

rule align_regions_by_species:
    input: "query_regions_by_species/{species}.fasta"
    output: "multiple_sequence_alignments_by_species/{species}.fasta"
    params: sge_opts="-l mfree=10G -pe serial 2", threads="2"
    shell: "mafft --auto --thread {params.threads} {input} > {output}"

rule combine_query_regions_by_species:
    input: _get_clones_for_species
    output: "query_regions_by_species/{species}.fasta"
    params: sge_opts=""
    shell: "sed 's/>/>{wildcards.species}_/' {input} > {output}"

rule extract_query_region_from_clone:
    input: fasta="original_sequences/{clone}.fasta", regions="merged_query_placements/{clone}.bed"
    output: "query_regions/{clone}.fasta"
    params: sge_opts=""
    shell: "bedtools getfasta -fi {input.fasta} -bed {input.regions} -fo {output} -s"

rule convert_bed_to_merged_bed:
    input: bed="query_placements/{clone}.bed", clone_index="original_sequences/{clone}.fasta.fai"
    output: "merged_query_placements/{clone}.bed"
    params: sge_opts="", merge_distance="10000", slop="5000"
    run:
        if os.stat(input["bed"]).st_size > 0:
            shell("bedtools merge -i {input.bed} -d {params.merge_distance} -s -c 4,5,6 -o distinct,mean,distinct | bedtools slop -i stdin -g {input.clone_index} -b {params.slop} > {output}")
        else:
            shell("touch {output}")

rule convert_psl_to_merged_bed:
    input: "psl_alignments/{clone}.psl"
    output: "query_placements/{clone}.bed"
    params: sge_opts=""
    shell: "pslToBed {input} {output}"

rule convert_lav_to_psl:
    input: "lav_alignments/{clone}.lav"
    output: "psl_alignments/{clone}.psl"
    params: sge_opts=""
    shell: "lavToPsl {input} {output}"

rule plot_dotplot_for_query_in_clone:
    input: "dotplots/{clone}.dotplot"
    output: "dotplots/{clone}.pdf"
    params: sge_opts=""
    run:
        try:
            shell("Rscript ~jlhudd/fasta_tools/lastz_dotplot.R {input} {output}")
        except:
            shell("touch {output}")

rule find_query_in_clone:
    input: clone="original_sequences/{clone}.fasta", query=config["query_sequence"]
    output: alignment="lav_alignments/{clone}.lav", dotplot="dotplots/{clone}.dotplot"
    params: sge_opts=""
    shell: "lastz {input.clone} {input.query} --gfextend --chain --gapped --notransition --format=lav --rdotplot={output.dotplot} > {output.alignment}"

rule index_clone:
    input: "original_sequences/{clone}.fasta"
    output: "original_sequences/{clone}.fasta.fai"
    params: sge_opts=""
    shell: "samtools faidx {input}"

rule get_clone:
    input: _get_sequence_for_clone
    output: "original_sequences/{clone}.fasta"
    params: sge_opts=""
    shell: "cat {input} > {output}"

rule clean:
    shell: "rm -rf query_* psl_alignments/ original_sequences/ multiple_sequence_alignments_by_species/ lav_alignments/ dotplots/ masked_alignments_by_species/ pairwise_identity/ dotplots.pdf merged_query_placements/"
