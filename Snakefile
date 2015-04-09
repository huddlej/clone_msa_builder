#
# Define configuration.
#

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
        "dotplots.pdf"
    params: sge_opts=""

rule merge_dotplots:
    input: expand("dotplots/{clone}.pdf", clone=CLONES)
    output: "dotplots.pdf"
    params: sge_opts=""
    shell: "gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile={output} {input}"

rule align_regions_by_species:
    input: "query_regions_by_species/{species}.fasta"
    output: "multiple_sequence_alignments_by_species/{species}.fasta"
    params: sge_opts="-l mfree=1G -pe serial 4", threads="4"
    shell: "mafft-linsi --thread {params.threads} {input} > {output}"

rule combine_query_regions_by_species:
    input: _get_clones_for_species
    output: "query_regions_by_species/{species}.fasta"
    params: sge_opts=""
    shell: "sed 's/>/>{wildcards.species}_/' {input} > {output}"

rule extract_query_region_from_clone:
    input: fasta="original_sequences/{clone}.fasta", regions="query_placements/{clone}.bed"
    output: "query_regions/{clone}.fasta"
    params: sge_opts=""
    shell: "bedtools getfasta -fi {input.fasta} -bed {input.regions} -fo {output} -s"

rule convert_psl_to_merged_bed:
    input: alignment="psl_alignments/{clone}.psl", clone_index="original_sequences/{clone}.fasta.fai"
    output: "query_placements/{clone}.bed"
    params: sge_opts="", merge_distance="10000", slop="5000"
    shell: "pslToBed {input.alignment} /dev/stdout | bedtools merge -i stdin -d {params.merge_distance} -s | bedtools slop -i stdin -g {input.clone_index} -b {params.slop} > {output}"

rule convert_lav_to_psl:
    input: "lav_alignments/{clone}.lav"
    output: "psl_alignments/{clone}.psl"
    params: sge_opts=""
    shell: "lavToPsl {input} {output}"

rule plot_dotplot_for_query_in_clone:
    input: "dotplots/{clone}.dotplot"
    output: "dotplots/{clone}.pdf"
    params: sge_opts=""
    shell: "Rscript ~jlhudd/fasta_tools/lastz_dotplot.R {input} {output}"

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
    shell: "rsync {input} {output}"

rule clean:
    shell: "rm -rf query_* psl_alignments/ original_sequences/ multiple_sequence_alignments_by_species/ lav_alignments/ dotplots/"
