# --- Step 2: precursor protein -> genome (splice-aware) --------------------
# miniprot aligns each precursor to GRCh38 and emits a GFF with the CDS/exon
# structure we need to convert amino-acid coordinates into genomic intervals.

rule miniprot_index:
    input:
        fa = RES + "/GRCh38.primary_assembly.fa",
    output:
        mpi = RES + "/GRCh38.mpi",
    conda:
        "../envs/core.yaml"
    threads: config["params"]["threads"]
    shell:
        "miniprot -t {threads} -d {output.mpi} {input.fa}"

rule miniprot_align:
    input:
        mpi = RES + "/GRCh38.mpi",
        faa = OUT + "/homology/precursors.faa",
    output:
        gff = OUT + "/genome_map/precursors.miniprot.gff",
    conda:
        "../envs/core.yaml"
    threads: config["params"]["threads"]
    shell:
        r"""
        mkdir -p {OUT}/genome_map
        miniprot -t {threads} --gff {input.mpi} {input.faa} > {output.gff}
        """
