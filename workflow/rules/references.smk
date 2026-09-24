# --- Reference acquisition -------------------------------------------------
# GRCh38 genome, Ensembl GFF3, and the UniProt human reference proteome.
# All downloads are cached under resources/ (gitignored).

rule download_genome:
    output:
        fa = RES + "/GRCh38.primary_assembly.fa",
    params:
        url = config["reference"]["genome_fasta_url"],
    threads: 1
    shell:
        r"""
        mkdir -p {RES}
        wget -q -O {output.fa}.gz "{params.url}"
        gunzip -f {output.fa}.gz
        """

rule index_genome:
    input:
        fa = RES + "/GRCh38.primary_assembly.fa",
    output:
        fai = RES + "/GRCh38.primary_assembly.fa.fai",
    conda:
        "../envs/core.yaml"
    shell:
        "samtools faidx {input.fa}"

rule download_gff3:
    output:
        gff = RES + "/GRCh38.annotation.gff3",
    params:
        url = config["reference"]["gff3_url"],
    shell:
        r"""
        mkdir -p {RES}
        wget -q -O {output.gff}.gz "{params.url}"
        gunzip -f {output.gff}.gz
        """

rule download_uniprot:
    output:
        fa = RES + "/uniprot_human.fasta",
    params:
        url = config["reference"]["uniprot_human_url"],
    shell:
        r"""
        mkdir -p {RES}
        wget -q -O {output.fa}.gz "{params.url}"
        gunzip -f {output.fa}.gz
        """

rule diamond_makedb:
    input:
        fa = RES + "/uniprot_human.fasta",
    output:
        db = RES + "/uniprot_human.dmnd",
    conda:
        "../envs/core.yaml"
    threads: config["params"]["threads"]
    shell:
        "diamond makedb --in {input.fa} -d {RES}/uniprot_human -p {threads}"
