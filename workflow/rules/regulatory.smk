# --- Cis-regulatory / TFBS layer -------------------------------------------
# Extract promoters for each AMP gene and its processing-enzyme gene, scan for
# TF binding motifs, and report the motifs SHARED between co-regulated pairs
# (a candidate common regulatory signature) plus a genome-mapped TFBS track.

rule extract_promoters:
    input:
        map    = OUT + "/homology/peptide_to_precursor.tsv",
        enzyme = OUT + "/methods/method4_enzyme.tsv",
        gff    = RES + "/GRCh38.annotation.gff3",
        fa     = RES + "/GRCh38.primary_assembly.fa",
        fai    = RES + "/GRCh38.primary_assembly.fa.fai",
    output:
        fa    = OUT + "/regulatory/promoters.fa",
        pairs = OUT + "/regulatory/pairs.tsv",
    conda:
        "../envs/core.yaml"
    params:
        up   = config["regulatory"]["promoter_upstream"],
        down = config["regulatory"]["promoter_downstream"],
    shell:
        r"""
        mkdir -p {OUT}/regulatory
        python workflow/scripts/extract_promoters.py \
            --map {input.map} --enzyme-tsv {input.enzyme} \
            --gff {input.gff} --genome {input.fa} \
            --upstream {params.up} --downstream {params.down} \
            --out-fasta {output.fa} --out-pairs {output.pairs}
        """

rule tfbs_scan:
    input:
        prom   = OUT + "/regulatory/promoters.fa",
        motifs = config["tf_motifs"],
        pairs  = OUT + "/regulatory/pairs.tsv",
    output:
        hits   = OUT + "/regulatory/tfbs_hits.tsv",
        shared = OUT + "/regulatory/shared_signatures.tsv",
        bed    = OUT + "/regulatory/tfbs.bed",
    conda:
        "../envs/core.yaml"
    params:
        mm = config["regulatory"]["tfbs_max_mismatch"],
    shell:
        r"""
        python workflow/scripts/shared_tfbs.py \
            --promoters {input.prom} --motifs {input.motifs} \
            --pairs {input.pairs} --max-mismatch {params.mm} \
            --out-hits {output.hits} --out-shared {output.shared} \
            --out-bed {output.bed}
        """
