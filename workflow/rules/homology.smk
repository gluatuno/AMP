# --- Step 1: mature peptide -> full-length precursor protein ---------------
# DIAMOND blastp of each mature peptide against the human proteome.
# full_sseq gives us the complete precursor sequence in one shot.

rule diamond_blastp:
    input:
        q  = OUT + "/prep/peptides.clean.fa",
        db = RES + "/uniprot_human.dmnd",
    output:
        tsv = OUT + "/homology/peptides_vs_human.tsv",
    conda:
        "../envs/core.yaml"
    threads: config["params"]["threads"]
    params:
        evalue = config["params"]["diamond_evalue"],
    shell:
        r"""
        mkdir -p {OUT}/homology
        diamond blastp \
            -q {input.q} -d {input.db} -p {threads} \
            --evalue {params.evalue} --max-target-seqs 5 --very-sensitive \
            --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore stitle full_sseq \
            -o {output.tsv}
        """

rule best_precursor:
    # Pick the best precursor per peptide, apply identity/coverage filters,
    # and emit (a) a mapping table and (b) the precursor protein FASTA.
    input:
        tsv = OUT + "/homology/peptides_vs_human.tsv",
    output:
        map = OUT + "/homology/peptide_to_precursor.tsv",
        fa  = OUT + "/homology/precursors.faa",
    conda:
        "../envs/core.yaml"
    params:
        min_pident = config["params"]["diamond_min_pident"],
        min_cov    = config["params"]["min_query_cov"],
    shell:
        r"""
        python workflow/scripts/best_hit_precursor.py \
            --blast {input.tsv} \
            --min-pident {params.min_pident} \
            --min-cov {params.min_cov} \
            --out-map {output.map} \
            --out-fasta {output.fa}
        """
