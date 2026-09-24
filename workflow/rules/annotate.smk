# --- Step 3: annotate precursor feature blocks -----------------------------
# deepsig predicts the signal peptide; define_regions.py then locates the
# mature peptide inside its precursor and derives the pro-region between them.

rule signalp:
    input:
        faa = OUT + "/homology/precursors.faa",
    output:
        tsv = OUT + "/annotate/signalp.tsv",
    conda:
        "../envs/annotate.yaml"
    threads: config["params"]["threads"]
    shell:
        r"""
        mkdir -p {OUT}/annotate
        deepsig -f {input.faa} -o {output.tsv} -k euk -t {threads}
        """

rule define_regions:
    # Produce per-precursor amino-acid coordinates for signal / pro / mature.
    input:
        precursors = OUT + "/homology/precursors.faa",
        peptides   = OUT + "/prep/peptides.clean.fa",
        map        = OUT + "/homology/peptide_to_precursor.tsv",
        signalp    = OUT + "/annotate/signalp.tsv",
    output:
        regions = OUT + "/annotate/regions_aa.tsv",
    conda:
        "../envs/core.yaml"
    params:
        motifs = ",".join(config["params"]["dibasic_motifs"]),
    shell:
        r"""
        python workflow/scripts/define_regions.py \
            --precursors {input.precursors} \
            --peptides {input.peptides} \
            --map {input.map} \
            --signalp {input.signalp} \
            --dibasic {params.motifs} \
            --out {output.regions}
        """
