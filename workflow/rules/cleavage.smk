# --- Cleavage-site signatures: discovery, prediction, genome annotation ----
# Finds protease "cutting site" signatures encoded in the precursor CDS, uses
# them to predict signal/pro/mature boundaries WITHOUT a known mature peptide,
# and projects them onto the genome as an annotatable cut-site track.

rule scan_cleavage:
    # Knowledge-based scan of precursors with the curated motif library.
    input:
        faa    = OUT + "/homology/precursors.faa",
        motifs = config["protease_motifs"],
    output:
        tsv = OUT + "/cleavage/candidate_sites.tsv",
    conda:
        "../envs/core.yaml"
    shell:
        r"""
        mkdir -p {OUT}/cleavage
        python workflow/scripts/scan_cleavage.py \
            --precursors {input.faa} --motifs {input.motifs} --out {output.tsv}
        """

rule discover_cutsites:
    # De-novo motif discovery from the known junctions (PFM always; STREME/MEME
    # when the meme env is used).
    input:
        faa     = OUT + "/homology/precursors.faa",
        regions = OUT + "/annotate/regions_aa.tsv",
    output:
        report = OUT + "/cleavage/denovo/motif_report.tsv",
    conda:
        "../envs/cleavage.yaml"
    shell:
        r"""
        python workflow/scripts/discover_cutsites.py \
            --precursors {input.faa} --regions {input.regions} \
            --outdir {OUT}/cleavage/denovo
        """

rule predict_boundaries:
    # Signature-only boundary prediction + precision/recall vs known boundaries.
    input:
        faa       = OUT + "/homology/precursors.faa",
        cand      = OUT + "/cleavage/candidate_sites.tsv",
        signalp   = OUT + "/annotate/signalp.tsv",
        truth     = OUT + "/annotate/regions_aa.tsv",
    output:
        regions  = OUT + "/cleavage/predicted_regions.tsv",
        accuracy = OUT + "/cleavage/boundary_accuracy.tsv",
    conda:
        "../envs/core.yaml"
    params:
        min_len = config["params"]["min_mature_len"],
        max_len = config["params"]["max_mature_len"],
        tol     = config["params"]["boundary_tol"],
    shell:
        r"""
        python workflow/scripts/predict_boundaries.py \
            --precursors {input.faa} --candidates {input.cand} \
            --signalp {input.signalp} --truth {input.truth} \
            --min-mature {params.min_len} --max-mature {params.max_len} \
            --tol {params.tol} \
            --out {output.regions} --accuracy {output.accuracy}
        """

rule cutsites_to_genome:
    # Project cleavage sites onto GRCh38 as a 1-bp BED track.
    input:
        cand = OUT + "/cleavage/candidate_sites.tsv",
        gff  = OUT + "/genome_map/precursors.miniprot.gff",
    output:
        bed = OUT + "/cleavage/cutsites.bed",
    conda:
        "../envs/core.yaml"
    shell:
        r"""
        python workflow/scripts/cutsites_to_genome.py \
            --candidates {input.cand} --gff {input.gff} --out {output.bed}
        """
