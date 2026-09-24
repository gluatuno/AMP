# --- Steps 4-5: extract genomic DNA for each figure "Method" ---------------
# extract_methods.py converts amino-acid region coordinates into genomic BED
# intervals through the miniprot CDS model, then bedtools getfasta pulls DNA.
#   Method 1  AMP (mature)
#   Method 2  Pro-peptide
#   Method 3  Pro-peptide + AMP
#   Method 5  Transporter..AMP locus span (full precursor CDS footprint)
# Method 4 (Enzyme) is a trans-acting locus resolved from the curated table
# against the Ensembl GFF3 and extracted with gffread.

rule methods_bed:
    input:
        regions = OUT + "/annotate/regions_aa.tsv",
        gff     = OUT + "/genome_map/precursors.miniprot.gff",
    output:
        m1 = OUT + "/methods/method1_amp.bed",
        m2 = OUT + "/methods/method2_propeptide.bed",
        m3 = OUT + "/methods/method3_pro_plus_amp.bed",
        m5 = OUT + "/methods/method5_locus.bed",
    conda:
        "../envs/core.yaml"
    shell:
        r"""
        mkdir -p {OUT}/methods
        python workflow/scripts/extract_methods.py \
            --regions {input.regions} \
            --gff {input.gff} \
            --outdir {OUT}/methods
        """

def _bed_for(m):
    return {
        "method1_amp":            OUT + "/methods/method1_amp.bed",
        "method2_propeptide":     OUT + "/methods/method2_propeptide.bed",
        "method3_pro_plus_amp":   OUT + "/methods/method3_pro_plus_amp.bed",
        "method5_locus":          OUT + "/methods/method5_locus.bed",
    }[m]

rule methods_fasta:
    input:
        bed = lambda wc: _bed_for(wc.mname),
        fa  = RES + "/GRCh38.primary_assembly.fa",
        fai = RES + "/GRCh38.primary_assembly.fa.fai",
    output:
        fa = OUT + "/methods/{mname}.fa",
    wildcard_constraints:
        mname = r"method[1235]_[a-z_]+",
    conda:
        "../envs/core.yaml"
    shell:
        # -s = strand-aware, -name = use the BED name column as the FASTA header
        "bedtools getfasta -s -name -fi {input.fa} -bed {input.bed} -fo {output.fa}"

rule method4_enzyme:
    # Trans-acting processing enzymes: resolve gene symbols from the curated
    # table, pull their loci from the Ensembl GFF3, and extract with gffread.
    input:
        map   = OUT + "/homology/peptide_to_precursor.tsv",
        gff   = RES + "/GRCh38.annotation.gff3",
        fa    = RES + "/GRCh38.primary_assembly.fa",
        fai   = RES + "/GRCh38.primary_assembly.fa.fai",
        table = config["enzyme_table"],
    output:
        fa  = OUT + "/methods/method4_enzyme.fa",
        tsv = OUT + "/methods/method4_enzyme.tsv",
    conda:
        "../envs/core.yaml"
    shell:
        r"""
        python workflow/scripts/enzyme_loci.py \
            --map {input.map} \
            --gff {input.gff} \
            --genome {input.fa} \
            --table {input.table} \
            --out-fasta {output.fa} \
            --out-tsv {output.tsv}
        """

# --- Report ----------------------------------------------------------------
rule report:
    input:
        map     = OUT + "/homology/peptide_to_precursor.tsv",
        regions = OUT + "/annotate/regions_aa.tsv",
        m1 = OUT + "/methods/method1_amp.fa",
        m2 = OUT + "/methods/method2_propeptide.fa",
        m3 = OUT + "/methods/method3_pro_plus_amp.fa",
        m4 = OUT + "/methods/method4_enzyme.tsv",
        m5 = OUT + "/methods/method5_locus.fa",
    output:
        summary = OUT + "/report/pipeline_summary.tsv",
        featmap = OUT + "/report/precursor_feature_map.tsv",
    conda:
        "../envs/core.yaml"
    shell:
        r"""
        mkdir -p {OUT}/report
        python workflow/scripts/make_report.py \
            --map {input.map} --regions {input.regions} \
            --methods-dir {OUT}/methods \
            --out-summary {output.summary} \
            --out-featmap {output.featmap}
        """
