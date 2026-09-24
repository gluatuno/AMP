# --- Step 0: prep & de-duplicate mature peptides ---------------------------
# Many APD entries are fragments/variants of the same parent (e.g. the LL-37
# series). We keep the full set but also produce a non-redundant cluster set,
# and a cleaned FASTA (single-line, no blank records).

rule clean_input:
    input:
        fa = config["input_peptides"],
    output:
        fa = OUT + "/prep/peptides.clean.fa",
    conda:
        "../envs/core.yaml"
    shell:
        r"""
        mkdir -p {OUT}/prep
        seqkit seq -w 0 {input.fa} | seqkit rmdup -s -o {output.fa} 2> {OUT}/prep/rmdup.log
        """

rule cluster_input:
    input:
        fa = OUT + "/prep/peptides.clean.fa",
    output:
        rep = OUT + "/prep/peptides.nr.fa",
        clstr = OUT + "/prep/peptides.nr.fa.clstr",
    conda:
        "../envs/core.yaml"
    threads: config["params"]["threads"]
    shell:
        # short peptides -> permissive word size; 90% identity clustering
        "cd-hit -i {input.fa} -o {output.rep} -c 0.9 -n 2 -T {threads} -d 0 -M 0"
