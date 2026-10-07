# --- Frog AMP -> genome training-set rules ---------------------------------
# Per-species fan-out over genomed frog species (SPECIES from Snakefile_frog).
# Reuses the core scripts (best_hit_precursor, define_regions, scan_cleavage,
# train/score cleavage, locus_features) against each species' own genome+proteome.

P = config["params"]

# --- AMP acquisition + species routing -------------------------------------
rule frog_fetch_apd:
    output:
        faa = OUT + "/amp/peptides.faa",
        tsv = OUT + "/amp/peptides.tsv",
    params:
        off_fa = config["apd"]["offline_fasta"],
        off_meta = config["apd"].get("offline_metadata", ""),
    shell:
        r"""
        mkdir -p {OUT}/amp
        python workflow/scripts/fetch_apd_frog.py \
            --offline-fasta {params.off_fa} --offline-meta "{params.off_meta}" \
            --out-faa {output.faa} --out-tsv {output.tsv}
        """

rule frog_assign_species:
    input:
        tsv = OUT + "/amp/peptides.tsv",
        sp  = config["species_table"],
        fam = config["family_table"],
    output:
        asn = OUT + "/amp/assignments.tsv",
        faas = expand(OUT + "/amp/by_species/{sp}.faa", sp=SPECIES.keys()),
    params:
        species_keys = " ".join(SPECIES.keys()),
    shell:
        r"""
        python workflow/scripts/assign_species.py \
            --peptides {input.tsv} --species-table {input.sp} \
            --family-table {input.fam} --outdir {OUT}/amp
        mkdir -p {OUT}/amp/by_species
        # ensure a FASTA exists for every genomed species (empty if none mapped)
        for s in {params.species_keys}; do touch {OUT}/amp/by_species/$s.faa; done
        """

# --- Reference acquisition per species (NCBI datasets) ---------------------
rule frog_download_genome:
    output:
        genome = RES + "/{sp}/genome.fa",
        prot   = RES + "/{sp}/proteome.faa",
        gff    = RES + "/{sp}/annotation.gff",
    params:
        acc = lambda wc: SPECIES[wc.sp]["acc"],
    conda: "../envs/frog.yaml"
    shell:
        r"""
        mkdir -p {RES}/{wildcards.sp}
        cd {RES}/{wildcards.sp}
        datasets download genome accession {params.acc} \
            --include genome,protein,gff3 --filename ds.zip
        unzip -o ds.zip -d ds >/dev/null
        d=ds/ncbi_dataset/data/{params.acc}
        cat $d/*_genomic.fna > genome.fa
        cat $d/protein.faa   > proteome.faa
        cat $d/genomic.gff   > annotation.gff
        samtools faidx genome.fa
        """

# --- Peptide -> precursor -> genome (per species) --------------------------
rule frog_precursor:
    input:
        pep  = OUT + "/amp/by_species/{sp}.faa",
        prot = RES + "/{sp}/proteome.faa",
    output:
        mapp = OUT + "/map/{sp}/peptide_to_precursor.tsv",
        faa  = OUT + "/map/{sp}/precursors.faa",
    conda: "../envs/frog.yaml"
    threads: P["threads"]
    shell:
        r"""
        mkdir -p {OUT}/map/{wildcards.sp}
        diamond makedb --in {input.prot} -d {OUT}/map/{wildcards.sp}/prot -p {threads}
        diamond blastp -q {input.pep} -d {OUT}/map/{wildcards.sp}/prot -p {threads} \
            --evalue {P[diamond_evalue]} --max-target-seqs 5 --very-sensitive \
            --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore stitle full_sseq \
            -o {OUT}/map/{wildcards.sp}/blast.tsv || true
        python workflow/scripts/best_hit_precursor.py \
            --blast {OUT}/map/{wildcards.sp}/blast.tsv \
            --min-pident {P[diamond_min_pident]} --min-cov {P[min_query_cov]} \
            --out-map {output.mapp} --out-fasta {output.faa}
        """

rule frog_genome_map:
    input:
        faa    = OUT + "/map/{sp}/precursors.faa",
        genome = RES + "/{sp}/genome.fa",
    output:
        gff = OUT + "/map/{sp}/precursors.miniprot.gff",
    conda: "../envs/frog.yaml"
    threads: P["threads"]
    shell:
        r"""
        miniprot -t {threads} --gff {input.genome} {input.faa} > {output.gff}
        """

rule frog_signalp:
    input: faa = OUT + "/map/{sp}/precursors.faa",
    output: tsv = OUT + "/map/{sp}/signalp.tsv",
    conda: "../envs/frog.yaml"
    threads: P["threads"]
    shell: "deepsig -f {input.faa} -o {output.tsv} -k euk -t {threads} || touch {output.tsv}"

rule frog_regions:
    input:
        prec = OUT + "/map/{sp}/precursors.faa",
        pep  = OUT + "/amp/by_species/{sp}.faa",
        mapp = OUT + "/map/{sp}/peptide_to_precursor.tsv",
        sig  = OUT + "/map/{sp}/signalp.tsv",
    output: reg = OUT + "/map/{sp}/regions_aa.tsv",
    conda: "../envs/frog.yaml"
    shell:
        r"""
        python workflow/scripts/define_regions.py \
            --precursors {input.prec} --peptides {input.pep} --map {input.mapp} \
            --signalp {input.sig} --dibasic KR,RR,RK,KK --out {output.reg}
        """

rule frog_cleavage:
    input:
        prec = OUT + "/map/{sp}/precursors.faa",
        reg  = OUT + "/map/{sp}/regions_aa.tsv",
        mot  = config.get("protease_motifs", "config/protease_motifs.tsv"),
    output: scored = OUT + "/map/{sp}/candidate_sites_scored.tsv",
    conda: "../envs/frog.yaml"
    shell:
        r"""
        d={OUT}/map/{wildcards.sp}
        python workflow/scripts/scan_cleavage.py --precursors {input.prec} \
            --motifs {input.mot} --out $d/cand.tsv
        python workflow/scripts/train_cleavage_model.py --precursors {input.prec} \
            --regions {input.reg} --outdir $d/model
        python workflow/scripts/score_sites.py --candidates $d/cand.tsv \
            --precursors {input.prec} --pwm $d/model/pwm.tsv \
            --model $d/model/model.pkl --out {output.scored}
        """

# --- Positive features + locus avoid-BED (per species) ---------------------
rule frog_locus_features:
    input:
        pep  = OUT + "/amp/peptides.tsv",
        reg  = OUT + "/map/{sp}/regions_aa.tsv",
        gff  = OUT + "/map/{sp}/precursors.miniprot.gff",
        scr  = OUT + "/map/{sp}/candidate_sites_scored.tsv",
        gen  = RES + "/{sp}/genome.fa",
    output:
        feat = OUT + "/features/{sp}.positives.tsv",
        bed  = OUT + "/features/{sp}.loci.bed",
    conda: "../envs/frog.yaml"
    params:
        species = lambda wc: SPECIES[wc.sp]["species"],
    shell:
        r"""
        mkdir -p {OUT}/features
        python workflow/scripts/locus_features.py \
            --peptides {input.pep} --regions {input.reg} --gff {input.gff} \
            --cleavage {input.scr} --genome {input.gen} \
            --species "{params.species}" --out {output.feat}
        # avoid-BED: miniprot mRNA spans, so negatives do not overlap AMP loci
        awk -F'\t' '$3=="mRNA"{{print $1"\t"$4-1"\t"$5}}' {input.gff} > {output.bed}
        """

rule frog_negatives:
    input:
        gen = RES + "/{sp}/genome.fa",
        bed = OUT + "/features/{sp}.loci.bed",
    output: neg = OUT + "/features/{sp}.negatives.tsv",
    conda: "../envs/frog.yaml"
    params:
        species = lambda wc: SPECIES[wc.sp]["species"],
        npos = config["model"]["negatives_per_positive"],
        seed = config["model"]["seed"],
    shell:
        r"""
        # negatives ~ (positives * npos); count positives from the loci bed
        n=$(( ($(wc -l < {OUT}/features/{wildcards.sp}.loci.bed) + 1) * {params.npos} ))
        python workflow/scripts/sample_negatives.py --genome {input.gen} \
            --avoid-bed {input.bed} --n $n --species "{params.species}" \
            --seed {params.seed} --out {output.neg}
        """

# --- Aggregate + matrix + model --------------------------------------------
rule frog_matrix:
    input:
        pos = expand(OUT + "/features/{sp}.positives.tsv", sp=SPECIES.keys()),
        neg = expand(OUT + "/features/{sp}.negatives.tsv", sp=SPECIES.keys()),
    output: mat = OUT + "/model/feature_matrix.tsv",
    conda: "../envs/frog.yaml"
    shell:
        r"""
        mkdir -p {OUT}/model
        # concatenate per-species tables (keep one header)
        head -1 {input.pos[0]} > {OUT}/model/positives.tsv
        for f in {input.pos}; do tail -n +2 $f >> {OUT}/model/positives.tsv; done
        head -1 {input.neg[0]} > {OUT}/model/negatives.tsv
        for f in {input.neg}; do tail -n +2 $f >> {OUT}/model/negatives.tsv; done
        python workflow/scripts/build_training_matrix.py \
            --positives {OUT}/model/positives.tsv \
            --negatives {OUT}/model/negatives.tsv --out {output.mat}
        """

rule frog_train:
    input: mat = OUT + "/model/feature_matrix.tsv",
    output: metrics = OUT + "/model/metrics.tsv",
    conda: "../envs/frog.yaml"
    params:
        clf = config["model"]["classifier"],
        cv  = config["model"]["cv_folds"],
        seed = config["model"]["seed"],
    shell:
        r"""
        python workflow/scripts/train_amp_genome_model.py \
            --matrix {input.mat} --outdir {OUT}/model \
            --classifier {params.clf} --cv {params.cv} --seed {params.seed}
        """
