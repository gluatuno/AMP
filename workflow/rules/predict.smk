# --- Context-aware AMP predictor rules -------------------------------------

D = config["dataset"]
M = config["model"]
PR = config["predict"]
EMB = M.get("embeddings", {}) or {}
USE_EMB = bool(EMB.get("enabled", False))

rule predict_make_dataset:
    input:
        records = D["positive_records"],
    output:
        tsv = OUT + "/dataset.tsv",
    conda: "../envs/predict.yaml"
    params:
        neg = D.get("neg_fasta", ""),
        shuf = D.get("shuffle_per_pos", 0),
        seed = M["seed"],
    shell:
        r"""
        mkdir -p {OUT}
        python workflow/scripts/make_dataset.py \
            --records {input.records} --neg-fasta "{params.neg}" \
            --shuffle-per-pos {params.shuf} --seed {params.seed} \
            --out {output.tsv}
        """

# --- Method 6: foundation-model embeddings -> feature table ----------------
# Produces dataset_features.tsv (dataset + m6p_*/m6g_* columns when enabled,
# else a straight copy). train_ablation picks up M6 automatically.
if USE_EMB:
    rule predict_m6_features:
        input:
            tsv = OUT + "/dataset.tsv",
        output:
            feat = OUT + "/dataset_features.tsv",
        conda: "../envs/predict.yaml"
        params:
            backend = EMB.get("backend", "stub"),
            pdim = EMB.get("protein_dim", 32),
            ddim = EMB.get("dna_dim", 32),
            locus = EMB.get("locus_fasta", ""),
        shell:
            r"""
            d={OUT}/m6; mkdir -p $d
            # protein embedding (M6p) over every row's mature sequence
            awk -F'\t' 'NR>1{{print ">"$1"\n"$3}}' {input.tsv} > $d/seqs.faa
            python workflow/scripts/embed_sequences.py --fasta $d/seqs.faa \
                --kind protein --backend {params.backend} --dim {params.pdim} \
                --out $d/m6p.tsv
            EMB_ARGS="$d/m6p.tsv"
            # optional DNA embedding (M6g) over a provided locus FASTA
            if [ -n "{params.locus}" ] && [ -f "{params.locus}" ]; then
                python workflow/scripts/embed_sequences.py --fasta {params.locus} \
                    --kind dna --backend {params.backend} --dim {params.ddim} \
                    --out $d/m6g.tsv
                EMB_ARGS="$EMB_ARGS $d/m6g.tsv"
            fi
            python workflow/scripts/merge_embeddings.py --dataset {input.tsv} \
                --emb $EMB_ARGS --out {output.feat}
            """
else:
    rule predict_passthrough_features:
        input:
            tsv = OUT + "/dataset.tsv",
        output:
            feat = OUT + "/dataset_features.tsv",
        shell:
            "cp {input.tsv} {output.feat}"

rule predict_train_ablation:
    input:
        tsv = OUT + "/dataset_features.tsv",
    output:
        abl = OUT + "/model/ablation.tsv",
        met = OUT + "/model/metrics.tsv",
        m1  = OUT + "/model/model_M1.pkl",
        full = OUT + "/model/model_full.pkl",
    conda: "../envs/predict.yaml"
    params:
        cv = M["cv_folds"],
        seed = M["seed"],
    shell:
        r"""
        python workflow/scripts/train_ablation.py \
            --dataset {input.tsv} --outdir {OUT}/model \
            --cv {params.cv} --seed {params.seed}
        """

rule predict_score:
    input:
        model = OUT + "/model/model_M1.pkl",
        fasta = PR["query_fasta"],
    output:
        tsv = OUT + "/predictions.tsv",
    conda: "../envs/predict.yaml"
    params:
        thr = PR.get("threshold", 0.5),
    shell:
        r"""
        python workflow/scripts/predict_amp.py \
            --model {input.model} --fasta {input.fasta} \
            --threshold {params.thr} --out {output.tsv}
        """
