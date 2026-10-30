# --- Context-aware AMP predictor rules -------------------------------------

D = config["dataset"]
M = config["model"]
PR = config["predict"]

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

rule predict_train_ablation:
    input:
        tsv = OUT + "/dataset.tsv",
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
