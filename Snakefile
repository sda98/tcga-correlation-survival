cat > Snakefile << 'EOF'
configfile: "config.yaml"

def require_genes():
    if "genes" not in config:
        raise ValueError(
            "Specify genes: snakemake --cores 5 --latency-wait 30 --config genes=GENE1,GENE2"
        )
    return config["genes"]

rule all:
    input:
        "results/expression_clean.tsv",
        "results/correlation_done.txt",
        "results/survival_done.txt"

rule data_prep:
    input:
        clinical=config["clinical_file"]
    output:
        expression="results/expression_clean.tsv",
        survival="results/survival_clean.tsv"
    log:
        "logs/data_prep.log"
    shell:
        "python scripts/01_data_prep.py 2> {log}"

rule correlation:
    input:
        expression="results/expression_clean.tsv"
    output:
        "results/correlation_done.txt"
    params:
        genes=lambda wildcards: require_genes()
    log:
        "logs/correlation.log"
    shell:
        "python scripts/02_correlation.py --genes {params.genes} 2> {log}"

rule survival:
    input:
        expression="results/expression_clean.tsv",
        clinical=config["clinical_file"]
    output:
        "results/survival_done.txt"
    params:
        genes=lambda wildcards: require_genes()
    log:
        "logs/survival.log"
    shell:
        "python scripts/03_survival.py --genes {params.genes} 2> {log}"
EOF
