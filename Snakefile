configfile: "config.yaml"

if "genes" not in config:
    raise ValueError(
        "Specify genes: snakemake --cores 1 --config genes=GENE,GENE"
    )

rule all:
    input:
        "results/expression_clean.tsv",
        "results/correlation_done.txt",
        "results/survival_pancancer.png",
        "results/survival_aml.png"

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
    log:
        "logs/correlation.log"
    shell:
        "python scripts/02_correlation.py "
        "--genes {config[genes]} "
        "2> {log}"

rule survival:
    input:
        expression="results/expression_clean.tsv",
        clinical=config["clinical_file"]
    output:
        pancancer="results/survival_pancancer.png",
        pancancer_gh="results/survival_pancancer_github.png",
        aml="results/survival_aml.png",
        aml_gh="results/survival_aml_github.png"
    log:
        "logs/survival.log"
    shell:
        "python scripts/03_survival.py "
        "--genes {config[genes]} "
        "2> {log}"
