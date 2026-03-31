configfile: "config.yaml"


rule all:
    input:
        "results/expression_clean.tsv",
        "results/correlation_pancancer.png",
        "results/correlation_aml.png",
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
        pancancer="results/correlation_pancancer.png",
        aml="results/correlation_aml.png",
        top_csv="results/{gene1}_{gene2}_top.csv".format(
            gene1=config["gene1"], gene2=config["gene2"]
        )
    log:
        "logs/correlation.log"
    shell:
        "python scripts/02_correlation.py 2> {log}"


rule survival:
    input:
        expression="results/expression_clean.tsv",
        clinical=config["clinical_file"]
    output:
        pancancer="results/survival_pancancer.png",
        aml="results/survival_aml.png"
    log:
        "logs/survival.log"
    shell:
        "python scripts/03_survival.py 2> {log}"