configfile: "config.yaml"

if "gene1" not in config or "gene2" not in config:
    raise ValueError(
        "Specify two genes: snakemake --cores 1 --config genes
    )

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
        pancancer_gh="results/correlation_pancancer_github.png",
        aml="results/correlation_aml.png",
        aml_gh="results/correlation_aml_github.png",
        top_csv="results/{gene1}_{gene2}_top.csv".format(
            gene1=config["gene1"], gene2=config["gene2"]
        )
    log:
        "logs/correlation.log"
    shell:
        "python scripts/02_correlation.py "
        "--gene1 {config[gene1]} "
        "--gene2 {config[gene2]} "
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
        "--gene1 {config[gene1]} "
        "--gene2 {config[gene2]} "
        "2> {log}"
