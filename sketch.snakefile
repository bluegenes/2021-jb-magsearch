"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: conda activate test-magsearch; snakemake -s sketch.snakefile -n
"""

import os
import re
import pandas as pd

configfile: "prep_magsearch.yml"
out_dir = config['outdir']
logs_dir = os.path.join(out_dir, "logs")
basename = config["basename"]
scaled = config["scaled"]
ksize = config["ksize"]
if not isinstance(scaled, list):
    config["scaled"] = [scaled]
if not isinstance(ksize, list):
    config["ksize"] = [ksize]

genome_info  = pd.read_csv(config["genome_csv"], header=0)
genome_info.set_index("name", inplace=True)

rule all:
    input: os.path.join(out_dir, f"{basename}.siglist.txt")

def make_param_str(ksizes, scaled):
    ks = [ f'k={k}' for k in ksizes ]
    ks = ",".join(ks)
    scaled = min(scaled) #take minimum value of scaled list
    return f"{ks},scaled={scaled},abund"

rule sourmash_sketch:
    input:
        lambda w: genome_info.at[w.accession, "filename"],
    output:
        os.path.join(out_dir, "signatures", "{accession}.sig"),
    params:
        sketch_params=make_param_str(config["ksize"], config["scaled"]),
        signame = lambda w: w.accession #genome_info.at[w.accession, "signame"]
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *1000,
        runtime=1200,
    log: os.path.join(logs_dir, "sourmash_sketch", "{accession}.sketch.log")
    benchmark: os.path.join(logs_dir, "sourmash_sketch", "{accession}.sketch.benchmark")
    conda: "envs/sourmash4.yml"
    group: "sketch"
    shell:
        """
        sourmash sketch dna -p {params.sketch_params} -o {output} --name {params.signame:q} {input} 2> {log}
        """

localrules: signames_to_file
rule signames_to_file:
    input: 
        sigs=expand(os.path.join(out_dir, "signatures", "{accession}.sig"), accession = genome_info.index),
    output: os.path.join(out_dir, f"{basename}.siglist.txt")
    run:
        with open(str(output), "w") as outF:
            for inF in input.sigs:
                full_filename = os.path.abspath(str(inF))
                #outF.write(str(inF) + "\n")
                outF.write(full_filename + "\n")

