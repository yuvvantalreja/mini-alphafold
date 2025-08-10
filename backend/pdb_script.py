#!/usr/bin/env python3
"""
Condensed ColabFold script for local protein structure prediction
Usage: python colabfold_minimal.py
"""

import os
import re
import hashlib
from pathlib import Path
from colabfold.download import download_alphafold_params
from colabfold.utils import setup_logging
from colabfold.batch import get_queries, run, set_model_type
import string

os.environ["JAX_PLATFORMS"] = "cpu"



def GeneratePDB(query_sequence):


    for i in query_sequence:
        if i not in string.ascii_letters:
            return -1

    def add_hash(x, y):
        return x + "_" + hashlib.sha1(y.encode()).hexdigest()[:5]

    # Input protein sequence
    #query_sequence = 'PIAQIHILEGRSDEQKETLIREVSEAISRSLDAPLTSVRVIITEMAKGHFGIGGELASK'
    jobname = 'test'

    # Clean inputs
    query_sequence = "".join(query_sequence.split())
    basejobname = "".join(jobname.split())
    basejobname = re.sub(r'\W+', '', basejobname)
    jobname = add_hash(basejobname, query_sequence)

    # Check if directory exists and create unique name if needed
    def check(folder):
        return not os.path.exists(folder)

    if not check(jobname):
        n = 0
        while not check(f"{jobname}_{n}"):
            n += 1
        jobname = f"{jobname}_{n}"

    # Create output directory
    os.makedirs(jobname, exist_ok=True)

    # Save query sequence
    queries_path = os.path.join(jobname, f"{jobname}.csv")
    with open(queries_path, "w") as text_file:
        text_file.write(f"id,sequence\n{jobname},{query_sequence}")

    print(f"jobname: {jobname}")
    print(f"sequence: {query_sequence}")
    print(f"length: {len(query_sequence.replace(':', ''))}")

    # Setup logging
    log_filename = os.path.join(jobname, "log.txt")
    setup_logging(Path(log_filename))

    # Get queries and determine if complex
    queries, is_complex = get_queries(queries_path)
    model_type = set_model_type(is_complex, "alphafold2_ptm")

    # Download AlphaFold parameters
    download_alphafold_params(model_type, Path("."))


    
    # Run prediction
    results = run(
        queries=queries,
        result_dir=jobname,
        use_templates=False,
        custom_template_path=None,
        num_relax=0,
        msa_mode="mmseqs2_uniref_env",
        model_type=model_type,
        num_models=1,
        num_recycles=1,
        relax_max_iterations=200,
        recycle_early_stop_tolerance=0.0,
        num_seeds=1,
        use_dropout=False,
        model_order=[1],
        is_complex=is_complex,
        data_dir=Path("."),
        keep_existing_results=False,
        rank_by="auto",
        pair_mode="unpaired_paired",
        pairing_strategy="greedy",
        stop_at_score=float(100),
        prediction_callback=None,
        dpi=200,
        zip_results=False,
        save_all=False,
        max_msa="64:128",
        use_cluster_profile=True,
        input_features_callback=None,
        save_recycles=False,
        user_agent="colabfold/local",
        calc_extra_ptm=False,
    )
    
    return jobname
    
# print(f"Structure prediction completed. PDB files saved in {jobname}/ directory.")
