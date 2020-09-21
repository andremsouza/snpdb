# %% [markdown]
#  # Experimentos SNPDB

# %% [markdown]
#  ## Imports e inicialização

# %%
from experiment_config import data_dir, results_fname, exps, nsnps_ids
from experiment_config import nsamples_ids, reset_db, generate_random_file
import json
import numpy as np
import os
import snpdb
import time

# Verifying existing FastQ files
print('Directory:', data_dir, '\n', os.listdir(data_dir), '\n')

# %% [markdown]
# # 2 - Análise de importação e consultas para as quantidades de indivíduos:
# - 10
# - 100
# - 1000
# - 10000
# - 100000
# - 1000000

# %%
# Import existing results
results: dict = {}
try:
    with open(results_fname, 'r') as f:
        results = json.load(f)
except Exception:
    results = {}

# %% [markdown]
# ## Definição de funções

# %%


def execute_experiment_two_files_bin(result: dict,
                                     experiment_id: str,
                                     compression_method: str,
                                     nsnps: int,
                                     nsamples: int,
                                     N: int = 1) -> None:
    """Executes experiment for file types with two files (map, sample)
    Warning: this function uses variables defined outside its scope

    Args:
        result (dict): Dictionary for experiment. Values will be assigned
                       to this dictionary "by reference".
    """
    nsnps_id = nsnps_ids[nsnps]
    nsamples_id = nsamples_ids[nsamples]
    print("Starting Experiment " + experiment_id + " (" + nsnps_id +
          " SNPs, " + nsamples_id + " individuals) with N = " + str(N) +
          "; Compression method: " + compression_method)

    # get filenames
    f_ext: dict = exps[experiment_id]['file_extensions']
    mfname = str(data_dir + 'out_' + nsnps_id + '_' + nsamples_id +
                 f_ext['map'])
    pfname = str(data_dir + 'out_' + nsnps_id + '_' + nsamples_id +
                 f_ext['ped'])
    ifname = str(data_dir + 'out_' + nsnps_id + '_' + nsamples_id +
                 f_ext['ids'])

    # Setting up result dictionary
    result['fsize'] = []  # file size
    result['dbsize'] = []  # doc size in db
    result['time'] = []  # insertion time
    result['summarize'] = []  # summarization example and time
    result['individuals_of_snps'] = []
    result['delete_individual'] = []

    # * Performing experiment N times and storing results
    for i in range(N):
        print("i: " + str(i))
        print("Resetting database...")
        reset_db(compression_method=compression_method)
        print("Database reset operation successful.")
        print("Generating input files...")
        t_map: float = 0.0
        t_sample: float = 0.0
        file_size: float = 0.0
        # * Generating input files
        # If less than 10000 samples, generate in one file
        # Else, generate blocks of up to 10000 samples
        n_blocks: int = int(np.ceil(nsamples / 10000))
        remaining_samples: int = nsamples
        start_sample: int = 1
        # Map file
        generate_random_file(filename=mfname,
                             file_type=f_ext['map'],
                             verbose=True,
                             n=nsnps)
        #  start_from_id=start_map)
        # Importing map file
        with open(mfname, 'rb') as mf:
            t_tmp: float = time.time()
            snpdb.insert_file(mf,
                              map_name=experiment_id + '_' + nsnps_id + '_' +
                              nsamples_id,
                              file_type=f_ext['map'])
            t_tmp = time.time() - t_tmp
        t_map += t_tmp
        for i in range(n_blocks):
            print("Block: " + str(i))
            nsamples_block = int(np.minimum(remaining_samples, 10000.0))
            # Samples file
            generate_random_file(filename=pfname,
                                 file_type=f_ext['ped'],
                                 verbose=True,
                                 n=nsamples_block,
                                 map_size=nsnps,
                                 start_from_id=start_sample)
            # Id map file
            generate_random_file(filename=ifname,
                                 file_type='.ids',
                                 verbose=True,
                                 n=nsamples_block,
                                 first_sample_id=start_sample)
            start_sample += nsamples_block
            remaining_samples -= nsamples_block

            # Importing sample file
            id_map: dict = {}
            # Linking samples to individuals in the database
            if ifname is not None:
                with open(ifname, "r") as f:
                    for line in f:
                        (sample, individual) = line.split()
                        id_map[sample] = individual
            with open(pfname, 'rb') as pf:
                t_tmp = time.time()
                snpdb.insert_file(pf,
                                  map_name=experiment_id + '_' + nsnps_id +
                                  '_' + nsamples_id,
                                  file_type=f_ext['ped'])
                t_tmp = time.time() - t_tmp
            t_sample += t_tmp
            with open(ifname, 'rb') as iff:
                t_tmp = time.time()
                snpdb.insert_file(iff,
                                  map_name=experiment_id + '_' + nsnps_id +
                                  '_' + nsamples_id,
                                  file_type='.ids')
                t_tmp = time.time() - t_tmp
            t_sample += t_tmp
        # Validating Statistics
        snpdb._db.command("validate", "fs.chunks", full=True)
        snpdb._db.command("validate", "fs.files", full=True)
        file_size = snpdb._db.command("collstats", "fs.chunks")["storageSize"]
        # file_size = snpdb._db["fs.files"].find_one()["length"]
        t: float = t_map + t_sample
        print("Imported file file\tTime: " + str(round(t, 3)) + "s\tSize:" +
              str(round(file_size / 1024**2, 2)) + "MB")

        # Appending generated file sizes
        result['fsize'].append(
            float(os.stat(mfname).st_size) + float(os.stat(pfname).st_size) +
            float(os.stat(ifname).st_size))
        # Appending stored document sizes from MongoDB
        result['dbsize'].append(file_size)
        # Appending insertion times
        result['time'].append(t_map + t_sample)

        # Writing partial results to file
        with open(results_fname, 'w') as f:
            json.dump(results,
                      f,
                      ensure_ascii=True,
                      check_circular=True,
                      allow_nan=True,
                      indent=1,
                      sort_keys=True)


# %% [markdown]
# ## 2.7 - Comparação de importação e exportação de dados no modo “bruto”
# Deseja-se importar e exportar dados no modo “bruto”, e comparar a
# performance dessas operações com as usuais, utilizando funções
# pré-existentes. Inicialmente, deseja-se utilizar para este experimento um
# arquivo 0125 de 1GB.

# %%
# ? Using 1GB 0125 file for experiment (100kSNPs/10k individuals)
experiment_ids = ['2.7']
for experiment_id in experiment_ids:
    if experiment_id.startswith('1'):
        N = 10  # Performing experiments with N loops
    else:
        N = 1  # Only 1 loop for experiments 2.*
    file_type = exps[experiment_id]['file_type']
    results.setdefault(experiment_id, {})
    for compression_method in exps[experiment_id]['compression_methods']:
        results[experiment_id].setdefault(compression_method, {})
        for nsnps in exps[experiment_id]['nsnps_list']:
            nsnps_id = nsnps_ids[nsnps]
            nsnps = int(nsnps)
            results[experiment_id][compression_method].setdefault(nsnps_id, {})
            for nsamples in exps[experiment_id]['nsamples_list']:
                nsamples_id = nsamples_ids[nsamples]
                nsamples = int(nsamples)
                results[experiment_id][compression_method][
                    nsnps_id].setdefault(nsamples_id, {})
                # ? Using pointers to results dictionary
                result = results[experiment_id][compression_method][nsnps_id][
                    nsamples_id]
                # Only execute experiment if empty
                if len(result) == 0:
                    execute_experiment_two_files_bin(result, experiment_id,
                                                     compression_method, nsnps,
                                                     nsamples, N)

# %% [markdown]
# # Saving results

# %%
# Writing results to file
with open(results_fname, 'w') as f:
    json.dump(results,
              f,
              ensure_ascii=True,
              check_circular=True,
              allow_nan=True,
              indent=1,
              sort_keys=True)
