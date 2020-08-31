# %% [markdown]
#  # Experimentos SNPDB

# %% [markdown]
#  ## Imports e inicialização

# %%
from experiment2_1 import reset_db
import json
import os
import readers
import snpdb
import time

data_dir: str = './data/'  # Data output directory
fastq_dir_1: str = '../1_FASTq_VCF/1_SRR10752699_MISSOURI/'  # data directory
fastq_dir_2: str = '../1_FASTq_VCF/2_VCF_bos_taurus/'  # data directory
graph_dir: str = './graphs/'  # Graph output directory

# Verifying existing FastQ files
print('Directory:', data_dir, '\n', os.listdir(data_dir), '\n')
print('Directory:', fastq_dir_1, '\n', os.listdir(fastq_dir_1), '\n')
print('Directory:', fastq_dir_2, '\n', os.listdir(fastq_dir_2), '\n')

# %% [markdown]
# # 2 - Análise de importação e consultas para as quantidades de indivíduos:
# - 10
# - 100
# - 1000
# - 10000
# - 100000
# - 1000000

# %%
N: int = 1  # Performing experiments with N loops
results: dict = {}  # Storing results in a dictionary
nspns_list: list = [100000, 1000000]
nsamples_list: list = [10, 100, 1000, 10000, 100000, 1000000]

# %% [markdown]
# ## 2.7 - Comparação de importação e exportação de dados no modo “bruto”
# Deseja-se importar e exportar dados no modo “bruto”, e comparar a
# performance dessas operações com as usuais, utilizando funções
# pré-existentes. Inicialmente, deseja-se utilizar para este experimento um
# arquivo 0125 de 1GB.

# %%
# ? Using 1GB 0125 file for experiment (100kSNPs/10k individuals)
experiment_id = '2.7'
results[experiment_id] = {}

for nsnps in [100000]:  # ! Verify snp count
    size_id = {100000: '100k', 1000000: '1m'}[nsnps]
    for nsamples in [10000]:  # ! Verify sample count
        count_id = {
            10: '10',
            100: '100',
            1000: '1k',
            10000: '10k',
            100000: '100k',
            1000000: '1m'
        }[nsamples]
        # Printing experiment start information
        print("Starting Experiment " + experiment_id + " (" + size_id +
              " SNPs, " + count_id + " individuals) with N = " + str(N))

        # start_map: int = 1  # starting snp id
        # start_sample: int = 1  # starting sample id
        t: float = 0.0  # timing variable

        # Filenames
        mfname = data_dir + 'out_' + size_id + '_' + count_id + '.0125map'
        pfname = data_dir + 'out_' + size_id + '_' + count_id + '.0125ped'
        ifname = data_dir + 'out_' + size_id + '_' + count_id + '.0125ids'

        # Storing results in a dictionary (saving to a JSON file)
        results[experiment_id][size_id] = {}
        results[experiment_id][size_id][count_id] = {}
        results[experiment_id][size_id][count_id]['fsize'] = []
        results[experiment_id][size_id][count_id]['dbsize'] = []
        # map result data
        results[experiment_id][size_id][count_id]['map'] = {}
        results[experiment_id][size_id][count_id]['map']['fsize'] = []
        results[experiment_id][size_id][count_id]['map']['dbsize'] = []
        results[experiment_id][size_id][count_id]['map']['time'] = []
        # sample result data
        results[experiment_id][size_id][count_id]['sample'] = {}
        results[experiment_id][size_id][count_id]['sample']['fsize'] = []
        results[experiment_id][size_id][count_id]['sample']['dbsize'] = []
        results[experiment_id][size_id][count_id]['sample']['time'] = []
        # binary file  result data
        results[experiment_id][size_id][count_id]['binary'] = {}
        results[experiment_id][size_id][count_id]['binary']['fsize'] = []
        results[experiment_id][size_id][count_id]['binary']['dbsize'] = []
        results[experiment_id][size_id][count_id]['binary']['time'] = []

        # * Performing experiment N times and storing results
        for i in range(N):
            print("i: " + str(i))
            print("Resetting database...")
            reset_db()
            print("Database reset operation successful.")
            # * Inserting input files into db
            print("Inserting input files into database...")
            # Importing map file
            t = time.time()
            snpdb.import_map(map_reader=readers.Z125MapReader(mfname),
                             map_name=experiment_id + size_id,
                             force_create_new=True,
                             force_use_existing=False,
                             report=False)
            t = time.time() - t
            results[experiment_id][size_id][count_id]['map']['time'].append(t)
            # Validating statistics
            snpdb._db.command("validate",
                              snpdb._config["MAPS_COLL"],
                              full=True)
            snpdb._db.command("validate",
                              snpdb._config["MAPSNPS_COLL"],
                              full=True)
            snpdb._db.command("validate",
                              snpdb._config["SNPS_COLL"],
                              full=True)
            # map_size = _MAPS + _MAPSNPS + _SNPS
            map_size = snpdb._db.command(
                "collstats",
                snpdb._config["MAPS_COLL"])['storageSize'] + snpdb._db.command(
                    "collstats", snpdb._config["MAPSNPS_COLL"]
                )['storageSize'] + snpdb._db.command(
                    "collstats", snpdb._config["SNPS_COLL"])['storageSize']
            print("Imported map file\tTime: " + str(round(t, 3)) + "s\tSize:" +
                  str(round(map_size / 1024**2, 2)) + "MB")
            # Importing sample file
            id_map = {}
            # Linking samples to individuals in the database
            if ifname is not None:
                with open(ifname, "r") as f:
                    for line in f:
                        (sample, individual) = line.split()
                        id_map[sample] = individual
            t = time.time()
            snpdb.import_samples(
                sample_reader=readers.Z125SampleReader(pfname),
                map_name=experiment_id + size_id,
                id_map=id_map,
                report=False)
            t = time.time() - t
            results[experiment_id][size_id][count_id]['sample']['time'].append(
                t)
            # Validating Statistics
            snpdb._db.command("validate",
                              snpdb._config["SAMPLES_COLL"],
                              full=True)
            snpdb._db.command("validate",
                              snpdb._config["SNPBLOCKS_COLL"],
                              full=True)
            snpdb._db.command("validate",
                              snpdb._config["INDIVIDUALS_COLL"],
                              full=True)
            # sample_size = _SAMPLES + _SNPBLOCKS + _INDS
            sample_size = snpdb._db.command(
                "collstats", snpdb._config["SAMPLES_COLL"]
            )['storageSize'] + snpdb._db.command(
                "collstats", snpdb._config["SNPBLOCKS_COLL"]
            )['storageSize'] + snpdb._db.command(
                "collstats", snpdb._config["INDIVIDUALS_COLL"])['storageSize']
            print("Imported samples file\tTime: " + str(round(t, 3)) +
                  "s\tSize:" + str(round(sample_size / 1024**2, 2)) + "MB")

            # Resetting database and inserting as binary files
            print("Resetting database...")
            reset_db()
            print("Database reset operation successful.")
            # * Inserting input files into db
            # Importing files in binary mode
            with open(mfname,
                      "rb") as mf, open(pfname,
                                        "rb") as pf, open(ifname, "rb") as iff:
                t = time.time()
                snpdb.insert_file(file=mf, individual_id="IND1")
                snpdb.insert_file(file=pf, individual_id="IND1")
                snpdb.insert_file(file=iff, individual_id="IND1")
                t = time.time() - t
                # Validating Statistics
                snpdb._db.command("validate", "fs.chunks", full=True)
                snpdb._db.command("validate", "fs.files", full=True)
                file_size = snpdb._db.command("collstats",
                                              "fs.chunks")["storageSize"]
                # file_size = snpdb._db["fs.files"].find_one()["length"]
                results[experiment_id][size_id][count_id]['binary'][
                    'time'].append(t)
                print("Imported 0125 file in binary mode\tTime: " +
                      str(round(t, 3)) + "s\tSize:" +
                      str(round(file_size / 1024**2, 2)) + "MB")

            # Appending generated file sizes
            results[experiment_id][size_id][count_id]['map']['fsize'].append(
                os.stat(mfname).st_size)
            results[experiment_id][size_id][count_id]['sample'][
                'fsize'].append(os.stat(pfname).st_size)
            results[experiment_id][size_id][count_id]['fsize'].append(
                os.stat(mfname).st_size + os.stat(pfname).st_size)
            # Appending stored document sizes from MongoDB
            results[experiment_id][size_id][count_id]['map']['dbsize'].append(
                map_size)
            results[experiment_id][size_id][count_id]['sample'][
                'dbsize'].append(sample_size)
            results[experiment_id][size_id][count_id]['dbsize'].append(
                map_size + sample_size)
            results[experiment_id][size_id][count_id]['binary'][
                'fsize'].append(
                    os.stat(mfname).st_size + os.stat(pfname).st_size +
                    os.stat(ifname).st_size)
            results[experiment_id][size_id][count_id]['binary'][
                'dbsize'].append(file_size)

# %% [markdown]
# # Armazenando resultados

# %%
# Saving results to a JSON file
jfname = data_dir + 'exp2_7results.json'
with open(jfname, 'w') as jf:
    json.dump(results, jf)
