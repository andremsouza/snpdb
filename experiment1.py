# %% [markdown]
#  # Experimentos SNPDB
#  ## Imports e inicialização

# %%
import importlib
import json
import numpy as np
import os
from PIL import Image
import snpdb
import subprocess
import readers
import testing.random_file_generator as rfgen
import time
# import writers

data_dir: str = './data/'  # Data output directory
fastq_dir_1: str = '../1_FASTq_VCF/1_SRR10752699_MISSOURI/'  # data directory
fastq_dir_2: str = '../1_FASTq_VCF/2_VCF_bos_taurus/'  # data directory
graph_dir: str = './graphs/'  # Graph output directory

# Verifying existing FastQ files
print('Directory:', data_dir, '\n', os.listdir(data_dir), '\n')
print('Directory:', fastq_dir_1, '\n', os.listdir(fastq_dir_1), '\n')
print('Directory:', fastq_dir_2, '\n', os.listdir(fastq_dir_2), '\n')

# %% [markdown]
# ## Definição de funções


# %%
def reset_db():
    """Reset MongoDB database for experiments."""
    snpdb._client.drop_database(snpdb._config["DB_NAME"])
    subprocess.run("""mongo --eval "load('./mongo_setup.js')" """,
                   stdout=subprocess.DEVNULL,
                   stderr=subprocess.STDOUT,
                   shell=True)
    importlib.reload(snpdb)


# %% [markdown]
#  ## 1 - Tempo de importação e espaço ocupado a partir da base zerada
#  A partir da base zerada, importar um arquivo de amostra e mapa de um animal
#
#  Considerando cada caso abaixo:
#  - A: 1 arquivo 0125 (incluindo Mapa)
#  - B: 1 arquivo Final Report (FR) (incluindo Mapa)
#  - C: 1 arquivo VCF
#  - D: 1 arquivo PLINK (incluindo Mapa)
#  - E: 1 arquivo FastQ (sequenciamento)
#  - F: arquivo de mídia
#  - G: todos os arquivos de A a F
#

# %%
N: int = 10  # Performing experiments with N loops
results: dict = {}  # Storing results in a dictionary

# %% [markdown]
# ## 1.A - Arquivo 0125

# %%
experiment_id: str = '1.A'
results[experiment_id] = {}

# %% [markdown]
# ### 100k SNPs

# %%
size_id: str = '100k'

# Printing experiment start information
print("Starting Experiment " + experiment_id + " (" + size_id +
      " SNPs) with N = " + str(N))

# Map and sample sizes
start_map: int = 1  # starting snp id
nsnps: int = 100000  # number of snps in map
start_sample: int = 1  # starting sample id
nsamples: int = 1  # number of samples in file # ! Definir número de amostras
t: float = 0  # timing variable

# Filenames
mfname: str = data_dir + 'out_' + size_id + '.0125map'
pfname: str = data_dir + 'out_' + size_id + '.0125ped'
ifname: str = data_dir + 'out_' + size_id + '.0125ids'

# Storing results in a dictionary (saving to a JSON file)
results[experiment_id][size_id] = {}
results[experiment_id][size_id]['fsize'] = []  # file size
results[experiment_id][size_id]['dbsize'] = []  # doc size in db
results[experiment_id][size_id]['map'] = {}  # map result data
results[experiment_id][size_id]['map']['fsize'] = []  # file size
results[experiment_id][size_id]['map']['dbsize'] = []  # doc size in db
results[experiment_id][size_id]['map']['time'] = []  # insertion time
results[experiment_id][size_id]['sample'] = {}  # sample result data
results[experiment_id][size_id]['sample']['fsize'] = []  # file size
results[experiment_id][size_id]['sample']['dbsize'] = []  # doc size in db
results[experiment_id][size_id]['sample']['time'] = []  # insertion time

# * Performing experiment N times and storing results
for i in range(N):
    print("i: " + str(i))
    print("Resetting database...")
    reset_db()
    print("Database reset operation successful.")
    print("Generating input files...")
    # * Generating input files
    with open(mfname, 'w') as mf, open(pfname, 'w') as pf, open(ifname,
                                                                'w') as iff:
        # Map file
        t = time.time()
        rfgen.random_0125_map(n=nsnps, outfile=mf, start_from_id=start_map)
        t = time.time() - t
        print("Generated map file\tTime: " + str(round(t, 3)) + "s\tSize: " +
              str(round(os.stat(mfname).st_size / (1024**2), 2)) + "MB")
        # Sample file
        t = time.time()
        rfgen.random_0125_samples(n=nsamples,
                                  map_size=nsnps,
                                  outfile=pf,
                                  start_from_id=start_sample)
        t = time.time() - t
        print("Generated samples file\tTime: " + str(round(t, 3)) +
              "s\tSize: " + str(round(os.stat(pfname).st_size /
                                      (1024**2), 2)) + "MB")
        # Id map file
        t = time.time()
        rfgen.id_mapping(n=nsamples, outfile=iff, first_sample_id=start_sample)
        t = time.time() - t
        print("Generated id mapping file\tTime: " + str(round(t, 3)) +
              "s\tSize: " + str(round(os.stat(ifname).st_size /
                                      (1024**2), 2)) + "MB")
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
    results[experiment_id][size_id]['map']['time'].append(t)
    # map_size = _MAPS + _MAPSNPS + _SNPS
    map_size = snpdb._db.command(
        "collstats",
        snpdb._config["MAPS_COLL"])['storageSize'] + snpdb._db.command(
            "collstats",
            snpdb._config["MAPSNPS_COLL"])['storageSize'] + snpdb._db.command(
                "collstats", snpdb._config["SNPS_COLL"])['storageSize']
    print("Imported map file\tTime: " + str(round(t, 3)) + "s\tSize:" +
          str(round(map_size / 1024**2, 2)) + "MB")
    # Importing sample file
    id_map: dict = {}
    # Linking samples to individuals in the database
    if ifname is not None:
        with open(ifname, "r") as f:
            for line in f:
                (sample, individual) = line.split()
                id_map[sample] = individual
    t = time.time()
    snpdb.import_samples(sample_reader=readers.Z125SampleReader(pfname),
                         map_name=experiment_id + size_id,
                         id_map=id_map,
                         report=False)
    t = time.time() - t
    results[experiment_id][size_id]['sample']['time'].append(t)
    # sample_size = _SAMPLES + _SNPBLOCKS + _INDS
    sample_size = snpdb._db.command(
        "collstats",
        snpdb._config["SAMPLES_COLL"])['storageSize'] + snpdb._db.command(
            "collstats", snpdb.
            _config["SNPBLOCKS_COLL"])['storageSize'] + snpdb._db.command(
                "collstats", snpdb._config["INDIVIDUALS_COLL"])['storageSize']
    print("Imported samples file\tTime: " + str(round(t, 3)) + "s\tSize:" +
          str(round(sample_size / 1024**2, 2)) + "MB")

    # Appending generated file sizes
    results[experiment_id][size_id]['map']['fsize'].append(
        os.stat(mfname).st_size)
    results[experiment_id][size_id]['sample']['fsize'].append(
        os.stat(pfname).st_size)
    results[experiment_id][size_id]['fsize'].append(
        os.stat(mfname).st_size + os.stat(pfname).st_size)
    # Appending stored document sizes from MongoDB
    results[experiment_id][size_id]['map']['dbsize'].append(map_size)
    results[experiment_id][size_id]['sample']['dbsize'].append(sample_size)
    results[experiment_id][size_id]['dbsize'].append(map_size + sample_size)

# %% [markdown]
# ### 1m SNPs

# %%
size_id = '1m'

# Printing experiment start information
print("Starting Experiment " + experiment_id + " (" + size_id +
      " SNPs) with N = " + str(N))

# Map and sample sizes
start_map = 1  # starting snp id
nsnps = 1000000  # number of snps in map
start_sample = 1  # starting sample id
nsamples = 1  # number of samples in file # ! Definir número de amostras
t = 0  # timing variable

# Filenames
mfname = data_dir + 'out_' + size_id + '.0125map'
pfname = data_dir + 'out_' + size_id + '.0125ped'
ifname = data_dir + 'out_' + size_id + '.0125ids'

# Storing results in a dictionary (saving to a JSON file)
results[experiment_id][size_id] = {}
results[experiment_id][size_id]['fsize'] = []  # file size
results[experiment_id][size_id]['dbsize'] = []  # doc size in db
results[experiment_id][size_id]['map'] = {}  # map result data
results[experiment_id][size_id]['map']['fsize'] = []  # file size
results[experiment_id][size_id]['map']['dbsize'] = []  # doc size in db
results[experiment_id][size_id]['map']['time'] = []  # insertion time
results[experiment_id][size_id]['sample'] = {}  # sample result data
results[experiment_id][size_id]['sample']['fsize'] = []  # file size
results[experiment_id][size_id]['sample']['dbsize'] = []  # doc size in db
results[experiment_id][size_id]['sample']['time'] = []  # insertion time

# * Performing experiment N times and storing results
for i in range(N):
    print("i: " + str(i))
    print("Resetting database...")
    reset_db()
    print("Database reset operation successful.")
    print("Generating input files...")
    # * Generating input files
    with open(mfname, 'w') as mf, open(pfname, 'w') as pf, open(ifname,
                                                                'w') as iff:
        # Map file
        t = time.time()
        rfgen.random_0125_map(n=nsnps, outfile=mf, start_from_id=start_map)
        t = time.time() - t
        print("Generated map file\tTime: " + str(round(t, 3)) + "s\tSize: " +
              str(round(os.stat(mfname).st_size / (1024**2), 2)) + "MB")
        # Sample file
        t = time.time()
        rfgen.random_0125_samples(n=nsamples,
                                  map_size=nsnps,
                                  outfile=pf,
                                  start_from_id=start_sample)
        t = time.time() - t
        print("Generated samples file\tTime: " + str(round(t, 3)) +
              "s\tSize: " + str(round(os.stat(pfname).st_size /
                                      (1024**2), 2)) + "MB")
        # Id map file
        t = time.time()
        rfgen.id_mapping(n=nsamples, outfile=iff, first_sample_id=start_sample)
        t = time.time() - t
        print("Generated id mapping file\tTime: " + str(round(t, 3)) +
              "s\tSize: " + str(round(os.stat(ifname).st_size /
                                      (1024**2), 2)) + "MB")
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
    results[experiment_id][size_id]['map']['time'].append(t)
    # map_size = _MAPS + _MAPSNPS + _SNPS
    map_size = snpdb._db.command(
        "collstats",
        snpdb._config["MAPS_COLL"])['storageSize'] + snpdb._db.command(
            "collstats",
            snpdb._config["MAPSNPS_COLL"])['storageSize'] + snpdb._db.command(
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
    snpdb.import_samples(sample_reader=readers.Z125SampleReader(pfname),
                         map_name=experiment_id + size_id,
                         id_map=id_map,
                         report=False)
    t = time.time() - t
    results[experiment_id][size_id]['sample']['time'].append(t)
    # sample_size = _SAMPLES + _SNPBLOCKS + _INDS
    sample_size = snpdb._db.command(
        "collstats",
        snpdb._config["SAMPLES_COLL"])['storageSize'] + snpdb._db.command(
            "collstats", snpdb.
            _config["SNPBLOCKS_COLL"])['storageSize'] + snpdb._db.command(
                "collstats", snpdb._config["INDIVIDUALS_COLL"])['storageSize']
    print("Imported samples file\tTime: " + str(round(t, 3)) + "s\tSize:" +
          str(round(sample_size / 1024**2, 2)) + "MB")

    # Appending generated file sizes
    results[experiment_id][size_id]['map']['fsize'].append(
        os.stat(mfname).st_size)
    results[experiment_id][size_id]['sample']['fsize'].append(
        os.stat(pfname).st_size)
    results[experiment_id][size_id]['fsize'].append(
        os.stat(mfname).st_size + os.stat(pfname).st_size)
    # Appending stored document sizes from MongoDB
    results[experiment_id][size_id]['map']['dbsize'].append(map_size)
    results[experiment_id][size_id]['sample']['dbsize'].append(sample_size)
    results[experiment_id][size_id]['dbsize'].append(map_size + sample_size)

# %% [markdown]
# ## 1.B - Arquivo Final Report (FR)

# %%
experiment_id = '1.B'
results[experiment_id] = {}

# %% [markdown]
# ### 100k SNPs

# %%
size_id = '100k'

# Printing experiment start information
print("Starting Experiment " + experiment_id + " (" + size_id +
      " SNPs) with N = " + str(N))

# Map and sample sizes
start_map = 1  # starting snp id
nsnps = 100000  # number of snps in map
start_sample = 1  # starting sample id
nsamples = 1  # number of samples in file # ! Definir número de amostras
t = 0  # timing variable

# Filenames
# ? FR -> one file for both map and samples
frfname: str = data_dir + 'out_' + size_id + '.fr'
ifname = data_dir + 'out_' + size_id + '.frids'

# Storing results in a dictionary (saving to a JSON file)
results[experiment_id][size_id] = {}
results[experiment_id][size_id]['fsize'] = []  # file size
results[experiment_id][size_id]['dbsize'] = []  # doc size in db
results[experiment_id][size_id]['map'] = {}  # map result data
results[experiment_id][size_id]['map']['dbsize'] = []  # doc size in db
results[experiment_id][size_id]['map']['time'] = []  # insertion time
results[experiment_id][size_id]['sample'] = {}  # sample result data
results[experiment_id][size_id]['sample']['dbsize'] = []  # doc size in db
results[experiment_id][size_id]['sample']['time'] = []  # insertion time

# * Performing experiment N times and storing results
for i in range(N):
    print("i: " + str(i))
    print("Resetting database...")
    reset_db()
    print("Database reset operation successful.")
    print("Generating input files...")
    # * Generating input files
    with open(frfname, 'w') as frf, open(ifname, 'w') as iff:
        t = time.time()
        rfgen.random_final_report(n=nsamples,
                                  map_size=nsnps,
                                  outfile=frf,
                                  start_snps_from_id=start_map,
                                  start_samples_from_id=start_sample)
        t = time.time() - t
        print("Generated Final Report file\tTime: " + str(round(t, 3)) +
              "s\tSize: " + str(round(os.stat(frfname).st_size /
                                      (1024**2), 2)) + "MB")
        # Id map file
        t = time.time()
        rfgen.id_mapping(n=nsamples, outfile=iff, first_sample_id=start_sample)
        t = time.time() - t
        print("Generated id mapping file\tTime: " + str(round(t, 3)) +
              "s\tSize: " + str(round(os.stat(ifname).st_size /
                                      (1024**2), 2)) + "MB")
    # * Inserting input files into db
    print("Inserting input files into database...")
    # Importing map file
    t = time.time()
    snpdb.import_map(map_reader=readers.FinalReportMapReader(frfname),
                     map_name=experiment_id + size_id,
                     force_create_new=True,
                     force_use_existing=False,
                     report=False)
    t = time.time() - t
    results[experiment_id][size_id]['map']['time'].append(t)
    # map_size = _MAPS + _MAPSNPS + _SNPS
    map_size = snpdb._db.command(
        "collstats",
        snpdb._config["MAPS_COLL"])['storageSize'] + snpdb._db.command(
            "collstats",
            snpdb._config["MAPSNPS_COLL"])['storageSize'] + snpdb._db.command(
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
        sample_reader=readers.FinalReportSampleReader(frfname),
        map_name=experiment_id + size_id,
        id_map=id_map,
        report=False)
    t = time.time() - t
    results[experiment_id][size_id]['sample']['time'].append(t)
    # sample_size = _SAMPLES + _SNPBLOCKS + _INDS
    sample_size = snpdb._db.command(
        "collstats",
        snpdb._config["SAMPLES_COLL"])['storageSize'] + snpdb._db.command(
            "collstats", snpdb.
            _config["SNPBLOCKS_COLL"])['storageSize'] + snpdb._db.command(
                "collstats", snpdb._config["INDIVIDUALS_COLL"])['storageSize']
    print("Imported samples file\tTime: " + str(round(t, 3)) + "s\tSize:" +
          str(round(sample_size / 1024**2, 2)) + "MB")

    # Appending generated file sizes
    results[experiment_id][size_id]['fsize'].append(os.stat(frfname).st_size)
    # Appending stored document sizes from MongoDB
    results[experiment_id][size_id]['map']['dbsize'].append(map_size)
    results[experiment_id][size_id]['sample']['dbsize'].append(sample_size)
    results[experiment_id][size_id]['dbsize'].append(map_size + sample_size)

# %% [markdown]
# ### 1m SNPs

# %%
size_id = '1m'

# Printing experiment start information
print("Starting Experiment " + experiment_id + " (" + size_id +
      " SNPs) with N = " + str(N))

# Map and sample sizes
start_map = 1  # starting snp id
nsnps = 1000000  # number of snps in map
start_sample = 1  # starting sample id
nsamples = 1  # number of samples in file # ! Definir número de amostras
t = 0  # timing variable

# Filenames
# ? FR -> one file for both map and samples
frfname = data_dir + 'out_' + size_id + '.fr'
ifname = data_dir + 'out_' + size_id + '.frids'

# Storing results in a dictionary (saving to a JSON file)
results[experiment_id][size_id] = {}
results[experiment_id][size_id]['fsize'] = []  # file size
results[experiment_id][size_id]['dbsize'] = []  # doc size in db
results[experiment_id][size_id]['map'] = {}  # map result data
results[experiment_id][size_id]['map']['dbsize'] = []  # doc size in db
results[experiment_id][size_id]['map']['time'] = []  # insertion time
results[experiment_id][size_id]['sample'] = {}  # sample result data
results[experiment_id][size_id]['sample']['dbsize'] = []  # doc size in db
results[experiment_id][size_id]['sample']['time'] = []  # insertion time

# * Performing experiment N times and storing results
for i in range(N):
    print("i: " + str(i))
    print("Resetting database...")
    reset_db()
    print("Database reset operation successful.")
    print("Generating input files...")
    # * Generating input files
    with open(frfname, 'w') as frf, open(ifname, 'w') as iff:
        t = time.time()
        rfgen.random_final_report(n=nsamples,
                                  map_size=nsnps,
                                  outfile=frf,
                                  start_snps_from_id=start_map,
                                  start_samples_from_id=start_sample)
        t = time.time() - t
        print("Generated Final Report file\tTime: " + str(round(t, 3)) +
              "s\tSize: " + str(round(os.stat(frfname).st_size /
                                      (1024**2), 2)) + "MB")
        # Id map file
        t = time.time()
        rfgen.id_mapping(n=nsamples, outfile=iff, first_sample_id=start_sample)
        t = time.time() - t
        print("Generated id mapping file\tTime: " + str(round(t, 3)) +
              "s\tSize: " + str(round(os.stat(ifname).st_size /
                                      (1024**2), 2)) + "MB")
    # * Inserting input files into db
    print("Inserting input files into database...")
    # Importing map file
    t = time.time()
    snpdb.import_map(map_reader=readers.FinalReportMapReader(frfname),
                     map_name=experiment_id + size_id,
                     force_create_new=True,
                     force_use_existing=False,
                     report=False)
    t = time.time() - t
    results[experiment_id][size_id]['map']['time'].append(t)
    # map_size = _MAPS + _MAPSNPS + _SNPS
    map_size = snpdb._db.command(
        "collstats",
        snpdb._config["MAPS_COLL"])['storageSize'] + snpdb._db.command(
            "collstats",
            snpdb._config["MAPSNPS_COLL"])['storageSize'] + snpdb._db.command(
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
        sample_reader=readers.FinalReportSampleReader(frfname),
        map_name=experiment_id + size_id,
        id_map=id_map,
        report=False)
    t = time.time() - t
    results[experiment_id][size_id]['sample']['time'].append(t)
    # sample_size = _SAMPLES + _SNPBLOCKS + _INDS
    sample_size = snpdb._db.command(
        "collstats",
        snpdb._config["SAMPLES_COLL"])['storageSize'] + snpdb._db.command(
            "collstats", snpdb.
            _config["SNPBLOCKS_COLL"])['storageSize'] + snpdb._db.command(
                "collstats", snpdb._config["INDIVIDUALS_COLL"])['storageSize']
    print("Imported samples file\tTime: " + str(round(t, 3)) + "s\tSize:" +
          str(round(sample_size / 1024**2, 2)) + "MB")

    # Appending generated file sizes
    results[experiment_id][size_id]['fsize'].append(os.stat(frfname).st_size)
    # Appending stored document sizes from MongoDB
    results[experiment_id][size_id]['map']['dbsize'].append(map_size)
    results[experiment_id][size_id]['sample']['dbsize'].append(sample_size)
    results[experiment_id][size_id]['dbsize'].append(map_size + sample_size)

# %% [markdown]
# ## 1.C - Arquivo VCF

# %%
experiment_id = '1.C'
results[experiment_id] = {}

# %% [markdown]
# ### 100k SNPs

# %%
size_id = '100k'

# Printing experiment start information
print("Starting Experiment " + experiment_id + " (" + size_id +
      " SNPs) with N = " + str(N))

# Map and sample sizes
start_map = 1  # starting snp id
nsnps = 100000  # number of snps in map
start_sample = 1  # starting sample id
nsamples = 1  # number of samples in file # ! Definir número de amostras
t = 0  # timing variable

# Filenames
# ? VCF -> one file for both map and samples
vcffname: str = data_dir + 'out_' + size_id + '.vcf'
ifname = data_dir + 'out_' + size_id + '.vcfids'

# Storing results in a dictionary (saving to a JSON file)
results[experiment_id][size_id] = {}
results[experiment_id][size_id]['fsize'] = []  # file size
results[experiment_id][size_id]['dbsize'] = []  # doc size in db
results[experiment_id][size_id]['map'] = {}  # map result data
results[experiment_id][size_id]['map']['dbsize'] = []  # doc size in db
results[experiment_id][size_id]['map']['time'] = []  # insertion time
results[experiment_id][size_id]['sample'] = {}  # sample result data
results[experiment_id][size_id]['sample']['dbsize'] = []  # doc size in db
results[experiment_id][size_id]['sample']['time'] = []  # insertion time

# * Performing experiment N times and storing results
for i in range(N):
    print("i: " + str(i))
    print("Resetting database...")
    reset_db()
    print("Database reset operation successful.")
    print("Generating input files...")
    # * Generating input files
    with open(vcffname, 'w') as vcf, open(ifname, 'w') as iff:
        t = time.time()
        rfgen.random_vcf(n=nsamples,
                         map_size=nsnps,
                         outfile=vcf,
                         start_snps_from_id=start_map,
                         start_samples_from_id=start_sample)
        t = time.time() - t
        print("Generated VCF file\tTime: " + str(round(t, 3)) + "s\tSize: " +
              str(round(os.stat(vcffname).st_size / (1024**2), 2)) + "MB")
        # Id map file
        t = time.time()
        rfgen.id_mapping(n=nsamples, outfile=iff, first_sample_id=start_sample)
        t = time.time() - t
        print("Generated id mapping file\tTime: " + str(round(t, 3)) +
              "s\tSize: " + str(round(os.stat(ifname).st_size /
                                      (1024**2), 2)) + "MB")
    # * Inserting input files into db
    print("Inserting input files into database...")
    # Importing map file
    t = time.time()
    snpdb.import_map(map_reader=readers.VcfMapReader(vcffname),
                     map_name=experiment_id + size_id,
                     force_create_new=True,
                     force_use_existing=False,
                     report=False)
    t = time.time() - t
    results[experiment_id][size_id]['map']['time'].append(t)
    # map_size = _MAPS + _MAPSNPS + _SNPS
    map_size = snpdb._db.command(
        "collstats",
        snpdb._config["MAPS_COLL"])['storageSize'] + snpdb._db.command(
            "collstats",
            snpdb._config["MAPSNPS_COLL"])['storageSize'] + snpdb._db.command(
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
    snpdb.import_samples(sample_reader=readers.VcfSampleReader(vcffname),
                         map_name=experiment_id + size_id,
                         id_map=id_map,
                         report=False)
    t = time.time() - t
    results[experiment_id][size_id]['sample']['time'].append(t)
    # sample_size = _SAMPLES + _SNPBLOCKS + _INDS
    sample_size = snpdb._db.command(
        "collstats",
        snpdb._config["SAMPLES_COLL"])['storageSize'] + snpdb._db.command(
            "collstats", snpdb.
            _config["SNPBLOCKS_COLL"])['storageSize'] + snpdb._db.command(
                "collstats", snpdb._config["INDIVIDUALS_COLL"])['storageSize']
    print("Imported samples file\tTime: " + str(round(t, 3)) + "s\tSize:" +
          str(round(sample_size / 1024**2, 2)) + "MB")

    # Appending generated file sizes
    results[experiment_id][size_id]['fsize'].append(os.stat(vcffname).st_size)
    # Appending stored document sizes from MongoDB
    results[experiment_id][size_id]['map']['dbsize'].append(map_size)
    results[experiment_id][size_id]['sample']['dbsize'].append(sample_size)
    results[experiment_id][size_id]['dbsize'].append(map_size + sample_size)

# %% [markdown]
# ### 1m SNPs

# %%
size_id = '1m'

# Printing experiment start information
print("Starting Experiment " + experiment_id + " (" + size_id +
      " SNPs) with N = " + str(N))

# Map and sample sizes
start_map = 1  # starting snp id
nsnps = 1000000  # number of snps in map
start_sample = 1  # starting sample id
nsamples = 1  # number of samples in file # ! Definir número de amostras
t = 0  # timing variable

# Filenames
# ? VCF -> one file for both map and samples
vcffname = data_dir + 'out_' + size_id + '.vcf'
ifname = data_dir + 'out_' + size_id + '.vcfids'

# Storing results in a dictionary (saving to a JSON file)
results[experiment_id][size_id] = {}
results[experiment_id][size_id]['fsize'] = []  # file size
results[experiment_id][size_id]['dbsize'] = []  # doc size in db
results[experiment_id][size_id]['map'] = {}  # map result data
results[experiment_id][size_id]['map']['dbsize'] = []  # doc size in db
results[experiment_id][size_id]['map']['time'] = []  # insertion time
results[experiment_id][size_id]['sample'] = {}  # sample result data
results[experiment_id][size_id]['sample']['dbsize'] = []  # doc size in db
results[experiment_id][size_id]['sample']['time'] = []  # insertion time

# * Performing experiment N times and storing results
for i in range(N):
    print("i: " + str(i))
    print("Resetting database...")
    reset_db()
    print("Database reset operation successful.")
    print("Generating input files...")
    # * Generating input files
    with open(vcffname, 'w') as vcf, open(ifname, 'w') as iff:
        t = time.time()
        rfgen.random_vcf(n=nsamples,
                         map_size=nsnps,
                         outfile=vcf,
                         start_snps_from_id=start_map,
                         start_samples_from_id=start_sample)
        t = time.time() - t
        print("Generated VCF file\tTime: " + str(round(t, 3)) + "s\tSize: " +
              str(round(os.stat(vcffname).st_size / (1024**2), 2)) + "MB")
        # Id map file
        t = time.time()
        rfgen.id_mapping(n=nsamples, outfile=iff, first_sample_id=start_sample)
        t = time.time() - t
        print("Generated id mapping file\tTime: " + str(round(t, 3)) +
              "s\tSize: " + str(round(os.stat(ifname).st_size /
                                      (1024**2), 2)) + "MB")
    # * Inserting input files into db
    print("Inserting input files into database...")
    # Importing map file
    t = time.time()
    snpdb.import_map(map_reader=readers.VcfMapReader(vcffname),
                     map_name=experiment_id + size_id,
                     force_create_new=True,
                     force_use_existing=False,
                     report=False)
    t = time.time() - t
    results[experiment_id][size_id]['map']['time'].append(t)
    # map_size = _MAPS + _MAPSNPS + _SNPS
    map_size = snpdb._db.command(
        "collstats",
        snpdb._config["MAPS_COLL"])['storageSize'] + snpdb._db.command(
            "collstats",
            snpdb._config["MAPSNPS_COLL"])['storageSize'] + snpdb._db.command(
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
    snpdb.import_samples(sample_reader=readers.VcfSampleReader(vcffname),
                         map_name=experiment_id + size_id,
                         id_map=id_map,
                         report=False)
    t = time.time() - t
    results[experiment_id][size_id]['sample']['time'].append(t)
    # sample_size = _SAMPLES + _SNPBLOCKS + _INDS
    sample_size = snpdb._db.command(
        "collstats",
        snpdb._config["SAMPLES_COLL"])['storageSize'] + snpdb._db.command(
            "collstats", snpdb.
            _config["SNPBLOCKS_COLL"])['storageSize'] + snpdb._db.command(
                "collstats", snpdb._config["INDIVIDUALS_COLL"])['storageSize']
    print("Imported samples file\tTime: " + str(round(t, 3)) + "s\tSize:" +
          str(round(sample_size / 1024**2, 2)) + "MB")

    # Appending generated file sizes
    results[experiment_id][size_id]['fsize'].append(os.stat(vcffname).st_size)
    # Appending stored document sizes from MongoDB
    results[experiment_id][size_id]['map']['dbsize'].append(map_size)
    results[experiment_id][size_id]['sample']['dbsize'].append(sample_size)
    results[experiment_id][size_id]['dbsize'].append(map_size + sample_size)

# %% [markdown]
# ## 1.D - Arquivo PLink

# %%
experiment_id = '1.D'
results[experiment_id] = {}

# %% [markdown]
# ### 100k SNPs

# %%
size_id = '100k'

# Printing experiment start information
print("Starting Experiment " + experiment_id + " (" + size_id +
      " SNPs) with N = " + str(N))

# Map and sample sizes
start_map = 1  # starting snp id
nsnps = 100000  # number of snps in map
start_sample = 1  # starting sample id
nsamples = 1  # number of samples in file # ! Definir número de amostras
t = 0  # timing variable

# Filenames
mfname = data_dir + 'out_' + size_id + '.plmap'
pfname = data_dir + 'out_' + size_id + '.plped'
ifname = data_dir + 'out_' + size_id + '.plids'

# Storing results in a dictionary (saving to a JSON file)
results[experiment_id][size_id] = {}
results[experiment_id][size_id]['fsize'] = []  # file size
results[experiment_id][size_id]['dbsize'] = []  # doc size in db
results[experiment_id][size_id]['map'] = {}  # map result data
results[experiment_id][size_id]['map']['fsize'] = []  # file size
results[experiment_id][size_id]['map']['dbsize'] = []  # doc size in db
results[experiment_id][size_id]['map']['time'] = []  # insertion time
results[experiment_id][size_id]['sample'] = {}  # sample result data
results[experiment_id][size_id]['sample']['fsize'] = []  # file size
results[experiment_id][size_id]['sample']['dbsize'] = []  # doc size in db
results[experiment_id][size_id]['sample']['time'] = []  # insertion time

# * Performing experiment N times and storing results
for i in range(N):
    print("i: " + str(i))
    print("Resetting database...")
    reset_db()
    print("Database reset operation successful.")
    print("Generating input files...")
    # * Generating input files
    with open(mfname, 'w') as mf, open(pfname, 'w') as pf, open(ifname,
                                                                'w') as iff:
        # Map file
        t = time.time()
        rfgen.random_plink_map(n=nsnps, outfile=mf, start_from_id=start_map)
        t = time.time() - t
        print("Generated map file\tTime: " + str(round(t, 3)) + "s\tSize: " +
              str(round(os.stat(mfname).st_size / (1024**2), 2)) + "MB")
        # Sample file
        t = time.time()
        rfgen.random_plink_samples(n=nsamples,
                                   map_size=nsnps,
                                   outfile=pf,
                                   start_from_id=start_sample)
        t = time.time() - t
        print("Generated samples file\tTime: " + str(round(t, 3)) +
              "s\tSize: " + str(round(os.stat(pfname).st_size /
                                      (1024**2), 2)) + "MB")
        # Id map file
        t = time.time()
        rfgen.id_mapping(n=nsamples, outfile=iff, first_sample_id=start_sample)
        t = time.time() - t
        print("Generated id mapping file\tTime: " + str(round(t, 3)) +
              "s\tSize: " + str(round(os.stat(ifname).st_size /
                                      (1024**2), 2)) + "MB")
    # * Inserting input files into db
    print("Inserting input files into database...")
    # Importing map file
    t = time.time()
    snpdb.import_map(map_reader=readers.PlinkMapReader(mfname),
                     map_name=experiment_id + size_id,
                     force_create_new=True,
                     force_use_existing=False,
                     report=False)
    t = time.time() - t
    results[experiment_id][size_id]['map']['time'].append(t)
    # map_size = _MAPS + _MAPSNPS + _SNPS
    map_size = snpdb._db.command(
        "collstats",
        snpdb._config["MAPS_COLL"])['storageSize'] + snpdb._db.command(
            "collstats",
            snpdb._config["MAPSNPS_COLL"])['storageSize'] + snpdb._db.command(
                "collstats", snpdb._config["SNPS_COLL"])['storageSize']
    print("Imported map file\tTime: " + str(round(t, 3)) + "s\tSize:" +
          str(round(map_size / 1024**2, 2)) + "MB")
    # Importing sample file
    id_map = {}
    if ifname is not None:
        with open(ifname, "r") as f:
            for line in f:
                (sample, individual) = line.split()
                id_map[sample] = individual
    t = time.time()
    snpdb.import_samples(sample_reader=readers.PlinkSampleReader(pfname),
                         map_name=experiment_id + size_id,
                         id_map=id_map,
                         report=False)
    t = time.time() - t
    results[experiment_id][size_id]['sample']['time'].append(t)
    # sample_size = _SAMPLES + _SNPBLOCKS + _INDS
    sample_size = snpdb._db.command(
        "collstats",
        snpdb._config["SAMPLES_COLL"])['storageSize'] + snpdb._db.command(
            "collstats", snpdb.
            _config["SNPBLOCKS_COLL"])['storageSize'] + snpdb._db.command(
                "collstats", snpdb._config["INDIVIDUALS_COLL"])['storageSize']
    print("Imported samples file\tTime: " + str(round(t, 3)) + "s\tSize:" +
          str(round(sample_size / 1024**2, 2)) + "MB")

    # Appending generated file sizes
    results[experiment_id][size_id]['map']['fsize'].append(
        os.stat(mfname).st_size)
    results[experiment_id][size_id]['sample']['fsize'].append(
        os.stat(pfname).st_size)
    results[experiment_id][size_id]['fsize'].append(
        os.stat(mfname).st_size + os.stat(pfname).st_size)
    # Appending stored document sizes from MongoDB
    results[experiment_id][size_id]['map']['dbsize'].append(map_size)
    results[experiment_id][size_id]['sample']['dbsize'].append(sample_size)
    results[experiment_id][size_id]['dbsize'].append(map_size + sample_size)

# %% [markdown]
# ### 1m SNPs

# %%
size_id = '1m'

# Printing experiment start information
print("Starting Experiment " + experiment_id + " (" + size_id +
      " SNPs) with N = " + str(N))

# Map and sample sizes
start_map = 1  # starting snp id
nsnps = 1000000  # number of snps in map
start_sample = 1  # starting sample id
nsamples = 1  # number of samples in file # ! Definir número de amostras
t = 0  # timing variable

# Filenames
mfname = data_dir + 'out_' + size_id + '.plmap'
pfname = data_dir + 'out_' + size_id + '.plped'
ifname = data_dir + 'out_' + size_id + '.plids'

# Storing results in a dictionary (saving to a JSON file)
results[experiment_id][size_id] = {}
results[experiment_id][size_id]['fsize'] = []  # file size
results[experiment_id][size_id]['dbsize'] = []  # doc size in db
results[experiment_id][size_id]['map'] = {}  # map result data
results[experiment_id][size_id]['map']['fsize'] = []  # file size
results[experiment_id][size_id]['map']['dbsize'] = []  # doc size in db
results[experiment_id][size_id]['map']['time'] = []  # insertion time
results[experiment_id][size_id]['sample'] = {}  # sample result data
results[experiment_id][size_id]['sample']['fsize'] = []  # file size
results[experiment_id][size_id]['sample']['dbsize'] = []  # doc size in db
results[experiment_id][size_id]['sample']['time'] = []  # insertion time

# * Performing experiment N times and storing results
for i in range(N):
    print("i: " + str(i))
    print("Resetting database...")
    reset_db()
    print("Database reset operation successful.")
    print("Generating input files...")
    # * Generating input files
    with open(mfname, 'w') as mf, open(pfname, 'w') as pf, open(ifname,
                                                                'w') as iff:
        # Map file
        t = time.time()
        rfgen.random_plink_map(n=nsnps, outfile=mf, start_from_id=start_map)
        t = time.time() - t
        print("Generated map file\tTime: " + str(round(t, 3)) + "s\tSize: " +
              str(round(os.stat(mfname).st_size / (1024**2), 2)) + "MB")
        # Sample file
        t = time.time()
        rfgen.random_plink_samples(n=nsamples,
                                   map_size=nsnps,
                                   outfile=pf,
                                   start_from_id=start_sample)
        t = time.time() - t
        print("Generated samples file\tTime: " + str(round(t, 3)) +
              "s\tSize: " + str(round(os.stat(pfname).st_size /
                                      (1024**2), 2)) + "MB")
        # Id map file
        t = time.time()
        rfgen.id_mapping(n=nsamples, outfile=iff, first_sample_id=start_sample)
        t = time.time() - t
        print("Generated id mapping file\tTime: " + str(round(t, 3)) +
              "s\tSize: " + str(round(os.stat(ifname).st_size /
                                      (1024**2), 2)) + "MB")
    # * Inserting input files into db
    print("Inserting input files into database...")
    # Importing map file
    t = time.time()
    snpdb.import_map(map_reader=readers.PlinkMapReader(mfname),
                     map_name=experiment_id + size_id,
                     force_create_new=True,
                     force_use_existing=False,
                     report=False)
    t = time.time() - t
    results[experiment_id][size_id]['map']['time'].append(t)
    # map_size = _MAPS + _MAPSNPS + _SNPS
    map_size = snpdb._db.command(
        "collstats",
        snpdb._config["MAPS_COLL"])['storageSize'] + snpdb._db.command(
            "collstats",
            snpdb._config["MAPSNPS_COLL"])['storageSize'] + snpdb._db.command(
                "collstats", snpdb._config["SNPS_COLL"])['storageSize']
    print("Imported map file\tTime: " + str(round(t, 3)) + "s\tSize:" +
          str(round(map_size / 1024**2, 2)) + "MB")
    # Importing sample file
    id_map = {}
    if ifname is not None:
        with open(ifname, "r") as f:
            for line in f:
                (sample, individual) = line.split()
                id_map[sample] = individual
    t = time.time()
    snpdb.import_samples(sample_reader=readers.PlinkSampleReader(pfname),
                         map_name=experiment_id + size_id,
                         id_map=id_map,
                         report=False)
    t = time.time() - t
    results[experiment_id][size_id]['sample']['time'].append(t)
    # sample_size = _SAMPLES + _SNPBLOCKS + _INDS
    sample_size = snpdb._db.command(
        "collstats",
        snpdb._config["SAMPLES_COLL"])['storageSize'] + snpdb._db.command(
            "collstats", snpdb.
            _config["SNPBLOCKS_COLL"])['storageSize'] + snpdb._db.command(
                "collstats", snpdb._config["INDIVIDUALS_COLL"])['storageSize']
    print("Imported samples file\tTime: " + str(round(t, 3)) + "s\tSize:" +
          str(round(sample_size / 1024**2, 2)) + "MB")

    # Appending generated file sizes
    results[experiment_id][size_id]['map']['fsize'].append(
        os.stat(mfname).st_size)
    results[experiment_id][size_id]['sample']['fsize'].append(
        os.stat(pfname).st_size)
    results[experiment_id][size_id]['fsize'].append(
        os.stat(mfname).st_size + os.stat(pfname).st_size)
    # Appending stored document sizes from MongoDB
    results[experiment_id][size_id]['map']['dbsize'].append(map_size)
    results[experiment_id][size_id]['sample']['dbsize'].append(sample_size)
    results[experiment_id][size_id]['dbsize'].append(map_size + sample_size)

# %% [markdown]
# ## 1.E - Arquivo FastQ
# Utilizando um arquivo FastQ disponível no servidor de execução.

# %%
experiment_id = '1.E'
results[experiment_id] = {}

# %%
# Printing experiment start information
print("Starting Experiment " + experiment_id + " with N = " + str(N) + " ...")

# Storing results in a dictionary (saving to a JSON file)
results[experiment_id]['fsize'] = []  # file size
results[experiment_id]['dbsize'] = []  # doc size in db
results[experiment_id]['time'] = []  # insertion time

for i in range(N):
    print("i: " + str(i))
    print("Resetting database...")
    reset_db()
    # Selecting file # ! Verify which files to use
    fqfname = fastq_dir_1 + 'SH.71992.AP.01.1.fastq'
    # * Inserting input files into db
    # Importing sequencing file
    with open(fqfname, "rb") as fqf:
        print(fqfname + "\tSize: " +
              str(round(os.stat(fqfname).st_size / (1024**2), 2)) + "MB")
        t = time.time()
        snpdb.insert_file(file=fqf, individual_id="IND1")
        t = time.time() - t
        file_size = snpdb._db.command("collstats", "fs.chunks")["storageSize"]
        # file_size = snpdb._db["fs.files"].find_one()["length"]
        print("Imported sequencing file\tTime: " + str(round(t, 3)) +
              "s\tSize:" + str(round(file_size / 1024**2, 2)) + "MB")

        # Appending experiment results
        results[experiment_id]['fsize'].append(os.stat(fqfname).st_size)
        results[experiment_id]['dbsize'].append(file_size)
        results[experiment_id]['time'].append(t)

# %% [markdown]
# ## 1.F - Arquivo de Mídia
# Gerando uma imagem aleatoriamente e inserindo na base de dados.

# %%
experiment_id = '1.F'
results[experiment_id] = {}

# %%
# Printing experiment start information
print("Starting Experiment " + experiment_id + " with N = " + str(N) + " ...")

im_res = (800, 600)  # Image resolution
imfname = data_dir + 'out_image.jpg'  # Image filename

# Storing results in a dictionary (saving to a JSON file)
results[experiment_id]['fsize'] = []  # file size
results[experiment_id]['dbsize'] = []  # doc size in db
results[experiment_id]['time'] = []  # insertion time

for i in range(N):
    print("i: " + str(i))
    print("Resetting database...")
    reset_db()
    print("Generating input file...")
    # * Generating input file
    t = time.time()
    im_arr = np.random.rand(im_res[0], im_res[1], 3) * 255
    im_out = Image.fromarray(im_arr.astype('uint8')).convert('RGB')
    im_out.save(imfname)
    t = time.time() - t
    print("Generated media file\tTime: " + str(round(t, 3)) + "s\tSize: " +
          str(round(os.stat(imfname).st_size / (1024**2), 2)) + "MB")
    # * Inserting input files into db
    # Importing media file
    with open(imfname, "rb") as imf:
        t = time.time()
        snpdb.insert_file(file=imf, individual_id="IND1")
        t = time.time() - t
        file_size = snpdb._db.command("collstats", "fs.chunks")["storageSize"]
        # file_size = snpdb._db["fs.files"].find_one()["length"]
        print("Imported media file\tTime: " + str(round(t, 3)) + "s\tSize:" +
              str(round(file_size / 1024**2, 2)) + "MB")

        # Appending experiment results
        results[experiment_id]['fsize'].append(os.stat(imfname).st_size)
        results[experiment_id]['dbsize'].append(file_size)
        results[experiment_id]['time'].append(t)

# %% [markdown]
# # Saving results

# %%
# Saving results to a JSON file
jfname = data_dir + 'exp1results.json'
with open(jfname, 'w') as jf:
    json.dump(results, jf)
