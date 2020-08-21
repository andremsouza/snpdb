# %% [markdown]
#  # Experimentos SNPDB

# %% [markdown]
#  ## Imports e inicialização

# %%
import importlib
import json
import os
import snpdb
import subprocess
import readers
import testing.random_file_generator as rfgen
import time
from typing import Callable, Tuple
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
def reset_db() -> None:
    """Reset MongoDB database for experiments."""
    snpdb._client.drop_database(snpdb._config["DB_NAME"])
    subprocess.run("""mongo --eval "load('./mongo_setup.js')" """,
                   stdout=subprocess.DEVNULL,
                   stderr=subprocess.STDOUT,
                   shell=True)
    importlib.reload(snpdb)


def generate_random_file(filename: str,
                         file_type: str,
                         verbose: bool = False,
                         **kwargs) -> Tuple[float, float]:
    """Generates random files for snpdb experiments

    Args:
        filename (str): File name or path for output.
        file_type (str): Type of file to be generated. Must be in the
            following values:
            '0125_map', '0125_samples', 'final_report', 'vcf', 'plink_map',
            'plink_samples', 'id_mapping'.
        verbose (bool, optional): If True, prints timing and size information
            about the generated file. Defaults to False.
        **kwargs (Any): Keyword arguments for file generation function.
            'outfile' may be omitted.

    Returns:
        Tuple[float, float]: Returns the function generation time, and its
            size on disk
    """
    func: Callable = {
        '0125_map': rfgen.random_0125_map,
        '0125_samples': rfgen.random_0125_samples,
        'final_report': rfgen.random_final_report,
        'vcf': rfgen.random_vcf,
        'plink_map': rfgen.random_plink_map,
        'plink_samples': rfgen.random_plink_samples,
        'id_mapping': rfgen.id_mapping
    }[file_type]
    t: float = time.time()
    with open(filename, 'w') as f:
        func(outfile=f, **kwargs)
    t = time.time() - t
    if verbose:
        print("Generated " + file_type + "\tTime: " + str(round(t, 3)) +
              "s\tSize: " +
              str(round(os.stat(filename).st_size / (1024**2), 2)) + "MB")
    return t, os.stat(filename).st_size


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
# nspns_list: list = [100000, 1000000]
nsnps_list: list = [100000]
nsamples_list: list = [10, 100, 1000, 10000, 100000, 1000000]

# %% [markdown]
# ## 2.1 - Análise de performance de casos específicos
# Analisando performance considerando os casos com menor e maior tempo de
# importação da seção 1, assim como o caso G (tempo de importação para todos
# os tipos de arquivos).

# Tipos de arquivos analisados, segundo resultados do experimento 1:
# - A - Arquivo 0125 {'100k': 5.727, '1m': 60.724}
# - B - Arquivo VCF {'100k': 9.894, '1m': 72.766}
# - C - Todos os tipos de arquivos, excluindo arquivos FastQ e de mídia
# ! Verificar inclusão de arquivos de sequenciamento e de mídia no caso C

# %% [markdown]
# ### 2.1.A - Arquivo 0125

# %%
experiment_id: str = '2.1.A'
results[experiment_id] = {}

for nsnps in nsnps_list:
    size_id = {100000: '100k', 1000000: '1m'}[nsnps]
    for nsamples in nsamples_list:
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
            id_map: dict = {}
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

# %% [markdown]
# ### 2.1.B - Arquivo VCF

# %%
experiment_id = '2.1.B'
results[experiment_id] = {}

for nsnps in nsnps_list:
    size_id = {100000: '100k', 1000000: '1m'}[nsnps]
    for nsamples in nsamples_list:
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

        # Map and sample sizes
        # start_map = 1  # starting snp id
        # start_sample = 1  # starting sample id
        t = 0  # timing variable

        # Filenames
        # ? VCF -> one file for both map and samples
        vcffname: str = data_dir + 'out_' + size_id + '_' + count_id + '.vcf'
        ifname = data_dir + 'out_' + size_id + '_' + count_id + '.vcfids'

        # Storing results in a dictionary (saving to a JSON file)
        results[experiment_id][size_id] = {}
        results[experiment_id][size_id][count_id] = {}
        results[experiment_id][size_id][count_id]['fsize'] = []  # file size
        results[experiment_id][size_id][count_id]['dbsize'] = []
        # map result data
        results[experiment_id][size_id][count_id]['map'] = {}
        results[experiment_id][size_id][count_id]['map']['dbsize'] = []
        results[experiment_id][size_id][count_id]['map']['time'] = []
        # sample result data
        results[experiment_id][size_id][count_id]['sample'] = {}
        results[experiment_id][size_id][count_id]['sample']['dbsize'] = []
        results[experiment_id][size_id][count_id]['sample']['time'] = []

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
            snpdb.import_map(map_reader=readers.VcfMapReader(vcffname),
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
                sample_reader=readers.VcfSampleReader(vcffname),
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

            # Appending generated file sizes
            results[experiment_id][size_id][count_id]['fsize'].append(
                os.stat(vcffname).st_size)
            # Appending stored document sizes from MongoDB
            results[experiment_id][size_id][count_id]['map']['dbsize'].append(
                map_size)
            results[experiment_id][size_id][count_id]['sample'][
                'dbsize'].append(sample_size)
            results[experiment_id][size_id][count_id]['dbsize'].append(
                map_size + sample_size)

# %% [markdown]
# ### 2.1.C - Todos os tipos de arquivos

# %%
experiment_id = '2.1.C'
results[experiment_id] = {}

fqfname: str = fastq_dir_1 + 'SH.71992.AP.01.1.fastq'
imfname: str = data_dir + 'out_image.jpg'
im_res: tuple = (800, 600)

for nsnps in nsnps_list:
    size_id = {100000: '100k', 1000000: '1m'}[nsnps]
    for nsamples in nsamples_list:
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

        # Map and sample sizes
        # start_map = 1  # starting snp id
        # start_sample = 1  # starting sample id
        t = 0  # timing variable
        t_tmp: float = 0  # temporary timing variable

        # Filenames
        mfnames: dict = {
            'Z125': data_dir + 'out_' + size_id + '_' + count_id + '.0125map',
            'PL': data_dir + 'out_' + size_id + '_' + count_id + '.plmap'
        }
        pfnames: dict = {
            'Z125': data_dir + 'out_' + size_id + '_' + count_id + '.0125ped',
            'PL': data_dir + 'out_' + size_id + '_' + count_id + '.plped'
        }
        frfname = data_dir + 'out_' + size_id + '_' + count_id + '.fr'
        vcffname = data_dir + 'out_' + size_id + '_' + count_id + '.vcf'
        ifnames: dict = {
            'Z125': data_dir + 'out_' + size_id + '_' + count_id + '.0125ids',
            'FR': data_dir + 'out_' + size_id + '_' + count_id + '.frids',
            'VCF': data_dir + 'out_' + size_id + '_' + count_id + '.vcfids',
            'PL': data_dir + 'out_' + size_id + '_' + count_id + '.plids'
        }

        # Storing results in a dictionary (saving to a JSON file)
        results[experiment_id][size_id] = {}
        results[experiment_id][size_id][count_id] = {}
        results[experiment_id][size_id][count_id]['fsize'] = []  # file size
        results[experiment_id][size_id][count_id]['dbsize'] = []
        results[experiment_id][size_id][count_id]['time'] = []

        # * Performing experiment N times and storing results
        for i in range(N):
            print("i: " + str(i))
            print("Resetting database...")
            reset_db()
            print("Database reset operation successful.")
            # * Inserting input files into db
            print("Inserting input files into database...")
            # * Z125 Files
            # Importing map file
            t_tmp = time.time()
            snpdb.import_map(map_reader=readers.Z125MapReader(mfnames['Z125']),
                             map_name=experiment_id + size_id + 'Z125',
                             force_create_new=True,
                             force_use_existing=False,
                             report=False)
            t += time.time() - t_tmp
            # Importing sample file
            id_map = {}
            # Linking samples to individuals in the database
            if ifnames.get('Z125') is not None:
                with open(ifnames['Z125'], "r") as f:
                    for line in f:
                        (sample, individual) = line.split()
                        id_map[sample] = individual
            t_tmp = time.time()
            snpdb.import_samples(sample_reader=readers.Z125SampleReader(
                pfnames['Z125']),
                                 map_name=experiment_id + size_id + 'Z125',
                                 id_map=id_map,
                                 report=False)
            t += time.time() - t_tmp
            # * FR Files
            # ! Falta gerar um arquivo FR (100k SNPs, 1m samples)
            # Importing map file
            t_tmp = time.time()
            snpdb.import_map(map_reader=readers.FinalReportMapReader(frfname),
                             map_name=experiment_id + size_id + 'FR',
                             force_create_new=True,
                             force_use_existing=False,
                             report=False)
            t += time.time() - t_tmp
            # Importing sample file
            id_map = {}
            # Linking samples to individuals in the database
            if ifnames.get('FR') is not None:
                with open(ifnames['FR'], "r") as f:
                    for line in f:
                        (sample, individual) = line.split()
                        id_map[sample] = individual
            t_tmp = time.time()
            snpdb.import_samples(
                sample_reader=readers.FinalReportSampleReader(frfname),
                map_name=experiment_id + size_id + 'FR',
                id_map=id_map,
                report=False)
            t += time.time() - t_tmp
            # * VCF Files
            # Importing map file
            t_tmp = time.time()
            snpdb.import_map(map_reader=readers.VcfMapReader(vcffname),
                             map_name=experiment_id + size_id + 'VCF',
                             force_create_new=True,
                             force_use_existing=False,
                             report=False)
            t += time.time() - t_tmp
            # Importing sample file
            id_map = {}
            # Linking samples to individuals in the database
            if ifnames.get('VCF') is not None:
                with open(ifnames['VCF'], "r") as f:
                    for line in f:
                        (sample, individual) = line.split()
                        id_map[sample] = individual
            t_tmp = time.time()
            snpdb.import_samples(
                sample_reader=readers.VcfSampleReader(vcffname),
                map_name=experiment_id + size_id + 'VCF',
                id_map=id_map,
                report=False)
            t += time.time() - t_tmp
            # * PLink Files
            # Importing map file
            t_tmp = time.time()
            snpdb.import_map(map_reader=readers.PlinkMapReader(mfnames['PL']),
                             map_name=experiment_id + size_id + 'PL',
                             force_create_new=True,
                             force_use_existing=False,
                             report=False)
            t += time.time() - t_tmp
            # Importing sample file
            id_map = {}
            # Linking samples to individuals in the database
            if ifnames.get('PL') is not None:
                with open(ifnames['PL'], "r") as f:
                    for line in f:
                        (sample, individual) = line.split()
                        id_map[sample] = individual
            t_tmp = time.time()
            snpdb.import_samples(sample_reader=readers.PlinkSampleReader(
                pfnames['PL']),
                                 map_name=experiment_id + size_id + 'PL',
                                 id_map=id_map,
                                 report=False)
            t += time.time() - t_tmp
            # * FastQ Files
            # Importing sequencing file
            with open(fqfname, "rb") as fqf:
                t_tmp = time.time()
                snpdb.insert_file(file=fqf, individual_id="IND1")
                t += time.time() - t_tmp
            # * Media Files
            # Importing media file
            with open(imfname, "rb") as imf:
                t_tmp = time.time()
                snpdb.insert_file(file=imf, individual_id="IND1")
                t += time.time() - t_tmp

            # Validating Statistics
            snpdb._db.command("validate",
                              snpdb._config["MAPS_COLL"],
                              full=True)
            snpdb._db.command("validate",
                              snpdb._config["MAPSNPS_COLL"],
                              full=True)
            snpdb._db.command("validate",
                              snpdb._config["SNPS_COLL"],
                              full=True)
            snpdb._db.command("validate",
                              snpdb._config["SAMPLES_COLL"],
                              full=True)
            snpdb._db.command("validate",
                              snpdb._config["SNPBLOCKS_COLL"],
                              full=True)
            snpdb._db.command("validate",
                              snpdb._config["INDIVIDUALS_COLL"],
                              full=True)
            snpdb._db.command("validate", "fs.chunks", full=True)
            snpdb._db.command("validate", "fs.files", full=True)

            # Appending experiment results
            fsize: int = sum(
                [os.stat(mfnames[i]).st_size for i in mfnames.keys()] +
                [os.stat(pfnames[i]).st_size for i in pfnames.keys()] +
                [os.stat(ifnames[i]).st_size for i in ifnames.keys()] + [
                    os.stat(frfname).st_size,
                    os.stat(vcffname).st_size,
                    os.stat(fqfname).st_size,
                    os.stat(imfname).st_size
                ])
            dbsize: float = snpdb._db.command("dbStats")["storageSize"]

            results[experiment_id][size_id][count_id]['fsize'].append(fsize)
            results[experiment_id][size_id][count_id]['dbsize'].append(dbsize)
            results[experiment_id][size_id][count_id]['time'].append(t)

# %% [markdown]
# # Armazenando resultados

# %%
# Saving results to a JSON file
jfname = data_dir + 'exp2_1results.json'
with open(jfname, 'w') as jf:
    json.dump(results, jf)
