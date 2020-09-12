# %% [markdown]
# # Experimentos SNPDB

# %% [markdown]
# ## Imports e inicialização

# %%
from experiment_config import exps, nsnps_ids, nsamples_ids
import importlib
import json
import numpy as np
import os
from PIL import Image
from pprint import pprint
import snpdb
import subprocess
import testing.random_file_generator as rfgen
import time
from typing import Callable, Tuple
# import writers

# Constants and global variables
data_dir: str = './data/'  # Data output directory
fastq_dir_1: str = '../1_FASTq_VCF/1_SRR10752699_MISSOURI/'  # data directory
fastq_dir_2: str = '../1_FASTq_VCF/2_VCF_bos_taurus/'  # data directory
graph_dir: str = './graphs/'  # Graph output directory
results_fname: str = data_dir + 'experiment_results.json'
results: dict = {}
try:
    with open(results_fname, 'r') as f:
        results = json.load(f)
except Exception:
    results = {}

# Verifying existing files
print('Directory:', data_dir)
pprint(os.listdir(data_dir))
print('Directory:', fastq_dir_1)
pprint(os.listdir(fastq_dir_1))
print('Directory:', fastq_dir_2)
pprint(os.listdir(fastq_dir_2))

# %% [markdown]
# ## Definição de funções


# %%
def reset_db(compression_method: str = 'snappy') -> None:
    """Reset MongoDB database for experiments.

    Args:
        compression_method (str, optional): Compression method for database
            after reset operation. May be either 'snappy' or 'zlib'.
            Defaults to 'snappy'.
    """
    subprocess.run(['/home/rodrigo/snpdb/reset_db.sh', compression_method],
                   stdout=subprocess.DEVNULL,
                   stderr=subprocess.STDOUT,
                   check=True)
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
        verbose (bool, optional): If True, prints timing and size
            information about the generated file. Defaults to False.
        **kwargs (Any): Keyword arguments for file generation function.
            'outfile' may be omitted.

    Returns:
        Tuple[float, float]: Returns the function generation time, and its
            size on disk
    """
    func: Callable = {
        '.0125map': rfgen.random_0125_map,
        '.0125ped': rfgen.random_0125_samples,
        '.fr': rfgen.random_final_report,
        '.vcf': rfgen.random_vcf,
        '.plmap': rfgen.random_plink_map,
        '.plped': rfgen.random_plink_samples,
        '.ids': rfgen.id_mapping
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


def execute_experiment_two_files(result: dict) -> None:
    """Executes experiment for file types with two files (map, sample)
    Warning: this function uses variables defined outside its scope

    Args:
        result (dict): Dictionary for experiment. Values will be assigned
                       to this dictionary "by reference".
    """
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
        reset_db()
        print("Database reset operation successful.")
        print("Generating input files...")
        t_map: float = 0.0
        t_sample: float = 0.0
        map_size: float = 0.0
        sample_size: float = 0.0
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
        t_tmp: float = time.time()
        snpdb.import_map(
            map_reader=exps[experiment_id]['readers']['map'](mfname),
            map_name=experiment_id + '_' + nsnps_id + '_' + nsamples_id,
            force_create_new=True,
            force_use_existing=False,
            report=False)
        t_tmp = time.time() - t_tmp
        t_map += t_tmp
        # Validating statistics
        snpdb._db.command("validate", snpdb._config["MAPS_COLL"], full=True)
        snpdb._db.command("validate", snpdb._config["MAPSNPS_COLL"], full=True)
        snpdb._db.command("validate", snpdb._config["SNPS_COLL"], full=True)
        # map_size = _MAPS + _MAPSNPS + _SNPS
        map_size = snpdb._db.command(
            "collstats",
            snpdb._config["MAPS_COLL"])['storageSize'] + snpdb._db.command(
                "collstats", snpdb.
                _config["MAPSNPS_COLL"])['storageSize'] + snpdb._db.command(
                    "collstats", snpdb._config["SNPS_COLL"])['storageSize']
        print("Imported map file\tTime: " + str(round(t_map, 3)) + "s\tSize:" +
              str(round(map_size / 1024**2, 2)) + "MB")
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
            t_tmp = time.time()
            snpdb.import_samples(
                sample_reader=exps[experiment_id]['readers']['ped'](pfname),
                map_name=experiment_id + '_' + nsnps_id + '_' + nsamples_id,
                id_map=id_map,
                report=False)
            t_tmp = time.time() - t_tmp
            t_sample += t_tmp
        # Validating Statistics
        snpdb._db.command("validate", snpdb._config["SAMPLES_COLL"], full=True)
        snpdb._db.command("validate",
                          snpdb._config["SNPBLOCKS_COLL"],
                          full=True)
        snpdb._db.command("validate",
                          snpdb._config["INDIVIDUALS_COLL"],
                          full=True)
        # sample_size = _SAMPLES + _SNPBLOCKS + _INDS
        sample_size = snpdb._db.command(
            "collstats",
            snpdb._config["SAMPLES_COLL"])['storageSize'] + snpdb._db.command(
                "collstats", snpdb._config["SNPBLOCKS_COLL"]
            )['storageSize'] + snpdb._db.command(
                "collstats", snpdb._config["INDIVIDUALS_COLL"])['storageSize']
        print("Imported samples file\tTime: " + str(round(t_sample, 3)) +
              "s\tSize:" + str(round(sample_size / 1024**2, 2)) + "MB")

        # Appending generated file sizes
        result['fsize'].append(
            float(os.stat(mfname).st_size) + float(os.stat(pfname).st_size))
        # Appending stored document sizes from MongoDB
        result['dbsize'].append(map_size + sample_size)
        # Appending insertion times
        result['time'].append(t_map + t_sample)

        # Executing additional steps

        # 2.2 Sumarização
        ind = np.random.choice(snpdb.find_individuals())
        t_tmp = time.time()
        summary = snpdb.summarize(ind)
        t_tmp = time.time() - t_tmp
        result['summarize'].append({
            'individual': ind,
            'summary': summary,
            'time': t_tmp
        })

        # TODO: 2.3 Exportação de sumarização para formatos originais

        # 2.4 Busca de indivíduos, dada uma lista de SNPs
        snp = np.random.choice(snpdb.find_snp())
        result['individuals_of_snps'].append({})
        result['individuals_of_snps'][-1]
        t_tmp = time.time()
        try:
            inds = snpdb.find_individuals_of_snps(id=snp['i'], )
            result['individuals_of_snps'][-1]['snp'] = snp
            result['individuals_of_snps'][-1]['individuals'] = inds
        except Exception as e:
            print("Warning: couldn't retrieve individuals from database", e)
        t_tmp = time.time() - t_tmp
        result['individuals_of_snps'][-1]['snp'] = t_tmp

        # TODO: 2.5 Exportação de indivíduos para formatos originais

        # 2.6 Remoção de todos os dados de um indivíduo
        ind = np.random.choice(snpdb.find_individuals())
        t_tmp = time.time()
        delete_results = snpdb.delete_individuals(id=ind['_id'])
        t_tmp = time.time() - t_tmp
        result['delete_individual'].append({
            'individual':
            ind,
            'deleted_count': [i.deleted_count for i in delete_results],
            'time':
            t_tmp
        })

        # Writing partial results to file
        with open(results_fname, 'w') as f:
            json.dump(results,
                      f,
                      ensure_ascii=True,
                      check_circular=True,
                      allow_nan=True,
                      indent=1,
                      sort_keys=True)


def execute_experiment_one_file(result: dict) -> None:
    """Executes experiment for file types with one file
    Warning: this function uses variables defined outside its scope

    Args:
        result (dict): Dictionary for experiment. Values will be assigned
                       to this dictionary "by reference".
    """
    print("Starting Experiment " + experiment_id + " (" + nsnps_id +
          " SNPs, " + nsamples_id + " individuals) with N = " + str(N) +
          "; Compression method: " + compression_method)

    # get filenames
    f_ext: dict = exps[experiment_id]['file_extensions']
    fname = data_dir + 'out_' + nsnps_id + '_' + nsamples_id + f_ext['ext']
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
        reset_db()
        print("Database reset operation successful.")
        print("Generating input files...")
        t_map = 0.0
        t_sample = 0.0
        map_size = 0.0
        sample_size = 0.0
        # * Generating input files
        # If less than 10000 samples, generate in one file
        # Else, generate blocks of up to 10000 samples
        n_blocks = int(np.ceil(nsamples / 10000))
        remaining_samples = nsamples
        start_sample = 1
        imported_map = False
        #  start_from_id=start_map)
        for i in range(n_blocks):
            print("Block: " + str(i))
            nsamples_block = int(np.minimum(remaining_samples, 10000.0))
            # Map/Samples file
            generate_random_file(filename=fname,
                                 file_type=f_ext['ext'],
                                 verbose=True,
                                 n=nsamples_block,
                                 map_size=nsnps,
                                 start_samples_from_id=start_sample)
            # Id map file
            generate_random_file(filename=ifname,
                                 file_type='.ids',
                                 verbose=True,
                                 n=nsamples_block,
                                 first_sample_id=start_sample)
            start_sample += nsamples_block
            remaining_samples -= nsamples_block
            # * Inserting input files into db
            print("Inserting input files into database...")
            if not imported_map:
                imported_map = True
                # Importing map file
                t_tmp = time.time()
                snpdb.import_map(
                    map_reader=exps[experiment_id]['readers']['map'](fname),
                    map_name=experiment_id + '_' + nsnps_id + '_' +
                    nsamples_id,
                    force_create_new=True,
                    force_use_existing=False,
                    report=False)
                t_tmp = time.time() - t_tmp
                t_map += t_tmp
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
                    "collstats", snpdb._config["MAPS_COLL"]
                )['storageSize'] + snpdb._db.command(
                    "collstats", snpdb._config["MAPSNPS_COLL"]
                )['storageSize'] + snpdb._db.command(
                    "collstats", snpdb._config["SNPS_COLL"])['storageSize']
                print("Imported map file\tTime: " + str(round(t_map, 3)) +
                      "s\tSize:" + str(round(map_size / 1024**2, 2)) + "MB")
            # Importing sample file
            id_map: dict = {}
            # Linking samples to individuals in the database
            if ifname is not None:
                with open(ifname, "r") as f:
                    for line in f:
                        (sample, individual) = line.split()
                        id_map[sample] = individual
            t_tmp = time.time()
            snpdb.import_samples(
                sample_reader=exps[experiment_id]['readers']['ped'](fname),
                map_name=experiment_id + '_' + nsnps_id + '_' + nsamples_id,
                id_map=id_map,
                report=False)
            t_tmp = time.time() - t_tmp
            t_sample += t_tmp
        # Validating Statistics
        snpdb._db.command("validate", snpdb._config["SAMPLES_COLL"], full=True)
        snpdb._db.command("validate",
                          snpdb._config["SNPBLOCKS_COLL"],
                          full=True)
        snpdb._db.command("validate",
                          snpdb._config["INDIVIDUALS_COLL"],
                          full=True)
        # sample_size = _SAMPLES + _SNPBLOCKS + _INDS
        sample_size = snpdb._db.command(
            "collstats",
            snpdb._config["SAMPLES_COLL"])['storageSize'] + snpdb._db.command(
                "collstats", snpdb._config["SNPBLOCKS_COLL"]
            )['storageSize'] + snpdb._db.command(
                "collstats", snpdb._config["INDIVIDUALS_COLL"])['storageSize']
        print("Imported samples file\tTime: " + str(round(t_sample, 3)) +
              "s\tSize:" + str(round(sample_size / 1024**2, 2)) + "MB")

        # Appending generated file sizes
        result['fsize'].append(float(os.stat(fname).st_size))
        # Appending stored document sizes from MongoDB
        result['dbsize'].append(map_size + sample_size)
        # Appending insertion times
        result['time'].append(t_map + t_sample)

        # Executing additional steps

        # 2.2 Sumarização
        ind = np.random.choice(snpdb.find_individuals())
        t_tmp = time.time()
        summary = snpdb.summarize(ind)
        t_tmp = time.time() - t_tmp
        result['summarize'].append({
            'individual': ind,
            'summary': summary,
            'time': t_tmp
        })

        # TODO: 2.3 Exportação de sumarização para formatos originais

        # 2.4 Busca de indivíduos, dada uma lista de SNPs
        snp = np.random.choice(snpdb.find_snp())
        t_tmp = time.time()
        inds = snpdb.find_individuals_of_snps(id=snp['i'])
        t_tmp = time.time() - t_tmp
        result['individuals_of_snps'].append({
            'snp': snp,
            'individuals': inds,
            'time': t_tmp
        })

        # TODO: 2.5 Exportação de indivíduos para formatos originais

        # 2.6 Remoção de todos os dados de um indivíduo
        ind = np.random.choice(snpdb.find_individuals())
        t_tmp = time.time()
        delete_results = snpdb.delete_individuals(id=ind['_id'])
        t_tmp = time.time() - t_tmp
        result['delete_individual'].append({
            'individual':
            ind,
            'deleted_count': [i.deleted_count for i in delete_results],
            'time':
            t_tmp
        })

        # Writing partial results to file
        with open(results_fname, 'w') as f:
            json.dump(results,
                      f,
                      ensure_ascii=True,
                      check_circular=True,
                      allow_nan=True,
                      indent=1,
                      sort_keys=True)


def execute_experiment_bin_file(result: dict) -> None:
    """Executes experiment for file types inserted in binary mode
    Warning: this function uses variables defined outside its scope

    Args:
        result (dict): Dictionary for experiment. Values will be assigned
                       to this dictionary "by reference".
    """
    print("Starting Experiment " + experiment_id + " with N = " + str(N) +
          "; Compression method: " + compression_method)

    # get filenames
    f_ext: dict = exps[experiment_id]['file_extensions']
    fname: str = ''
    if file_type == 'FastQ':
        fname = fastq_dir_1 + 'SH.71992.AP.01.1.fastq'
    elif file_type == 'Media':
        fname = data_dir + 'out' + f_ext['ext']

    # Setting up result dictionary
    result['fsize'] = []  # file size
    result['dbsize'] = []  # doc size in db
    result['time'] = []  # insertion time

    # * Performing experiment N times and storing results
    for i in range(N):
        print("i: " + str(i))
        print("Resetting database...")
        reset_db()
        print("Database reset operation successful.")
        if file_type == 'FastQ':
            print(fname + "\tSize: " +
                  str(round(os.stat(fname).st_size / (1024**2), 2)) + "MB")
            # * Inserting input files into db
            print("Inserting input files into database...")
            # Importing map file
            # Importing sequencing file
            with open(fname, "rb") as fqf:
                t = time.time()
                snpdb.insert_file(file=fqf, individual_id=0)
                t = time.time() - t
            # Validating Statistics
            snpdb._db.command("validate", "fs.chunks", full=True)
            snpdb._db.command("validate", "fs.files", full=True)
            file_size = snpdb._db.command("collstats",
                                          "fs.chunks")["storageSize"]
            # file_size = snpdb._db["fs.files"].find_one()["length"]
            print("Imported sequencing file\tTime: " + str(round(t, 3)) +
                  "s\tSize:" + str(round(file_size / 1024**2, 2)) + "MB")
        elif file_type == 'Media':
            im_res = (800, 600)  # Image resolution
            print("Generating input file...")
            # * Generating input file
            t = time.time()
            im_arr = np.random.rand(im_res[0], im_res[1], 3) * 255
            im_out = Image.fromarray(im_arr.astype('uint8')).convert('RGB')
            im_out.save(fname)
            t = time.time() - t
            print("Generated media file\tTime: " + str(round(t, 3)) +
                  "s\tSize: " +
                  str(round(os.stat(fname).st_size / (1024**2), 2)) + "MB")
            # * Inserting input files into db
            # Importing media file
            with open(fname, "rb") as imf:
                t = time.time()
                snpdb.insert_file(file=imf, individual_id=0)
                t = time.time() - t
            # Validating Statistics
            snpdb._db.command("validate", "fs.chunks", full=True)
            snpdb._db.command("validate", "fs.files", full=True)
            file_size = snpdb._db.command("collstats",
                                          "fs.chunks")["storageSize"]
            # file_size = snpdb._db["fs.files"].find_one()["length"]
            print("Imported media file\tTime: " + str(round(t, 3)) +
                  "s\tSize:" + str(round(file_size / 1024**2, 2)) + "MB")

        # Appending generated file sizes
        result['fsize'].append(float(os.stat(fname).st_size))
        # Appending stored document sizes from MongoDB
        result['dbsize'].append(file_size)
        # Appending insertion times
        result['time'].append(t)

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
#  ## 1 - Tempo de importação e espaço ocupado a partir da base zerada
#  A partir da base zerada, importar um arquivo de amostra e mapa de um animal
#
#  Considerando cada caso abaixo:
#  - A: 1 arquivo 0125 (incluindo Mapa)
#  - B: 1 arquivo PLINK (incluindo Mapa)
#  - C: 1 arquivo Final Report (FR) (incluindo Mapa)
#  - D: 1 arquivo VCF
#  - E: 1 arquivo FastQ (sequenciamento)
#  - F: arquivo de mídia
#  - G: todos os arquivos de A a F

# %% [markdown]
# # 2 - Análise de importação e consultas para as quantidades de indivíduos:
# - 10
# - 100
# - 1000
# - 10000
# - 100000
# - 1000000

# %% [markdown]
# ## 2.1 - Análise de performance de casos específicos
# Analisando performance considerando os casos com menor e maior tempo de
# importação da seção 1, assim como o caso G (tempo de importação para todos
# os tipos de arquivos).

# Tipos de arquivos analisados, segundo resultados do experimento 1:
# - 2.1.A - Arquivo 0125 {'100k': 5.727, '1m': 60.724}
# - 2.1.B - Arquivo VCF {'100k': 9.894, '1m': 72.766}
# - 2.1.C - Todos os tipos de arquivos, excluindo arquivos FastQ e de mídia

# %% [markdown]
# ## 1.A - Arquivo 0125
# ## 1.B - Arquivo PLINK
# ## 2.A - Arquivo 0125

# %%
# ? Execute experiment only if not exists in results
# ? Grouping experiments with similar formats
experiment_ids: set = {'1.A', '1.B', '2.A'}
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
                    execute_experiment_two_files(result)

# %% [markdown]
# ## 1.C - Arquivo FR
# ## 1.D - Arquivo VCF
# ## 2.B - Arquivo VCF

# %%
# ? Execute experiment only if not exists in results
# ? Grouping experiments with similar formats
experiment_ids = {'1.C', '1.D', '2.B'}
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
                    execute_experiment_one_file(result)

# %% [markdown]
# ## 1.E - Arquivo FastQ
# ## 1.F - Arquivo de mídia

# %%
# ? Execute experiment only if not exists in results
# ? Grouping experiments with similar formats
experiment_ids = {'1.E', '1.F'}
for experiment_id in experiment_ids:
    if experiment_id.startswith('1'):
        N = 10  # Performing experiments with N loops
    else:
        N = 1  # Only 1 loop for experiments 2.*
    file_type = exps[experiment_id]['file_type']
    results.setdefault(experiment_id, {})
    for compression_method in exps[experiment_id]['compression_methods']:
        results[experiment_id].setdefault(compression_method, {})
        # ? Using pointers to results dictionary
        result = results[experiment_id][compression_method]
        # Only execute experiment if empty
        if len(result) == 0:
            execute_experiment_bin_file(result)

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
