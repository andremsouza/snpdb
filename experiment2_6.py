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
import testing.random_file_generator as rfgen
import time
from typing import Callable, Tuple

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
nspns_list: list = [100000, 1000000]
nsamples_list: list = [10, 100, 1000, 10000, 100000, 1000000]

# %% [markdown]
# ## 2.6 - Remoção de todos os dados armazenados de um indivíduo
# Deverão ser implementadas funções para remoção de indivíduos específicos.

# %%
# TODO: Deletion function for snpdb


def delete_individuals(id=None, tatoo=None, sample_map=None, sample_id=None):
    """Searches and deletes all data from individuals in the database.

    Returns a list with the status of the delete operations.

    Parameters
    ----------
    id=None             Match only the individual with given internal id.
    tatoo=None          Match only individuals which contain tatoo among
                        their alternate IDs.
    sample_map=None     Match only individuals that have data on the specified
                        map.
    sample_id=None      Match only individuals that have the specified sample
                        id under some map.
    """
    result: list = []
    individuals: list = snpdb.find_individuals(id, tatoo, sample_map,
                                               sample_id)
    if len(individuals):
        query_samples: dict = {
            "$or": [{
                snpdb._config["SNPBLOCKS_MAP_ATTR"]: map_name,
                snpdb._config["SNPBLOCKS_SAMPLE_ATTR"]: sample_id
            } for map_name, sample_id in [
                zip(s[snpdb._config["SAMPLES_MAP_ATTR"]], s[
                    snpdb._config["SAMPLES_ID_ATTR"]]) for s in [
                        ind[snpdb._config["INDIVIDUALS_SAMPLE_LIST_ATTR"]]
                        for ind in individuals
                    ]
            ]]
        }
        query_individuals: dict = {
            "$or": [{
                "_id": id
            } for id in individuals["_id"]]
        }
        result.append(snpdb._SNPBLOCKS.deleteMany(query_samples))
        result.append(snpdb._SAMPLES.deleteMany(query_samples))
        result.append(snpdb._INDIVIDUALS.deleteMany(query_individuals))
    return result


# %% [markdown]
# # Armazenando resultados

# %%
# Saving results to a JSON file
jfname = data_dir + 'exp2_6results.json'
with open(jfname, 'w') as jf:
    json.dump(results, jf)
