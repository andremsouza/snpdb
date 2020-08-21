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
# ## 2.2 - Sumarização
# Dado um indivíduo, analisar o tempo de consulta para sumarização de
# informação relacionadas a esteindivíduo (a definir):
# - Quantidade de amostras de cada "tipo de arquivo"(0125, FR, VCF, PLINK),
# descrição;
# - Quantidade de marcadores (SNPs) de cada "tipo de arquivo";
# - Quantidade de arquivos de sequenciamento (FastQ), descrição;
# - Quantidade de arquivos de mídia, descrição.
#
# Já foi implementada uma função de sumarização de dados de um indivíduo,
# contendo:
# - Arquivos associados a este indivíduo;
# - Mapas e amostras associados a este indivíduo, separados por formato de
# arquivo de dados.

# %%


def summarize(individual: dict):
    pass


# %% [markdown]
# # Armazenando resultados

# %%
# Saving results to a JSON file
jfname = data_dir + 'exp2_2results.json'
with open(jfname, 'w') as jf:
    json.dump(results, jf)
