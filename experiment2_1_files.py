# %% [markdown]
#  # Experimentos SNPDB

# %% [markdown]
#  ## Imports e inicialização

# %%
import multiprocessing
import numpy as np
import os
from PIL import Image
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
#  -> None:
#     """Reset MongoDB database for experiments."""
#     snpdb._client.drop_database(snpdb._config["DB_NAME"])
#     subprocess.run("""mongo --eval "load('./mongo_setup.js')" """,
#                    stdout=subprocess.DEVNULL,
#                    stderr=subprocess.STDOUT,
#                    shell=True)
#     importlib.reload(snpdb)


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
# ## Geração de arquivos de dados

# %%
N: int = 1  # Performing experiments with N loops
# nspns_list: list = [100000, 1000000]
nsnps_list: list = [100000]
nsamples_list: list = [10, 100, 1000, 10000, 100000, 1000000]

# multiprocessing variables
jobs: list = []
p = multiprocessing.Process()

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

for nsnps in nsnps_list:
    size_id: str = {100000: '100k', 1000000: '1m'}[nsnps]
    for nsamples in nsamples_list:
        count_id: str = {
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
        t: float = 0  # timing variable

        # Filenames
        mfname: str = data_dir + 'out_' + size_id + '_' + count_id + '.0125map'
        pfname: str = data_dir + 'out_' + size_id + '_' + count_id + '.0125ped'
        ifname: str = data_dir + 'out_' + size_id + '_' + count_id + '.0125ids'

        for i in range(N):
            print("Generating input files...")
            print("i: " + str(i))
            # * Generating input files
            # Map file
            p = multiprocessing.Process(target=generate_random_file,
                                        kwargs={
                                            'filename': mfname,
                                            'file_type': '0125_map',
                                            'verbose': True,
                                            'n': nsnps
                                        })
            jobs.append(p)
            p.start()
            #  start_from_id=start_map)
            # Samples file
            p = multiprocessing.Process(target=generate_random_file,
                                        kwargs={
                                            'filename': pfname,
                                            'file_type': '0125_samples',
                                            'verbose': True,
                                            'n': nsamples,
                                            'map_size': nsnps
                                        })
            jobs.append(p)
            p.start()
            #  start_from_id=start_sample)
            # Id map file
            p = multiprocessing.Process(target=generate_random_file,
                                        kwargs={
                                            'filename': ifname,
                                            'file_type': 'id_mapping',
                                            'verbose': True,
                                            'n': nsamples
                                        })
            jobs.append(p)
            p.start()
            #  first_sample_id=start_sample)

# %% [markdown]
# ### 2.1.B - Arquivo VCF

# %%
experiment_id = '2.1.B'

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
        t = 0.0  # timing variable

        # Filenames
        # ? VCF -> one file for both map and samples
        vcffname: str = data_dir + 'out_' + size_id + '_' + count_id + '.vcf'
        ifname = data_dir + 'out_' + size_id + '_' + count_id + '.vcfids'

        # * Performing experiment N times and storing results
        for i in range(N):
            print("i: " + str(i))
            print("Generating input files...")
            # * Generating input files
            p = multiprocessing.Process(target=generate_random_file,
                                        kwargs={
                                            'filename': vcffname,
                                            'file_type': 'vcf',
                                            'verbose': True,
                                            'n': nsamples,
                                            'map_size': nsnps
                                        })
            jobs.append(p)
            p.start()
            #  start_snps_from_id=start_map,
            #  start_samples_from_id=start_sample)
            # Id map file
            p = multiprocessing.Process(target=generate_random_file,
                                        kwargs={
                                            'filename': ifname,
                                            'file_type': 'id_mapping',
                                            'verbose': True,
                                            'n': nsamples
                                        })
            jobs.append(p)
            p.start()
            #  first_sample_id=start_sample)

# %% [markdown]
# ### 2.1.C - Todos os tipos de arquivos

# %%
experiment_id = '2.1.C'

im_res: tuple = (800, 600)
fqfname: str = fastq_dir_1 + 'SH.71992.AP.01.1.fastq'
imfname: str = data_dir + 'out_image.jpg'

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
        t = 0.0  # timing variable
        t_tmp: float = 0.0  # temporary timing variable

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

        # * Performing experiment N times and storing results
        for i in range(N):
            print("i: " + str(i))
            # * Generating input files
            # generate_random_file(filename=mfnames['Z125'],
            #                      file_type='0125_map',
            #                      verbose=True,
            #                      n=nsnps)
            # #  start_from_id=start_map)
            # generate_random_file(filename=pfnames['Z125'],
            #                      file_type='0125_samples',
            #                      verbose=True,
            #                      n=nsamples,
            #                      map_size=nsnps)
            # #  start_from_id=start_sample)
            # generate_random_file(filename=ifnames['Z125'],
            #                      file_type='id_mapping',
            #                      verbose=True,
            #                      n=nsamples)
            # #  first_sample_id=start_sample)
            p = multiprocessing.Process(target=generate_random_file,
                                        kwargs={
                                            'filename': frfname,
                                            'file_type': 'final_report',
                                            'verbose': True,
                                            'n': nsamples,
                                            'map_size': nsnps
                                        })
            jobs.append(p)
            p.start()
            #  start_snps_from_id=start_map,
            #  start_samples_from_id=start_sample)
            p = multiprocessing.Process(target=generate_random_file,
                                        kwargs={
                                            'filename': ifnames['FR'],
                                            'file_type': 'id_mapping',
                                            'verbose': True,
                                            'n': nsamples
                                        })
            jobs.append(p)
            p.start()
            #  first_sample_id=start_sample)
            # generate_random_file(filename=vcffname,
            #                      file_type='vcf',
            #                      verbose=True,
            #                      n=nsamples,
            #                      map_size=nsnps)
            # #  start_snps_from_id=start_map,
            # #  start_samples_from_id=start_sample)
            # # Id map file
            # generate_random_file(filename=ifnames['VCF'],
            #                      file_type='id_mapping',
            #                      verbose=True,
            #                      n=nsamples)
            # #  first_sample_id=start_sample)
            p = multiprocessing.Process(target=generate_random_file,
                                        kwargs={
                                            'filename': mfnames['PL'],
                                            'file_type': 'plink_map',
                                            'verbose': True,
                                            'n': nsnps
                                        })
            jobs.append(p)
            p.start()
            #  start_from_id=start_map)
            p = multiprocessing.Process(target=generate_random_file,
                                        kwargs={
                                            'filename': pfnames['PL'],
                                            'file_type': 'plink_samples',
                                            'verbose': True,
                                            'n': nsamples,
                                            'map_size': nsnps
                                        })
            jobs.append(p)
            p.start()
            #  start_from_id=start_sample)
            p = multiprocessing.Process(target=generate_random_file,
                                        kwargs={
                                            'filename': ifnames['PL'],
                                            'file_type': 'id_mapping',
                                            'verbose': True,
                                            'n': nsamples
                                        })
            jobs.append(p)
            p.start()
            #  first_sample_id=start_sample)
            im_arr = np.random.rand(im_res[0], im_res[1], 3) * 255
            im_out = Image.fromarray(im_arr.astype('uint8')).convert('RGB')
            im_out.save(imfname)

# %% [markdown]
# # Synchronizing processes

# %%
for p in jobs:
    p.join()
