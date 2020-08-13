# %% [markdown]
# # Experimentos SNPDB (gráficos)

# %% [markdown]
# ## Imports e inicialização

# %%
import itertools
import json
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

data_dir: str = './data/'  # data output directory
graph_dir: str = './graphs/'  # graph output directory

# %% [markdown]
# ## Gráficos para análise de resultados
# - Para cada tipo de arquivo, foram medidos os tempos de execução para
# inserção, levando em consideração arquivos com 100 mil e 1 milhão de
# marcadores (SNPs).
# - Os arquivos (com exceção do FastQ) foram gerados aleatoriamente.
# - Para cada tipo de arquivo, o experimento foi executado 10 vezes,
# visando aumentar a precisão dos resultados, diminuindo variações
# provenientes de ruídos.
# - Para cada indivíduo, foi gerada uma única mostra para as quantidades de
# marcadores analisadas. Verificar se é necessário um número maior de
# amostras por arquivo de mapa.
# - Nos gráficos abaixo, cada tipo de arquivo / tamanho do mapa está
# devidamente etiquetado. Ainda são possíveis melhorias na apresentação dos
# resultados.
# - O MongoDB apresenta dois métodos de compressão de dados: *snappy* (padrão)
# e *zlib*. O *snappy* é mais eficiente em custo de processamento, enquanto
# o *zlib* apresenta um maior fator de compressão. Os gráficos foram dividos
# em duas seções, para comparação dos dois métodos de compressão.

# %% [markdown]
# ### Abrindo e processando arquivo de resultados (snappy)

# %%
with open(graph_dir + 'exp1results_snappy.json', 'r') as f:
    results = json.load(f)

names: dict = {
    '1.A': 'Z125',
    '1.B': 'FR',
    '1.C': 'VCF',
    '1.D': 'PLINK',
    '1.E': 'FastQ',
    '1.F': 'Media'
}
# sizes: dict = {'100k': 100000, '1m': 1000000}

snps_ids: list = ['1.A', '1.B', '1.C', '1.D']
snps_sizes: list = ['100k', '1m']
files_ids: list = ['1.E', '1.F']
col_names: list = []
times: list = []
fsizes: list = []
dbsizes: list = []
# file_type: list = []

# Processing data into a DataFrame
for i in snps_ids:
    for j in snps_sizes:
        col_names.append(names[i] + '_' + j)

        # n_snps
        times.append([
            x + y for x, y in zip(results[i][j]['map']['time'], results[i][j]
                                  ['sample']['time'])
        ])
        fsizes.append(results[i][j]['fsize'])
        dbsizes.append(results[i][j]['dbsize'])

for i in files_ids:
    col_names.append(names[i])
    times.append(results[i]['time'])
    fsizes.append(results[i]['fsize'])
    dbsizes.append(results[i]['dbsize'])

col_names = ['times_' + i
             for i in col_names] + ['fsizes_' + i for i in col_names
                                    ] + ['dbsizes_' + i for i in col_names]

df = pd.DataFrame(data=np.array(times + fsizes + dbsizes).T,
                  columns=col_names,
                  dtype=np.float,
                  copy=True)

del snps_ids, snps_sizes, files_ids, times, fsizes, dbsizes

# %% [markdown]
# ### Tempos de execução

# %%
sns.set()
plt.figure(figsize=(10, 5))
g = sns.barplot(data=df.iloc[:, :10], saturation=0.5, orient='h', ci=None)
# g.set(yscale="log")
for idx, i in enumerate(df.iloc[:, :10].mean(axis=0)):
    g.text(i + 10, idx, round(i, 3), color='black', ha="center")
plt.title("Comparação de desempenho de inserções de arquivos (N=10)")
plt.ylabel("Tipo de arquivo / tamanho")
plt.xlabel("Tempo médio de execução (s)")
plt.savefig(graph_dir + './times_snappy.png')

# %% [markdown]
# ### Comparação de tamanhos de arquivos antes e depois das inserções

# %%
cols = list(itertools.chain.from_iterable(zip(range(10, 18), range(20, 28))))
cols += [19, 29]
df2 = df.iloc[:, cols].apply(lambda x: x / 1024**2)

sns.set()
plt.figure(figsize=(10, 5))
g = sns.barplot(data=df2, saturation=0.5, orient='h', ci=None)
# g.set(yscale="log")
for idx, i in enumerate(df2.mean(axis=0)):
    g.text(i + 3, idx, round(i, 2), color='black', ha="center")
plt.title("Comparação de armazenamento de arquivos (N=10)")
plt.ylabel("Tipo de arquivo antes e depois da inserção")
plt.xlabel("Tamanho médio dos arquivos de entrada (MB)")
plt.savefig(graph_dir + './files_snappy.png')

# %%
cols = [18, 28]
df2 = df.iloc[:, cols].apply(lambda x: x / 1024**2)

sns.set()
plt.figure(figsize=(10, 5))
g = sns.barplot(data=df2, saturation=0.5, orient='h', ci=None)
# g.set(yscale="log")
for idx, i in enumerate(df2.mean(axis=0)):
    g.text(i + 2000, idx, round(i, 2), color='black', ha="center")
plt.title("Comparação de armazenamento de arquivos (N=10)")
plt.ylabel("Tipo de arquivo antes e depois da inserção")
plt.xlabel("Tamanho médio dos arquivos de entrada (MB)")
plt.savefig(graph_dir + './fastq_snappy.png')

# %% [markdown]
# ## Abrindo e processando arquivo de resultados (zlib)

# %%
with open(graph_dir + 'exp1results_zlib.json', 'r') as f:
    results = json.load(f)

names = {
    '1.A': 'Z125',
    '1.B': 'FR',
    '1.C': 'VCF',
    '1.D': 'PLINK',
    '1.E': 'FastQ',
    '1.F': 'Media'
}
# sizes: dict = {'100k': 100000, '1m': 1000000}

snps_ids = ['1.A', '1.B', '1.C', '1.D']
snps_sizes = ['100k', '1m']
files_ids = ['1.E', '1.F']
col_names = []
times = []
fsizes = []
dbsizes = []
# file_type: list = []

# Processing data into a DataFrame
for i in snps_ids:
    for j in snps_sizes:
        col_names.append(names[i] + '_' + j)

        # n_snps
        times.append([
            x + y for x, y in zip(results[i][j]['map']['time'], results[i][j]
                                  ['sample']['time'])
        ])
        fsizes.append(results[i][j]['fsize'])
        dbsizes.append(results[i][j]['dbsize'])

for i in files_ids:
    col_names.append(names[i])
    times.append(results[i]['time'])
    fsizes.append(results[i]['fsize'])
    dbsizes.append(results[i]['dbsize'])

col_names = ['times_' + i
             for i in col_names] + ['fsizes_' + i for i in col_names
                                    ] + ['dbsizes_' + i for i in col_names]

df = pd.DataFrame(data=np.array(times + fsizes + dbsizes).T,
                  columns=col_names,
                  dtype=np.float,
                  copy=True)

del snps_ids, snps_sizes, files_ids, times, fsizes, dbsizes

# %% [markdown]
# ### Tempos de execução

# %%
sns.set()
plt.figure(figsize=(10, 5))
g = sns.barplot(data=df.iloc[:, :10], saturation=0.5, orient='h', ci=None)
# g.set(yscale="log")
for idx, i in enumerate(df.iloc[:, :10].mean(axis=0)):
    g.text(i + 10, idx, round(i, 3), color='black', ha="center")
plt.title("Comparação de desempenho de inserções de arquivos (N=10)")
plt.ylabel("Tipo de arquivo / tamanho")
plt.xlabel("Tempo médio de execução (s)")
plt.savefig(graph_dir + './times_zlib.png')

# %% [markdown]
# #### Comparação de tamanhos de arquivos antes e depois das inserções

# %%
cols = list(itertools.chain.from_iterable(zip(range(10, 18), range(20, 28))))
cols += [19, 29]
df2 = df.iloc[:, cols].apply(lambda x: x / 1024**2)

sns.set()
plt.figure(figsize=(10, 5))
g = sns.barplot(data=df2, saturation=0.5, orient='h', ci=None)
# g.set(yscale="log")
for idx, i in enumerate(df2.mean(axis=0)):
    g.text(i + 3, idx, round(i, 2), color='black', ha="center")
plt.title("Comparação de armazenamento de arquivos (N=10)")
plt.ylabel("Tipo de arquivo antes e depois da inserção")
plt.xlabel("Tamanho médio dos arquivos de entrada (MB)")
plt.savefig(graph_dir + './files_zlib.png')

# %%
cols = [18, 28]
df2 = df.iloc[:, cols].apply(lambda x: x / 1024**2)

sns.set()
plt.figure(figsize=(10, 5))
g = sns.barplot(data=df2, saturation=0.5, orient='h', ci=None)
# g.set(yscale="log")
for idx, i in enumerate(df2.mean(axis=0)):
    g.text(i + 2000, idx, round(i, 2), color='black', ha="center")
plt.title("Comparação de armazenamento de arquivos (N=10)")
plt.ylabel("Tipo de arquivo antes e depois da inserção")
plt.xlabel("Tamanho médio dos arquivos de entrada (MB)")
plt.savefig(graph_dir + './fastq_zlib.png')

# %% [markdown]
# ## Comentários
# - O tempo de execução resultante aparenta condizer com o apresentado na
# monografia do Rodrigo.
# - Quanto ao espaço necessário de armazenamento, o método de compressão zlib
# apresentou um fator de compressão relativamente maior, quando comparado ao
# snappy. Contudo, em alguns casos, o espaço necessário para armazenamento
# pós-importação ainda é maior que o tamanho original dos arquivos de dados.

# %%
