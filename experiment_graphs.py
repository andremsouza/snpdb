# %% [markdown]
# # Experimentos SNPDB (gráficos)

# %% [markdown]
# ## Imports e inicialização

# %%
from experiment_config import graph_dir, bin_file_types
from experiment_config import results_fname, exps, nsnps_ids, nsamples_ids
import json
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pprint import pp
import seaborn as sns
import snpdb

# Import results from output of experiment.py
results: dict = {}
try:
    with open(results_fname, 'r') as f:
        results = json.load(f)
except Exception:
    results = {}

# Storing data in a Pandas DataFrame
col_names: list = [
    'experiment_id', 'compression_method', 'file_type', 'nsnps', 'nsamples',
    'fsize', 'dbsize', 'time'
]
df: pd.DataFrame = pd.DataFrame(columns=col_names)

# %% [markdown]
# ### Abrindo e processando arquivos de resultados

# %%

# Reading experiment result files and building DataFrame
for experiment_id in exps:
    if experiment_id not in results:
        continue
    for compression_method in exps[experiment_id]['compression_methods']:
        for nsnps in exps[experiment_id]['nsnps_list']:
            for nsamples in exps[experiment_id]['nsamples_list']:
                if exps[experiment_id]['file_type'] in bin_file_types:
                    nsnps = np.nan
                    nsamples = np.nan
                    nsnps_id = None
                    nsamples_id = None
                else:
                    nsnps_id = nsnps_ids[nsnps]
                    nsamples_id = nsamples_ids[nsamples]

                # Getting result values lists from json
                fsizes: list = []
                dbsizes: list = []
                times: list = []
                if exps[experiment_id]['file_type'] in bin_file_types:
                    fsizes = results[experiment_id][compression_method][
                        'fsize']
                    dbsizes = results[experiment_id][compression_method][
                        'dbsize']
                    times = results[experiment_id][compression_method]['time']
                    delete_times = [0.0] * len(fsizes)
                else:
                    fsizes = results[experiment_id][compression_method][
                        nsnps_id][nsamples_id]['fsize']
                    dbsizes = results[experiment_id][compression_method][
                        nsnps_id][nsamples_id]['dbsize']
                    times = results[experiment_id][compression_method][
                        nsnps_id][nsamples_id]['time']
                # Mounting rows and appending to DataFrame
                rows = [{
                    'experiment_id': experiment_id,
                    'compression_method': compression_method,
                    'file_type': exps[experiment_id]['file_type'],
                    'nsnps': nsnps,
                    'nsamples': nsamples,
                    'fsize': float(i),
                    'dbsize': float(j),
                    'time': float(k),
                } for i, j, k in zip(fsizes, dbsizes, times)]
                df = df.append(rows, ignore_index=True, verify_integrity=True)

# %% [markdown]
# ## Experimento 1 - Tempo de importação/espaço ocupado a partir da base zerada
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
# ### Tempos de execução

# %%
df_melted = pd.melt(df[df['experiment_id'].str.startswith('1')],
                    id_vars=[
                        'experiment_id', 'compression_method', 'file_type',
                        'nsnps', 'nsamples'
                    ],
                    value_vars=['time'])
df_melted = df_melted.fillna(value='Binary file')
sns.set(style="whitegrid",
        palette=sns.color_palette("muted", n_colors=6, desat=1.0))
snsplot = sns.catplot(
    x='file_type',
    y='value',
    hue='nsnps',
    data=df_melted,
    col='compression_method',
    col_order=['snappy', 'zlib'],
    kind='bar',
    height=8,
    aspect=1,
)

# Labelling bars of graph
for idx, g in enumerate(snsplot.axes[0]):
    c_method = {0: 'snappy', 1: 'zlib'}[idx]
    df_tmp = df_melted[df_melted['compression_method'] == c_method]
    for i, file_type in enumerate(df_tmp['file_type'].unique()):
        for j, nsnps in enumerate(df_tmp['nsnps'].unique()):
            x = i - 0.4 + (j * (0.8 / 3.0))
            if (file_type in bin_file_types and nsnps
                    == 'Binary file') or (file_type not in bin_file_types
                                          and nsnps != 'Binary file'):
                y = df_tmp[df_tmp['file_type'] == file_type][
                    df_tmp['nsnps'] == nsnps]['value'].mean() * 1.05
                label = "{:.2f}".format(
                    df_tmp[df_tmp['file_type'] == file_type][
                        df_tmp['nsnps'] == nsnps]['value'].mean())
            else:
                y = 0.0
                label = ""
            g.text(x, y, label)

snsplot._legend.set_title("# de SNPs")
# snsplot.set(xscale='log')
snsplot.set(yscale='log')
snsplot = snsplot.set_axis_labels("Tipo de arquivo",
                                  "Tempo de execução log10 (s)")
# .set_titles(template="Método de compressão: {row_name} - {col_name}")
snsplot.savefig(graph_dir + 'experiment1_times.png')
plt.draw()

# %% [markdown]
# ### Comparação de tamanhos de arquivos antes/depois das inserções (100k SNPs)

# %%
df_melted = pd.melt(df[df['experiment_id'].str.startswith('1')],
                    id_vars=[
                        'experiment_id', 'compression_method', 'file_type',
                        'nsnps', 'nsamples'
                    ],
                    value_vars=['fsize', 'dbsize'])
df_melted = df_melted.fillna(
    value='Binary file')[~df_melted['nsnps'].isin([1000000.0])]
df_melted['value'] /= 1024**2
sns.set(style="whitegrid",
        palette=sns.color_palette("muted", n_colors=6, desat=1.0))
snsplot = sns.catplot(
    x='file_type',
    y='value',
    hue='variable',
    data=df_melted,
    col='compression_method',
    col_order=['snappy', 'zlib'],
    kind='bar',
    height=8,
    aspect=1,
)

# Labelling bars of graph
for idx, g in enumerate(snsplot.axes[0]):
    c_method = {0: 'snappy', 1: 'zlib'}[idx]
    df_tmp = df_melted[df_melted['compression_method'] == c_method]
    for i, file_type in enumerate(df_tmp['file_type'].unique()):
        for j, variable in enumerate(df_tmp['variable'].unique()):
            x = i - 0.35 + (j * (0.8 / 2.0))
            y = df_tmp[df_tmp['file_type'] == file_type][
                df_tmp['variable'] == variable]['value'].mean() * 1.1
            label = "{:.2f}".format(df_tmp[df_tmp['file_type'] == file_type][
                df_tmp['variable'] == variable]['value'].mean())
            g.text(x, y, label)

snsplot._legend.set_title("")
# snsplot.set(xscale='log')
snsplot.set(yscale='log')
snsplot = snsplot.set_axis_labels(
    "Tipo de arquivo", "Tamanho do arquivo log10 (MB)").set_titles(
        template=" 100k SNPs - {col_var} = {col_name}", )
snsplot.savefig(graph_dir + 'experiment1_sizes_100k.png')
plt.draw()

# %% [markdown]
# ### Comparação de tamanhos de arquivos antes/depois das inserções (1m SNPs)

# %%
df_melted = pd.melt(df[df['experiment_id'].str.startswith('1')],
                    id_vars=[
                        'experiment_id', 'compression_method', 'file_type',
                        'nsnps', 'nsamples'
                    ],
                    value_vars=['fsize', 'dbsize'])
df_melted = df_melted.fillna(
    value='Binary file')[~df_melted['nsnps'].isin([100000.0])]
df_melted['value'] /= 1024**2
sns.set(style="whitegrid",
        palette=sns.color_palette("muted", n_colors=6, desat=1.0))
snsplot = sns.catplot(
    x='file_type',
    y='value',
    hue='variable',
    data=df_melted,
    col='compression_method',
    col_order=['snappy', 'zlib'],
    kind='bar',
    height=8,
    aspect=1,
)

# Labelling bars of graph
for idx, g in enumerate(snsplot.axes[0]):
    c_method = {0: 'snappy', 1: 'zlib'}[idx]
    df_tmp = df_melted[df_melted['compression_method'] == c_method]
    for i, file_type in enumerate(df_tmp['file_type'].unique()):
        for j, variable in enumerate(df_tmp['variable'].unique()):
            x = i - 0.35 + (j * (0.8 / 2.0))
            y = df_tmp[df_tmp['file_type'] == file_type][
                df_tmp['variable'] == variable]['value'].mean() * 1.1
            label = "{:.2f}".format(df_tmp[df_tmp['file_type'] == file_type][
                df_tmp['variable'] == variable]['value'].mean())
            g.text(x, y, label)

snsplot._legend.set_title("")
# snsplot.set(xscale='log')
snsplot.set(yscale='log')
snsplot = snsplot.set_axis_labels(
    "Tipo de arquivo", "Tamanho do arquivo log10 (MB)").set_titles(
        template=" 1m SNPs - {col_var} = {col_name}", )
snsplot.savefig(graph_dir + 'experiment1_sizes_1m.png')
plt.draw()

# %% [markdown]
# ## Comentários
# - O tempo de execução resultante aparenta condizer com o apresentado na
# monografia do Rodrigo.
# - Quanto ao espaço necessário de armazenamento, o método de compressão zlib
# apresentou um fator de compressão relativamente maior, quando comparado ao
# snappy. Contudo, em alguns casos, o espaço necessário para armazenamento
# pós-importação ainda é maior que o tamanho original dos arquivos de dados.

# %% [markdown]
# ## Experimento 2.1 - Tempos de importação para grandes volumes de dados
# Analisando performance considerando os casos com menor e maior tempo de
# importação da seção 1,assim como o caso F (tempo de importação para todos os
# tipos de arquivos).
# - Caso menor: 0125
# - Caso maior: VCF

# %% [markdown]
# ### Tempos de execução

# %%

df_melted = pd.melt(df[df['experiment_id'].str.startswith('2')],
                    id_vars=[
                        'experiment_id', 'compression_method', 'file_type',
                        'nsnps', 'nsamples'
                    ],
                    value_vars=['time'])
sns.set(style="whitegrid",
        palette=sns.color_palette("muted", n_colors=6, desat=1.0))
snsplot = sns.lmplot(data=df_melted[df_melted['variable'] == 'time'],
                     y='value',
                     x='nsamples',
                     col='file_type',
                     order=1,
                     height=8,
                     aspect=1,
                     legend_out=False,
                     truncate=False,
                     sharex=True,
                     sharey=False)

# Labelling points of graph
for idx, g in enumerate(snsplot.axes[0]):
    file_type = {0: '0125', 1: 'VCF'}[idx]
    df_tmp = df_melted[df_melted['file_type'] == file_type]
    for row_idx, row in df_tmp.loc[:, ['value', 'nsamples']].iterrows():
        x = row['nsamples'] * 0.6
        y = row['value'] + 0.05 * df_tmp['value'].max()
        label = "{:.3f}".format(row['value'])
        g.text(x, y, label)

# snsplot._legend.set_title("Tipo de arquivo")
snsplot.set(xscale='log')
# snsplot.set(yscale='log')
snsplot = snsplot.set_axis_labels("# de amostras",
                                  "Tempo de execução (s)").set_titles(
                                      template="Tipo de arquivo = {col_name}")
snsplot.savefig(graph_dir + 'experiment2_1_times.png')
plt.draw()

# %% [markdown]
# ### Tamanhos de arquivos antes e depois da inserção

# %%
df_melted = pd.melt(df[df['experiment_id'].str.startswith('2')],
                    id_vars=[
                        'experiment_id', 'compression_method', 'file_type',
                        'nsnps', 'nsamples'
                    ],
                    value_vars=['fsize', 'dbsize'])
df_melted['value'] /= 1024**2
sns.set(style="whitegrid",
        palette=sns.color_palette("muted", n_colors=6, desat=1.0))
snsplot = sns.lmplot(data=df_melted[df_melted['variable'].isin(
    ['fsize', 'dbsize'])],
                     y='value',
                     x='nsamples',
                     hue='variable',
                     col='file_type',
                     order=1,
                     height=8,
                     aspect=1,
                     legend_out=True,
                     truncate=False,
                     sharex=True,
                     sharey=False)
snsplot._legend.set_title("")
snsplot.set(xscale='log')
snsplot = snsplot.set_axis_labels("# de amostras",
                                  "Tamanho do arquivo (MB)").set_titles(
                                      template="Tipo de arquivo = {col_name}")
snsplot.savefig(graph_dir + 'experiment2_1_sizes.png')
plt.draw()

# %% [markdown]
# ### Comentários
# - Os gráficos são referentes aos resultados parciais do experimento 2.1,
# i.e., com execuções com 100 mil marcadores de SNPs, até 1 milhão de amostras
# para o tipo de arquivo 0125 e até 100 mil amostras para o tipo de arquivo
# VCF;
# - Para cada gráfico, foi feita uma regressão linear, e foi traçada a
# linha resultante;
# - Considerando o gráfico de tamanhos de arquivo, as legendas "fsize" e
# "dbsize" referem-se, respectivamente, aos tamanhos dos arquivos de dados
# antes e depois da inserção na base de dados.

# %% [markdown]
# ## Experimento 2.2 - Sumarização
# Análise de sumarização de dados, dado um indivíduo presente na base de dados.

# %% [markdown]
# ### Utilizando biblioteca snpdb para sumarização de um indíviduo

# %%
inds = snpdb.find_individuals()
ind = np.random.choice(inds, size=2, replace=True)
summary = snpdb.summarize(ind[0])
pp(summary)

# %% [markdown]
# ## Experimento 2.3 - Exportação de sumarização para formatos originais
# - Em progresso

# %%

# %% [markdown]
# ## Experimento 2.4 - Busca de indivíduos, dada uma lista de SNPs
# Dada uma lista de marcadores, deseja-se buscar todos os animais que os
# possuem, nos diferentes tipos dearquivos.

# %% [markdown]
# ### Utilizando biblioteca snpdb para busca de indivíduos

# %%
snp = np.random.choice(snpdb.find_snp(), size=1)
inds = snpdb.find_individuals_of_snps(id=snp[0]['i'])
pp(inds)

# %% [markdown]
# ## Experimento 2.5 - Exportação de indivíduos para formatos originais
# - Em progresso

# %%

# %% [markdown]
# ## Experimento 2.6 - Remoção de todos os dados de um indivíduo
# Dada uma lista de indivíduos, remover todos os dados associados

# %% [markdown]
# ### Utilizando biblioteca snpdb para deleção

# %%
ind = np.random.choice(inds, size=1)[0]
delete_results = snpdb.delete_individuals(id=ind['_id'])
pp(delete_results)
print("Deleted counts: ", end='')
pp([i.deleted_count for i in delete_results])

# %% [markdown]
# ## Experimento 2.7 - Comparação de importação/exportação no modo “bruto”
# - Em progresso: Script está feito, mas aguardando finalização do experimento
# 2.1

# %%

# %%
