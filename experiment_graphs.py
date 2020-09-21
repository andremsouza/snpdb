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
# from pprint import pp
import seaborn as sns
# import snpdb

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

                # Getting result values lists from json for insert operation
                fsizes: list = []
                dbsizes: list = []
                times: list = []
                summarize_times: list = []
                individuals_of_snps_times: list = []
                delete_individual_times: list = []
                if exps[experiment_id]['file_type'] in bin_file_types:
                    fsizes = results[experiment_id][compression_method][
                        'fsize']
                    dbsizes = results[experiment_id][compression_method][
                        'dbsize']
                    times = results[experiment_id][compression_method]['time']
                    summarize_times = [0.0] * len(fsizes)
                    individuals_of_snps_times = [0.0] * len(fsizes)
                    delete_individual_times = [0.0] * len(fsizes)
                else:
                    fsizes = results[experiment_id][compression_method][
                        nsnps_id][nsamples_id]['fsize']
                    dbsizes = results[experiment_id][compression_method][
                        nsnps_id][nsamples_id]['dbsize']
                    times = results[experiment_id][compression_method][
                        nsnps_id][nsamples_id]['time']
                    # Getting results values for summarize operation
                    summarize_times = [
                        summ['time']
                        for summ in results[experiment_id][compression_method]
                        [nsnps_id][nsamples_id]['summarize']
                    ]
                    for ind in results[experiment_id][compression_method][
                            nsnps_id][nsamples_id]['individuals_of_snps']:
                        try:
                            # Fix KeyError in results
                            individuals_of_snps_times.append(float(
                                ind['time']))
                        except KeyError:
                            individuals_of_snps_times.append(float(ind['snp']))
                    delete_individual_times = [
                        ind['time']
                        for ind in results[experiment_id][compression_method]
                        [nsnps_id][nsamples_id]['delete_individual']
                    ]

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
                    'summarize': float(l),
                    'individuals_of_snps': float(m),
                    'delete_individual': float(n)
                } for i, j, k, l, m, n in zip(
                    fsizes, dbsizes, times, summarize_times,
                    individuals_of_snps_times, delete_individual_times)]
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
# importação da seção 1
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
snsplot = sns.lmplot(data=df_melted,
                     y='value',
                     x='nsamples',
                     row='file_type',
                     col='compression_method',
                     col_order=['snappy', 'zlib'],
                     hue='nsnps',
                     order=1,
                     height=8,
                     aspect=1,
                     legend_out=False,
                     truncate=False,
                     sharex=True,
                     sharey=False)

# # Labelling points of graph
# for idx, g in enumerate(snsplot.axes[0]):
#     file_type = {0: '0125', 1: 'VCF'}[idx]
#     df_tmp = df_melted[df_melted['file_type'] == file_type]
#     for row_idx, row in df_tmp.loc[:, ['value', 'nsamples']].iterrows():
#         x = row['nsamples'] * 0.6
#         y = row['value'] + 0.05 * df_tmp['value'].max()
#         label = "{:.3f}".format(row['value'])
#         g.text(x, y, label)

# snsplot._legend.set_title("Tipo de arquivo")
# snsplot.set(xscale='log')
# snsplot.set(yscale='log')
snsplot = snsplot.set_axis_labels(
    "# de amostras", "Tempo de execução (s)"
).set_titles(
    template="Tipo de arquivo = {row_name}; Método de compressão = {col_name}")
snsplot.savefig(graph_dir + 'experiment2_1_times.png')
plt.draw()

# %% [markdown]
# ### Tamanhos de arquivos antes e depois da inserção (100k SNPs)

# %%
df_melted = pd.melt(df[df['experiment_id'].str.startswith('2')],
                    id_vars=[
                        'experiment_id', 'compression_method', 'file_type',
                        'nsnps', 'nsamples'
                    ],
                    value_vars=['fsize', 'dbsize'])
df_melted = df_melted[df_melted['nsnps'] == 100000.0]
df_melted['value'] /= 1024**2
sns.set(style="whitegrid",
        palette=sns.color_palette("muted", n_colors=6, desat=1.0))
snsplot = sns.lmplot(data=df_melted[df_melted['variable'].isin(
    ['fsize', 'dbsize'])],
                     y='value',
                     x='nsamples',
                     hue='variable',
                     row='file_type',
                     col='compression_method',
                     col_order=['snappy', 'zlib'],
                     order=1,
                     height=8,
                     aspect=1,
                     legend_out=True,
                     truncate=False,
                     sharex=True,
                     sharey=False)
snsplot._legend.set_title("")
# snsplot.set(xscale='log')
snsplot = snsplot.set_axis_labels(
    "# de amostras", "Tamanho do arquivo (MB)"
).set_titles(
    template="Tipo de arquivo = {row_name}; Método de compressão = {col_name}")
snsplot.savefig(graph_dir + 'experiment2_1_sizes_100k.png')
plt.draw()

# %% [markdown]
# ### Tamanhos de arquivos antes e depois da inserção (1m SNPs)

# %%
df_melted = pd.melt(df[df['experiment_id'].str.startswith('2')],
                    id_vars=[
                        'experiment_id', 'compression_method', 'file_type',
                        'nsnps', 'nsamples'
                    ],
                    value_vars=['fsize', 'dbsize'])
df_melted = df_melted[df_melted['nsnps'] == 1000000.0]
df_melted['value'] /= 1024**2
sns.set(style="whitegrid",
        palette=sns.color_palette("muted", n_colors=6, desat=1.0))
snsplot = sns.lmplot(data=df_melted[df_melted['variable'].isin(
    ['fsize', 'dbsize'])],
                     y='value',
                     x='nsamples',
                     hue='variable',
                     row='file_type',
                     col='compression_method',
                     col_order=['snappy', 'zlib'],
                     order=1,
                     height=8,
                     aspect=1,
                     legend_out=True,
                     truncate=False,
                     sharex=True,
                     sharey=False)
snsplot._legend.set_title("")
# snsplot.set(xscale='log')
snsplot = snsplot.set_axis_labels(
    "# de amostras", "Tamanho do arquivo (MB)"
).set_titles(
    template="Tipo de arquivo = {row_name}; Método de compressão = {col_name}")
snsplot.savefig(graph_dir + 'experiment2_1_sizes_1m.png')
plt.draw()

# %% [markdown]
# ### Comentários
# - Os gráficos são referentes aos resultados do experimento 2.1,
# i.e., com execuções com 100 mil e 1 milhão de marcadores de SNPs,
# até 100 mil de amostras e até 10 mil amostras para o tipo de arquivo VCF
# - Para cada gráfico, foi feita uma regressão linear, e foi traçada a
# linha resultante, e seu intervale de confiança com C=0.95;
# - Considerando o gráfico de tamanhos de arquivo, as legendas "fsize" e
# "dbsize" referem-se, respectivamente, aos tamanhos dos arquivos de dados
# antes e depois da inserção na base de dados.
# - Verificar possíveis melhorias para a apresentação dos gráficos dessa seção

# %% [markdown]
# ## Experimento 2.2 - Sumarização
# Análise de sumarização de dados, dado um indivíduo presente na base de dados.

# %%
df_melted = pd.melt(df[df['experiment_id'].str.startswith('2')],
                    id_vars=[
                        'experiment_id', 'compression_method', 'file_type',
                        'nsnps', 'nsamples'
                    ],
                    value_vars=['summarize'])
sns.set(style="whitegrid",
        palette=sns.color_palette("muted", n_colors=6, desat=1.0))
snsplot = sns.lmplot(data=df_melted,
                     y='value',
                     x='nsamples',
                     row='file_type',
                     col='compression_method',
                     col_order=['snappy', 'zlib'],
                     hue='nsnps',
                     order=1,
                     height=8,
                     aspect=1,
                     legend_out=False,
                     truncate=False,
                     sharex=True,
                     sharey=False)

# # Labelling points of graph
# for idx, g in enumerate(snsplot.axes[0]):
#     file_type = {0: '0125', 1: 'VCF'}[idx]
#     df_tmp = df_melted[df_melted['file_type'] == file_type]
#     for row_idx, row in df_tmp.loc[:, ['value', 'nsamples']].iterrows():
#         x = row['nsamples'] * 0.6
#         y = row['value'] + 0.05 * df_tmp['value'].max()
#         label = "{:.3f}".format(row['value'])
#         g.text(x, y, label)

# snsplot._legend.set_title("Tipo de arquivo")
# snsplot.set(xscale='log')
# snsplot.set(yscale='log')
snsplot = snsplot.set_axis_labels(
    "# de amostras", "Tempo de execução (s)"
).set_titles(
    template="Tipo de arquivo = {row_name}; Método de compressão = {col_name}")
snsplot.savefig(graph_dir + 'experiment2_1_times.png')
plt.draw()

# %% [markdown]
# ## Experimento 2.3 - Exportação de sumarização para formatos originais
# - Em progresso

# %%
# TODO

# %% [markdown]
# ## Experimento 2.4 - Busca de indivíduos, dada uma lista de SNPs
# Dada uma lista de marcadores, deseja-se buscar todos os animais que os
# possuem, nos diferentes tipos dearquivos.

# %%
df_melted = pd.melt(df[df['experiment_id'].str.startswith('2')],
                    id_vars=[
                        'experiment_id', 'compression_method', 'file_type',
                        'nsnps', 'nsamples'
                    ],
                    value_vars=['individuals_of_snps'])
sns.set(style="whitegrid",
        palette=sns.color_palette("muted", n_colors=6, desat=1.0))
snsplot = sns.lmplot(data=df_melted,
                     y='value',
                     x='nsamples',
                     row='file_type',
                     col='compression_method',
                     col_order=['snappy', 'zlib'],
                     hue='nsnps',
                     order=1,
                     height=8,
                     aspect=1,
                     legend_out=False,
                     truncate=False,
                     sharex=True,
                     sharey=False)

# # Labelling points of graph
# for idx, g in enumerate(snsplot.axes[0]):
#     file_type = {0: '0125', 1: 'VCF'}[idx]
#     df_tmp = df_melted[df_melted['file_type'] == file_type]
#     for row_idx, row in df_tmp.loc[:, ['value', 'nsamples']].iterrows():
#         x = row['nsamples'] * 0.6
#         y = row['value'] + 0.05 * df_tmp['value'].max()
#         label = "{:.3f}".format(row['value'])
#         g.text(x, y, label)

# snsplot._legend.set_title("Tipo de arquivo")
# snsplot.set(xscale='log')
# snsplot.set(yscale='log')
snsplot = snsplot.set_axis_labels(
    "# de amostras", "Tempo de execução (s)"
).set_titles(
    template="Tipo de arquivo = {row_name}; Método de compressão = {col_name}")
snsplot.savefig(graph_dir + 'experiment2_1_times.png')
plt.draw()

# %% [markdown]
# ## Experimento 2.5 - Exportação de indivíduos para formatos originais
# - Em progresso

# %%
# TODO

# %% [markdown]
# ## Experimento 2.6 - Remoção de todos os dados de um indivíduo
# Dada uma lista de indivíduos, remover todos os dados associados

# %%
df_melted = pd.melt(df[df['experiment_id'].str.startswith('2')],
                    id_vars=[
                        'experiment_id', 'compression_method', 'file_type',
                        'nsnps', 'nsamples'
                    ],
                    value_vars=['delete_individual'])
sns.set(style="whitegrid",
        palette=sns.color_palette("muted", n_colors=6, desat=1.0))
snsplot = sns.lmplot(data=df_melted,
                     y='value',
                     x='nsamples',
                     row='file_type',
                     col='compression_method',
                     col_order=['snappy', 'zlib'],
                     hue='nsnps',
                     order=1,
                     height=8,
                     aspect=1,
                     legend_out=False,
                     truncate=False,
                     sharex=True,
                     sharey=False)

# # Labelling points of graph
# for idx, g in enumerate(snsplot.axes[0]):
#     file_type = {0: '0125', 1: 'VCF'}[idx]
#     df_tmp = df_melted[df_melted['file_type'] == file_type]
#     for row_idx, row in df_tmp.loc[:, ['value', 'nsamples']].iterrows():
#         x = row['nsamples'] * 0.6
#         y = row['value'] + 0.05 * df_tmp['value'].max()
#         label = "{:.3f}".format(row['value'])
#         g.text(x, y, label)

# snsplot._legend.set_title("Tipo de arquivo")
# snsplot.set(xscale='log')
# snsplot.set(yscale='log')
snsplot = snsplot.set_axis_labels(
    "# de amostras", "Tempo de execução (s)"
).set_titles(
    template="Tipo de arquivo = {row_name}; Método de compressão = {col_name}")
snsplot.savefig(graph_dir + 'experiment2_1_times.png')
plt.draw()

# %% [markdown]
# ## Experimento 2.7 - Comparação de importação/exportação no modo “bruto”
# - Em progresso: Script está feito, mas aguardando finalização do experimento
# 2.1

# %%
# TODO

# %%
