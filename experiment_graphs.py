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
    'fsize', 'dbsize', 'time', 'summarize', 'individuals_of_snps',
    'delete_individual', 'export_time_0125', 'export_time_plink',
    'export_time_bin'
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
                export_time_0125: list = []
                export_time_plink: list = []
                export_time_bin: list = []
                if exps[experiment_id]['file_type'] in bin_file_types:
                    fsizes = results[experiment_id][compression_method][
                        'fsize']
                    dbsizes = results[experiment_id][compression_method][
                        'dbsize']
                    times = results[experiment_id][compression_method]['time']
                    summarize_times = [0.0] * len(fsizes)
                    individuals_of_snps_times = [0.0] * len(fsizes)
                    delete_individual_times = [0.0] * len(fsizes)
                    export_time_0125 = [0.0] * len(fsizes)
                    export_time_plink = [0.0] * len(fsizes)
                    export_time_bin = [0.0] * len(fsizes)
                elif experiment_id == '2.7':
                    fsizes = results[experiment_id][compression_method][
                        nsnps_id][nsamples_id]['fsize']
                    dbsizes = results[experiment_id][compression_method][
                        nsnps_id][nsamples_id]['dbsize']
                    times = results[experiment_id][compression_method][
                        nsnps_id][nsamples_id]['time']
                    summarize_times = [0.0] * len(fsizes)
                    individuals_of_snps_times = [0.0] * len(fsizes)
                    delete_individual_times = [0.0] * len(fsizes)
                    export_time_0125 = [0.0] * len(fsizes)
                    export_time_plink = [0.0] * len(fsizes)
                    export_time_bin = [0.0] * len(fsizes)
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
                    if exps[experiment_id]['file_type'] == 'ALL':
                        export_time_0125 = [
                            i['Z125']['map'] + i['Z125']['samples']
                            for i in results[experiment_id][compression_method]
                            [nsnps_id][nsamples_id]['export']
                        ]
                        export_time_plink = [
                            i['PLINK']['map'] + i['PLINK']['samples']
                            for i in results[experiment_id][compression_method]
                            [nsnps_id][nsamples_id]['export']
                        ]
                        export_time_bin = [
                            i
                            for i in results[experiment_id][compression_method]
                            [nsnps_id][nsamples_id]['export_bin']
                        ]
                    else:
                        export_time_0125 = [0.0] * len(fsizes)
                        export_time_plink = [0.0] * len(fsizes)
                        export_time_bin = [0.0] * len(fsizes)

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
                    'delete_individual': float(n),
                    'export_time_0125': float(o),
                    'export_time_plink': float(p),
                    'export_time_bin': float(q)
                } for i, j, k, l, m, n, o, p, q in zip(
                    fsizes, dbsizes, times, summarize_times,
                    individuals_of_snps_times, delete_individual_times,
                    export_time_0125, export_time_plink, export_time_bin)]
                df = df.append(rows, ignore_index=True, verify_integrity=True)

# %% [markdown]
# #### Criando colunas concatenadas
df['file_type_nsnps'] = df['file_type'].map(str) + '_' + df['nsnps'].map(str)

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
df_melted['nsnps'][df_melted['nsnps'] == 1000000] = nsnps_ids[float(1000000)]
df_melted['nsnps'][df_melted['nsnps'] == 100000] = nsnps_ids[float(100000)]
# df_melted['nsnps'] = df_melted['nsnps'].str.replace('1000000', '1m')
# df_melted['nsnps'] = df_melted['nsnps'].str.replace('100000', '100k')
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
# ### Comparação de tamanhos de arquivos antes/depois das inserções

# %%
df_melted = pd.melt(df[df['experiment_id'].str.startswith('1')],
                    id_vars=[
                        'experiment_id', 'compression_method', 'file_type',
                        'nsnps', 'nsamples'
                    ],
                    value_vars=['fsize', 'dbsize'])
df_melted = df_melted.fillna(value='bin')
df_melted['variable'] = df_melted['variable'].map(
    str) + '_' + df_melted['nsnps'].map(str)
df_melted['variable'] = df_melted['variable'].str.replace(
    '100000.0', '100k').str.replace('1000000.0', '1m')
df_melted['value'] /= 1024**2
sns.set(style="whitegrid",
        palette=sns.color_palette("muted", n_colors=6, desat=1.0))
snsplot = sns.catplot(
    x='file_type',
    y='value',
    hue='variable',
    hue_order=[
        'fsize_100k', 'dbsize_100k', 'fsize_1m', 'dbsize_1m', 'fsize_bin',
        'dbsize_bin'
    ],
    data=df_melted,
    col='compression_method',
    col_order=['snappy', 'zlib'],
    kind='bar',
    height=8,
    aspect=1,
)

# Labelling bars of graph
for idx, g in enumerate(snsplot.axes[0]):
    for p in g.patches:
        if not np.isnan(p.get_height()):
            g.text(p.get_x() + p.get_width() / 2.,
                   p.get_height(),
                   '%.1f' % float(p.get_height()),
                   fontsize=8,
                   color='black',
                   ha='center',
                   va='bottom')

snsplot._legend.set_title("size_nsnps")
# snsplot.set(xscale='log')
snsplot.set(yscale='log')
snsplot = snsplot.set_axis_labels("Tipo de arquivo",
                                  "Tamanho do arquivo log10 (MB)").set_titles(
                                      template="{col_var} = {col_name}", )
snsplot.savefig(graph_dir + 'experiment1_sizes.png')
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
# - Caso médio: PLINK
# - Caso maior: VCF

# %% [markdown]
# ### Tempos de importação

# %%
df_melted = pd.melt(df[df['experiment_id'].str.startswith('2')],
                    id_vars=[
                        'experiment_id', 'compression_method', 'file_type',
                        'nsnps', 'nsamples'
                    ],
                    value_vars=['time'])
sns.set(style="whitegrid",
        palette=sns.color_palette("muted", n_colors=6, desat=1.0))

fig, ax = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(16, 9))
ax[0].set_xscale('log')
ax[0].set_yscale('log')
sns.lineplot(data=df_melted[df['compression_method'] == 'snappy'],
             x='nsamples',
             y='value',
             hue='file_type',
             hue_order=['0125', 'PLINK', 'VCF', 'ALL'],
             palette='tab10',
             style='nsnps',
             ci=None,
             legend='full',
             figure=fig,
             ax=ax[0])
sns.lineplot(data=df_melted[df['compression_method'] == 'zlib'],
             x='nsamples',
             y='value',
             hue='file_type',
             hue_order=['0125', 'PLINK', 'VCF', 'ALL'],
             palette='tab10',
             style='nsnps',
             ci=None,
             legend='full',
             figure=fig,
             ax=ax[1])

ax[0].set_title('compression_method = snappy')
ax[1].set_title('compression_method = zlib')
ax[0].set_xlabel('# de amostras (log10)')
ax[1].set_xlabel('# de amostras (log10)')
ax[0].set_ylabel('Tempo de importação (log10) (s)')
ax[1].set_ylabel('Tempo de importação (log10) (s)')

for idx, g in enumerate(ax):
    c_methods = {0: 'snappy', 1: 'zlib'}
    for nsnps in df_melted['nsnps'].unique():
        for file_type in df_melted['file_type'].unique():
            df_tmp = df_melted[df_melted['file_type'] == file_type][
                df_melted['compression_method'] == c_methods[idx]][
                    df_melted['nsnps'] == nsnps]
            g.text(df_tmp['nsamples'].max(),
                   df_tmp['value'].max(),
                   '%.1f' % float(df_tmp['value'].max()),
                   fontsize=8,
                   color='black',
                   ha='center',
                   va='bottom')
fig.savefig(graph_dir + 'experiment2_1_times.png')
plt.draw()

# %% [markdown]
# ### Tamanhos de arquivos antes e depois da inserção (100k)

# %%
df_melted = pd.melt(df[df['experiment_id'].str.startswith('2')],
                    id_vars=[
                        'experiment_id', 'compression_method', 'file_type',
                        'nsnps', 'nsamples'
                    ],
                    value_vars=['fsize', 'dbsize'])
df_melted = df_melted.fillna(value='bin')
df_melted = df_melted[df_melted['nsnps'] == 100000.0]
# df_melted['variable'] = df_melted['variable'].map(
#     str) + '_' + df_melted['nsnps'].map(str)
# df_melted['variable'] = df_melted['variable'].str.replace(
#     '100000.0', '100k').str.replace('1000000.0', '1m')
df_melted['value'] /= 1024**2
sns.set(style="whitegrid",
        palette=sns.color_palette("muted", n_colors=6, desat=1.0))

fig, ax = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(16, 9))
ax[0].set_xscale('log')
ax[0].set_yscale('log')
sns.lineplot(data=df_melted[df['compression_method'] == 'snappy'],
             x='nsamples',
             y='value',
             hue='file_type',
             hue_order=['0125', 'PLINK', 'VCF', 'ALL'],
             palette='tab10',
             style='variable',
             style_order=['fsize', 'dbsize'],
             ci=None,
             legend='full',
             figure=fig,
             ax=ax[0])
sns.lineplot(data=df_melted[df['compression_method'] == 'zlib'],
             x='nsamples',
             y='value',
             hue='file_type',
             hue_order=['0125', 'PLINK', 'VCF', 'ALL'],
             palette='tab10',
             style='variable',
             style_order=['fsize', 'dbsize'],
             ci=None,
             legend='full',
             figure=fig,
             ax=ax[1])

ax[0].set_title('compression_method = snappy')
ax[1].set_title('compression_method = zlib')
ax[0].set_xlabel('# de amostras (log10)')
ax[1].set_xlabel('# de amostras (log10)')
ax[0].set_ylabel('Tamanho (log10) (MB)')
ax[1].set_ylabel('Tamanho (log10) (MB)')

# for idx, g in enumerate(ax):
#     c_methods = {0: 'snappy', 1: 'zlib'}
#     for nsnps in df_melted['nsnps'].unique():
#         for file_type in df_melted['file_type'].unique():
#             df_tmp = df_melted[df_melted['file_type'] == file_type][
#                 df_melted['compression_method'] == c_methods[idx]][
#                     df_melted['nsnps'] == nsnps]
#             g.text(df_tmp['nsamples'].max(),
#                    df_tmp['value'].max(),
#                    '%.1f' % float(df_tmp['value'].max()),
#                    fontsize=8,
#                    color='black',
#                    ha='center',
#                    va='bottom')
fig.savefig(graph_dir + 'experiment2_1_sizes_100k.png')
plt.draw()

# %% [markdown]
# ### Tamanhos de arquivos antes e depois da inserção (1m)

# %%
df_melted = pd.melt(df[df['experiment_id'].str.startswith('2')],
                    id_vars=[
                        'experiment_id', 'compression_method', 'file_type',
                        'nsnps', 'nsamples'
                    ],
                    value_vars=['fsize', 'dbsize'])
df_melted = df_melted.fillna(value='bin')
df_melted = df_melted[df_melted['nsnps'] == 1000000.0]
# df_melted['variable'] = df_melted['variable'].map(
#     str) + '_' + df_melted['nsnps'].map(str)
# df_melted['variable'] = df_melted['variable'].str.replace(
#     '100000.0', '100k').str.replace('1000000.0', '1m')
df_melted['value'] /= 1024**2
sns.set(style="whitegrid",
        palette=sns.color_palette("muted", n_colors=6, desat=1.0))

fig, ax = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(16, 9))
ax[0].set_xscale('log')
ax[0].set_yscale('log')
sns.lineplot(data=df_melted[df['compression_method'] == 'snappy'],
             x='nsamples',
             y='value',
             hue='file_type',
             hue_order=['0125', 'PLINK', 'VCF', 'ALL'],
             palette='tab10',
             style='variable',
             style_order=['fsize', 'dbsize'],
             ci=None,
             legend='full',
             figure=fig,
             ax=ax[0])
sns.lineplot(data=df_melted[df['compression_method'] == 'zlib'],
             x='nsamples',
             y='value',
             hue='file_type',
             hue_order=['0125', 'PLINK', 'VCF', 'ALL'],
             palette='tab10',
             style='variable',
             style_order=['fsize', 'dbsize'],
             ci=None,
             legend='full',
             figure=fig,
             ax=ax[1])

ax[0].set_title('compression_method = snappy')
ax[1].set_title('compression_method = zlib')
ax[0].set_xlabel('# de amostras (log10)')
ax[1].set_xlabel('# de amostras (log10)')
ax[0].set_ylabel('Tamanho (log10) (MB)')
ax[1].set_ylabel('Tamanho (log10) (MB)')

# for idx, g in enumerate(ax):
#     c_methods = {0: 'snappy', 1: 'zlib'}
#     for nsnps in df_melted['nsnps'].unique():
#         for file_type in df_melted['file_type'].unique():
#             df_tmp = df_melted[df_melted['file_type'] == file_type][
#                 df_melted['compression_method'] == c_methods[idx]][
#                     df_melted['nsnps'] == nsnps]
#             g.text(df_tmp['nsamples'].max(),
#                    df_tmp['value'].max(),
#                    '%.1f' % float(df_tmp['value'].max()),
#                    fontsize=8,
#                    color='black',
#                    ha='center',
#                    va='bottom')
fig.savefig(graph_dir + 'experiment2_1_sizes_1m.png')
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

fig, ax = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(16, 9))
ax[0].set_xscale('log')
ax[0].set_yscale('log')
sns.lineplot(data=df_melted[df['compression_method'] == 'snappy'],
             x='nsamples',
             y='value',
             hue='file_type',
             hue_order=['0125', 'PLINK', 'VCF', 'ALL'],
             palette='tab10',
             style='nsnps',
             ci=None,
             legend='full',
             figure=fig,
             ax=ax[0])
sns.lineplot(data=df_melted[df['compression_method'] == 'zlib'],
             x='nsamples',
             y='value',
             hue='file_type',
             hue_order=['0125', 'PLINK', 'VCF', 'ALL'],
             palette='tab10',
             style='nsnps',
             ci=None,
             legend='full',
             figure=fig,
             ax=ax[1])

ax[0].set_title('compression_method = snappy')
ax[1].set_title('compression_method = zlib')
ax[0].set_xlabel('# de amostras (log10)')
ax[1].set_xlabel('# de amostras (log10)')
ax[0].set_ylabel('Tempo de execução (log10) (s)')
ax[1].set_ylabel('Tempo de execução (log10) (s)')

for idx, g in enumerate(ax):
    c_methods = {0: 'snappy', 1: 'zlib'}
    for nsnps in df_melted['nsnps'].unique():
        for file_type in df_melted['file_type'].unique():
            df_tmp = df_melted[df_melted['file_type'] == file_type][
                df_melted['compression_method'] == c_methods[idx]][
                    df_melted['nsnps'] == nsnps]
            g.text(df_tmp['nsamples'].max(),
                   df_tmp['value'].max(),
                   '%.1f' % float(df_tmp['value'].max()),
                   fontsize=8,
                   color='black',
                   ha='center',
                   va='bottom')
fig.savefig(graph_dir + 'experiment2_summarize_times.png')
plt.draw()

# %% [markdown]
# ## Experimento 2.3 - Exportação de dados para formatos originais
# - Exportação de um par mapa/amostra de um indivíduo, escolhido
# aleatóriamente, para os formatos 0125 e PLINK.

# %%
df_melted = pd.melt(df[df['file_type'] == 'ALL'],
                    id_vars=[
                        'experiment_id', 'compression_method', 'file_type',
                        'nsnps', 'nsamples'
                    ],
                    value_vars=['export_time_0125', 'export_time_plink'])
sns.set(style="whitegrid",
        palette=sns.color_palette("muted", n_colors=6, desat=1.0))

fig, ax = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(16, 9))
# ax[0].set_xscale('log')
# ax[0].set_yscale('log')
sns.lineplot(data=df_melted[df['compression_method'] == 'snappy'],
             x='nsamples',
             y='value',
             hue='variable',
             palette='tab10',
             style='nsnps',
             ci=None,
             legend='full',
             figure=fig,
             ax=ax[0])
sns.lineplot(data=df_melted[df['compression_method'] == 'zlib'],
             x='nsamples',
             y='value',
             hue='variable',
             palette='tab10',
             style='nsnps',
             ci=None,
             legend='full',
             figure=fig,
             ax=ax[1])

ax[0].set_title('compression_method = snappy')
ax[1].set_title('compression_method = zlib')
ax[0].set_xlabel('# de amostras')
ax[1].set_xlabel('# de amostras')
ax[0].set_ylabel('Tempo de execução (s)')
ax[1].set_ylabel('Tempo de execução (s)')

# for idx, g in enumerate(ax):
#     c_methods = {0: 'snappy', 1: 'zlib'}
#     for nsnps in df_melted['nsnps'].unique():
#         for file_type in df_melted['file_type'].unique():
#             df_tmp = df_melted[df_melted['file_type'] == file_type][
#                 df_melted['compression_method'] == c_methods[idx]][
#                     df_melted['nsnps'] == nsnps]
#             g.text(df_tmp['nsamples'].max(),
#                    df_tmp['value'].max(),
#                    '%.1f' % float(df_tmp['value'].max()),
#                    fontsize=8,
#                    color='black',
#                    ha='center',
#                    va='bottom')
fig.savefig(graph_dir + 'experiment2_export_times.png')
plt.draw()

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

fig, ax = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(16, 9))
ax[0].set_xscale('log')
ax[0].set_yscale('log')
sns.lineplot(data=df_melted[df['compression_method'] == 'snappy'],
             x='nsamples',
             y='value',
             hue='file_type',
             hue_order=['0125', 'PLINK', 'VCF', 'ALL'],
             palette='tab10',
             style='nsnps',
             ci=None,
             legend='full',
             figure=fig,
             ax=ax[0])
sns.lineplot(data=df_melted[df['compression_method'] == 'zlib'],
             x='nsamples',
             y='value',
             hue='file_type',
             hue_order=['0125', 'PLINK', 'VCF', 'ALL'],
             palette='tab10',
             style='nsnps',
             ci=None,
             legend='full',
             figure=fig,
             ax=ax[1])

ax[0].set_title('compression_method = snappy')
ax[1].set_title('compression_method = zlib')
ax[0].set_xlabel('# de amostras (log10)')
ax[1].set_xlabel('# de amostras (log10)')
ax[0].set_ylabel('Tempo de execução (log10) (s)')
ax[1].set_ylabel('Tempo de execução (log10) (s)')

for idx, g in enumerate(ax):
    c_methods = {0: 'snappy', 1: 'zlib'}
    for nsnps in df_melted['nsnps'].unique():
        for file_type in df_melted['file_type'].unique():
            df_tmp = df_melted[df_melted['file_type'] == file_type][
                df_melted['compression_method'] == c_methods[idx]][
                    df_melted['nsnps'] == nsnps]
            g.text(df_tmp['nsamples'].max(),
                   df_tmp['value'].max(),
                   '%.1f' % float(df_tmp['value'].max()),
                   fontsize=8,
                   color='black',
                   ha='center',
                   va='bottom')
fig.savefig(graph_dir + 'experiment2_indiv_search_times.png')
plt.draw()

# %% [markdown]
# ## Experimento 2.5 - Exportação de dados brutos/binários
# - Exportação de todos os arquivos armazenados em modo "bruto" de um indivíduo
# aleatoriamente selecionado.
# - Para cada experimento, é armazena um arquivo FastQ e um arquivo de mídia.

# %%
df_melted = pd.melt(df[df['file_type'] == 'ALL'],
                    id_vars=[
                        'experiment_id', 'compression_method', 'file_type',
                        'nsnps', 'nsamples'
                    ],
                    value_vars=['export_time_bin'])
sns.set(style="whitegrid",
        palette=sns.color_palette("muted", n_colors=6, desat=1.0))

fig, ax = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(16, 9))
# ax[0].set_xscale('log')
# ax[0].set_yscale('log')
sns.lineplot(data=df_melted[df['compression_method'] == 'snappy'],
             x='nsamples',
             y='value',
             hue='nsnps',
             palette='tab10',
             ci=None,
             legend='full',
             figure=fig,
             ax=ax[0])
sns.lineplot(data=df_melted[df['compression_method'] == 'zlib'],
             x='nsamples',
             y='value',
             hue='nsnps',
             palette='tab10',
             ci=None,
             legend='full',
             figure=fig,
             ax=ax[1])

ax[0].set_title('compression_method = snappy')
ax[1].set_title('compression_method = zlib')
ax[0].set_xlabel('# de amostras')
ax[1].set_xlabel('# de amostras')
ax[0].set_ylabel('Tempo de execução (s)')
ax[1].set_ylabel('Tempo de execução (s)')

for idx, g in enumerate(ax):
    c_methods = {0: 'snappy', 1: 'zlib'}
    for nsnps in df_melted['nsnps'].unique():
        for file_type in df_melted['file_type'].unique():
            df_tmp = df_melted[df_melted['file_type'] == file_type][
                df_melted['compression_method'] == c_methods[idx]][
                    df_melted['nsnps'] == nsnps]
            g.text(df_tmp['nsamples'].max(),
                   df_tmp['value'].max(),
                   '%.1f' % float(df_tmp['value'].max()),
                   fontsize=8,
                   color='black',
                   ha='center',
                   va='bottom')
fig.savefig(graph_dir + 'experiment2_export_times.png')
plt.draw()

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

fig, ax = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(16, 9))
ax[0].set_xscale('log')
ax[0].set_yscale('log')
sns.lineplot(data=df_melted[df['compression_method'] == 'snappy'],
             x='nsamples',
             y='value',
             hue='file_type',
             hue_order=['0125', 'PLINK', 'VCF', 'ALL'],
             palette='tab10',
             style='nsnps',
             ci=None,
             legend='full',
             figure=fig,
             ax=ax[0])
sns.lineplot(data=df_melted[df['compression_method'] == 'zlib'],
             x='nsamples',
             y='value',
             hue='file_type',
             hue_order=['0125', 'PLINK', 'VCF', 'ALL'],
             palette='tab10',
             style='nsnps',
             ci=None,
             legend='full',
             figure=fig,
             ax=ax[1])

ax[0].set_title('compression_method = snappy')
ax[1].set_title('compression_method = zlib')
ax[0].set_xlabel('# de amostras (log10)')
ax[1].set_xlabel('# de amostras (log10)')
ax[0].set_ylabel('Tempo de execução (log10) (s)')
ax[1].set_ylabel('Tempo de execução (log10) (s)')

# for idx, g in enumerate(ax):
#     c_methods = {0: 'snappy', 1: 'zlib'}
#     for nsnps in df_melted['nsnps'].unique():
#         for file_type in df_melted['file_type'].unique():
#             df_tmp = df_melted[df_melted['file_type'] == file_type][
#                 df_melted['compression_method'] == c_methods[idx]][
#                     df_melted['nsnps'] == nsnps]
#             g.text(df_tmp['nsamples'].max(),
#                    df_tmp['value'].max(),
#                    '%.1f' % float(df_tmp['value'].max()),
#                    fontsize=8,
#                    color='black',
#                    ha='center',
#                    va='bottom')
fig.savefig(graph_dir + 'experiment2_remove_times.png')
plt.draw()

# %% [markdown]
# ## Experimento 2.7 - Comparação de importação/exportação no modo “bruto”

# %%
# %% [markdown]
# ### Tempos de execução

# %%
df_melted = pd.melt(df[df['experiment_id'].isin(
    ['2.A', '2.7'])][df['nsnps'] == 100000.0][df['nsamples'] == 10000.0],
                    id_vars=[
                        'experiment_id', 'compression_method', 'file_type',
                        'nsnps', 'nsamples'
                    ],
                    value_vars=['fsize', 'dbsize'])
sns.set(style="whitegrid",
        palette=sns.color_palette("muted", n_colors=6, desat=1.0))
snsplot = sns.catplot(
    x='compression_method',
    y='value',
    hue='experiment_id',
    data=df_melted,
    kind='bar',
    height=8,
    aspect=1,
)

snsplot._legend.set_title("Experimento")
# snsplot.set(xscale='log')
snsplot.set(yscale='log')
snsplot = snsplot.set_axis_labels("Método de compressão",
                                  "Tempo de execução log10 (s)")
# .set_titles(template="Método de compressão: {row_name} - {col_name}")
snsplot.savefig(graph_dir + 'experiment2_7_times.png')
plt.draw()

# %% [markdown]
# ### Comparação de tamanhos de arquivos antes/depois das inserções (100k SNPs)

# %%
df_melted = pd.melt(df[df['experiment_id'].isin(
    ['2.A', '2.7'])][df['nsnps'] == 100000.0][df['nsamples'] == 10000.0],
                    id_vars=[
                        'experiment_id', 'compression_method', 'file_type',
                        'nsnps', 'nsamples'
                    ],
                    value_vars=['fsize', 'dbsize'])
df_melted['value'] /= 1024**2
sns.set(style="whitegrid",
        palette=sns.color_palette("muted", n_colors=6, desat=1.0))
snsplot = sns.catplot(
    x='compression_method',
    y='value',
    hue='variable',
    data=df_melted,
    col='experiment_id',
    kind='bar',
    height=8,
    aspect=1,
)

# Labelling bars of graph
for idx, g in enumerate(snsplot.axes[0]):
    e_id = {0: '2.A', 1: '2.7'}[idx]
    df_tmp = df_melted[df_melted['experiment_id'] == e_id]
    for i, compression_method in enumerate(
            df_tmp['compression_method'].unique()):
        for j, variable in enumerate(df_tmp['variable'].unique()):
            x = i - 0.35 + (j * (0.8 / 2.0))
            y = df_tmp[df_tmp['compression_method'] == compression_method][
                df_tmp['variable'] == variable]['value'].mean() * 1.01
            label = "{:.2f}".format(
                df_tmp[df_tmp['compression_method'] == compression_method][
                    df_tmp['variable'] == variable]['value'].mean())
            g.text(x, y, label)

snsplot._legend.set_title("")
# snsplot.set(xscale='log')
snsplot.set(yscale='log')
snsplot = snsplot.set_axis_labels("Método de compressão",
                                  "Tamanho do arquivo log10 (MB)").set_titles(
                                      template="{col_var} = {col_name}", )
snsplot.savefig(graph_dir + 'experiment2_7_sizes.png')
plt.draw()

# %%

# %%

# %%

# %%

# %%

# %%
