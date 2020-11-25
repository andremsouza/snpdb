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
    with open(results_fname, "r") as f:
        results = json.load(f)
except Exception:
    results = {}

# Storing data in a Pandas DataFrame
col_names: list = [
    "experiment_id",
    "compression_method",
    "file_type",
    "nsnps",
    "nsamples",
    "fsize",
    "dbsize",
    "time",
    "summarize",
    "individuals_of_snps",
    "delete_individual",
    "export_time_0125",
    "export_time_plink",
    "export_time_bin",
]
df: pd.DataFrame = pd.DataFrame(columns=col_names)

# %% [markdown]
# ### Abrindo e processando arquivos de resultados

# %%

# Reading experiment result files and building DataFrame
for experiment_id in exps:
    if experiment_id not in results:
        continue
    for compression_method in exps[experiment_id]["compression_methods"]:
        for nsnps in exps[experiment_id]["nsnps_list"]:
            for nsamples in exps[experiment_id]["nsamples_list"]:
                if exps[experiment_id]["file_type"] in bin_file_types:
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
                if exps[experiment_id]["file_type"] in bin_file_types:
                    fsizes = results[experiment_id][compression_method]["fsize"]
                    dbsizes = results[experiment_id][compression_method]["dbsize"]
                    times = results[experiment_id][compression_method]["time"]
                    summarize_times = [0.0] * len(fsizes)
                    individuals_of_snps_times = [0.0] * len(fsizes)
                    delete_individual_times = [0.0] * len(fsizes)
                    export_time_0125 = [0.0] * len(fsizes)
                    export_time_plink = [0.0] * len(fsizes)
                    export_time_bin = [0.0] * len(fsizes)
                elif experiment_id == "2.7":
                    fsizes = results[experiment_id][compression_method][nsnps_id][
                        nsamples_id
                    ]["fsize"]
                    dbsizes = results[experiment_id][compression_method][nsnps_id][
                        nsamples_id
                    ]["dbsize"]
                    times = results[experiment_id][compression_method][nsnps_id][
                        nsamples_id
                    ]["time"]
                    summarize_times = [0.0] * len(fsizes)
                    individuals_of_snps_times = [0.0] * len(fsizes)
                    delete_individual_times = [0.0] * len(fsizes)
                    export_time_0125 = [0.0] * len(fsizes)
                    export_time_plink = [0.0] * len(fsizes)
                    export_time_bin = [0.0] * len(fsizes)
                else:
                    fsizes = results[experiment_id][compression_method][nsnps_id][
                        nsamples_id
                    ]["fsize"]
                    dbsizes = results[experiment_id][compression_method][nsnps_id][
                        nsamples_id
                    ]["dbsize"]
                    times = results[experiment_id][compression_method][nsnps_id][
                        nsamples_id
                    ]["time"]
                    # Getting results values for summarize operation
                    summarize_times = [
                        summ["time"]
                        for summ in results[experiment_id][compression_method][
                            nsnps_id
                        ][nsamples_id]["summarize"]
                    ]
                    for ind in results[experiment_id][compression_method][nsnps_id][
                        nsamples_id
                    ]["individuals_of_snps"]:
                        try:
                            # Fix KeyError in results
                            individuals_of_snps_times.append(float(ind["time"]))
                        except KeyError:
                            individuals_of_snps_times.append(float(ind["snp"]))
                    delete_individual_times = [
                        ind["time"]
                        for ind in results[experiment_id][compression_method][nsnps_id][
                            nsamples_id
                        ]["delete_individual"]
                    ]
                    if exps[experiment_id]["file_type"] == "ALL":
                        export_time_0125 = [
                            i["Z125"]["map"] + i["Z125"]["samples"]
                            for i in results[experiment_id][compression_method][
                                nsnps_id
                            ][nsamples_id]["export"]
                        ]
                        export_time_plink = [
                            i["PLINK"]["map"] + i["PLINK"]["samples"]
                            for i in results[experiment_id][compression_method][
                                nsnps_id
                            ][nsamples_id]["export"]
                        ]
                        export_time_bin = [
                            i
                            for i in results[experiment_id][compression_method][
                                nsnps_id
                            ][nsamples_id]["export_bin"]
                        ]
                    else:
                        export_time_0125 = [0.0] * len(fsizes)
                        export_time_plink = [0.0] * len(fsizes)
                        export_time_bin = [0.0] * len(fsizes)

                # Mounting rows and appending to DataFrame
                rows = [
                    {
                        "experiment_id": experiment_id,
                        "compression_method": compression_method,
                        "file_type": exps[experiment_id]["file_type"],
                        "nsnps": nsnps,
                        "nsamples": nsamples,
                        "fsize": float(i),
                        "dbsize": float(j),
                        "time": float(k),
                        "summarize": float(l),
                        "individuals_of_snps": float(m),
                        "delete_individual": float(n),
                        "export_time_0125": float(o),
                        "export_time_plink": float(p),
                        "export_time_bin": float(q),
                    }
                    for i, j, k, l, m, n, o, p, q in zip(
                        fsizes,
                        dbsizes,
                        times,
                        summarize_times,
                        individuals_of_snps_times,
                        delete_individual_times,
                        export_time_0125,
                        export_time_plink,
                        export_time_bin,
                    )
                ]
                df = df.append(rows, ignore_index=True, verify_integrity=True)

# %%
# Filtering DataFrame to only 'zlib' experiments
df = df[df["compression_method"] == "zlib"]

# %% [markdown]
# #### Criando colunas concatenadas
df["file_type_nsnps"] = df["file_type"].map(str) + "_" + df["nsnps"].map(str)

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
df_melted = pd.melt(
    df[df["experiment_id"].str.startswith("1")],
    id_vars=["experiment_id", "compression_method", "file_type", "nsnps", "nsamples"],
    value_vars=["time"],
)
df_melted = df_melted.fillna(value="Binary file")
df_melted.loc[df_melted["nsnps"] == 1000000, "nsnps"] = nsnps_ids[float(1000000)]
df_melted.loc[df_melted["nsnps"] == 100000, "nsnps"] = nsnps_ids[float(100000)]
sns.set(style="whitegrid", palette=sns.color_palette("muted", n_colors=6, desat=1.0))
snsplot = sns.catplot(
    x="file_type",
    y="value",
    hue="nsnps",
    data=df_melted,
    # col="compression_method",
    # col_order=["zlib"],
    kind="bar",
    legend_out=False,
    height=6,
)

snsplot._legend.set_title("# SNPs")
# snsplot.set(xscale='log')
snsplot.set(yscale="log")
snsplot = snsplot.set_axis_labels("File type", "Import time (s)")
# .set_titles(template="Método de compressão: {row_name} - {col_name}")
snsplot.savefig(graph_dir + "experiment1_times.png")
plt.draw()

# %% [markdown]
# ### Comparação de tamanhos de arquivos antes/depois das inserções

# %%
df_melted = pd.melt(
    df[df["experiment_id"].str.startswith("1")],
    id_vars=["experiment_id", "compression_method", "file_type", "nsnps", "nsamples"],
    value_vars=["fsize", "dbsize"],
)
df_melted = df_melted.fillna(value="bin")
df_melted["variable"] = (
    df_melted["variable"].map(str) + "_" + df_melted["nsnps"].map(str)
)
df_melted["variable"] = (
    df_melted["variable"].str.replace("100000.0", "100k").str.replace("1000000.0", "1m")
)
df_melted["value"] /= 1024 ** 2
sns.set(style="whitegrid", palette=sns.color_palette("muted", n_colors=6, desat=1.0))
snsplot = sns.catplot(
    x="file_type",
    y="value",
    hue="variable",
    hue_order=[
        "fsize_100k",
        "dbsize_100k",
        "fsize_1m",
        "dbsize_1m",
        "fsize_bin",
        "dbsize_bin",
    ],
    data=df_melted,
    kind="bar",
    legend_out=False,
    height=6,
)

snsplot._legend.set_title("size_nsnps")
# snsplot.set(xscale='log')
snsplot.set(yscale="log")
snsplot = snsplot.set_axis_labels("File type", "Size (MiB)").set_titles(
    template="{col_var} = {col_name}",
)
snsplot.savefig(graph_dir + "experiment1_sizes.png")
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
df_melted = pd.melt(
    df[df["experiment_id"].str.startswith("2")],
    id_vars=["experiment_id", "compression_method", "file_type", "nsnps", "nsamples"],
    value_vars=["time"],
)
sns.set(style="whitegrid", palette=sns.color_palette("muted", n_colors=6, desat=1.0))
fig = plt.figure(figsize=(6.4, 6.4))
ax = sns.lineplot(
    data=df_melted,
    x="nsamples",
    y="value",
    hue="file_type",
    hue_order=["0125", "PLINK", "VCF", "ALL"],
    palette="tab10",
    style="nsnps",
    ci=None,
    legend="auto",
    markers=True,
)

ax.set_xscale("log")
ax.set_yscale("log")
# ax.set_title("compression_method = zlib")
ax.set_xlabel("# of samples")
ax.set_ylabel("Import time (s)")

# for nsnps in df_melted["nsnps"].unique():
#     for file_type in df_melted["file_type"].unique():
#         df_tmp = df_melted[df_melted["file_type"] == file_type][
#             df_melted["nsnps"] == nsnps
#         ]
#         ax.text(
#             df_tmp["nsamples"].max(),
#             df_tmp["value"].max(),
#             "%.1f" % float(df_tmp["value"].max()),
#             fontsize=8,
#             color="black",
#             ha="center",
#             va="bottom",
#         )
fig.savefig(graph_dir + "experiment2_1_times.png")
plt.draw()

# %% [markdown]
# ### Tamanhos de arquivos antes e depois da inserção (100k)

# %%
df_melted = pd.melt(
    df[df["experiment_id"].str.startswith("2")],
    id_vars=["experiment_id", "compression_method", "file_type", "nsnps", "nsamples"],
    value_vars=["fsize", "dbsize"],
)
df_melted = df_melted.fillna(value="bin")
df_melted = df_melted[df_melted["nsnps"] == 100000.0]
# df_melted['variable'] = df_melted['variable'].map(
#     str) + '_' + df_melted['nsnps'].map(str)
# df_melted['variable'] = df_melted['variable'].str.replace(
#     '100000.0', '100k').str.replace('1000000.0', '1m')
df_melted["value"] /= 1024 ** 2
sns.set(style="whitegrid", palette=sns.color_palette("muted", n_colors=6, desat=1.0))

fig = plt.figure(figsize=(6.4, 6.4))
ax = sns.lineplot(
    data=df_melted,
    x="nsamples",
    y="value",
    hue="file_type",
    hue_order=["0125", "PLINK", "VCF", "ALL"],
    palette="tab10",
    style="variable",
    style_order=["fsize", "dbsize"],
    ci=None,
    legend="auto",
    markers=True,
)

ax.set_xscale("log")
ax.set_yscale("log")
# ax.set_title("compression_method = zlib")
ax.set_xlabel("# of samples")
ax.set_ylabel("Size (MiB)")

# for nsnps in df_melted["nsnps"].unique():
#     for file_type in df_melted["file_type"].unique():
#         df_tmp = df_melted[df_melted["file_type"] == file_type][
#             df_melted["nsnps"] == nsnps
#         ]
#         ax.text(
#             df_tmp["nsamples"].max(),
#             df_tmp["value"].max(),
#             "%.1f" % float(df_tmp["value"].max()),
#             fontsize=8,
#             color="black",
#             ha="center",
#             va="bottom",
#         )
fig.savefig(graph_dir + "experiment2_1_sizes_100k.png")
plt.draw()

# %% [markdown]
# ### Tamanhos de arquivos antes e depois da inserção (1m)

# %%
df_melted = pd.melt(
    df[df["experiment_id"].str.startswith("2")],
    id_vars=["experiment_id", "compression_method", "file_type", "nsnps", "nsamples"],
    value_vars=["fsize", "dbsize"],
)
df_melted = df_melted.fillna(value="bin")
df_melted = df_melted[df_melted["nsnps"] == 1000000.0]
# df_melted['variable'] = df_melted['variable'].map(
#     str) + '_' + df_melted['nsnps'].map(str)
# df_melted['variable'] = df_melted['variable'].str.replace(
#     '100000.0', '100k').str.replace('1000000.0', '1m')
df_melted["value"] /= 1024 ** 2
sns.set(style="whitegrid", palette=sns.color_palette("muted", n_colors=6, desat=1.0))

fig = plt.figure(figsize=(6.4, 6.4))
ax = sns.lineplot(
    data=df_melted,
    x="nsamples",
    y="value",
    hue="file_type",
    hue_order=["0125", "PLINK", "VCF", "ALL"],
    palette="tab10",
    style="variable",
    style_order=["fsize", "dbsize"],
    ci=None,
    legend="auto",
    markers=True,
)

ax.set_xscale("log")
ax.set_yscale("log")
# ax.set_title("compression_method = zlib")
ax.set_xlabel("# of samples")
ax.set_ylabel("Size (MiB)")

# for nsnps in df_melted["nsnps"].unique():
#     for file_type in df_melted["file_type"].unique():
#         df_tmp = df_melted[df_melted["file_type"] == file_type][
#             df_melted["nsnps"] == nsnps
#         ]
#         g.text(
#             df_tmp["nsamples"].max(),
#             df_tmp["value"].max(),
#             "%.1f" % float(df_tmp["value"].max()),
#             fontsize=8,
#             color="black",
#             ha="center",
#             va="bottom",
#         )
fig.savefig(graph_dir + "experiment2_1_sizes_1m.png")
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
df_melted = pd.melt(
    df[df["experiment_id"].str.startswith("2")],
    id_vars=["experiment_id", "compression_method", "file_type", "nsnps", "nsamples"],
    value_vars=["summarize"],
)
sns.set(style="whitegrid", palette=sns.color_palette("muted", n_colors=6, desat=1.0))

fig = plt.figure(figsize=(6.4, 6.4))
ax = sns.lineplot(
    data=df_melted,
    x="nsamples",
    y="value",
    hue="file_type",
    hue_order=["0125", "PLINK", "VCF", "ALL"],
    palette="tab10",
    style="nsnps",
    ci=None,
    legend="auto",
    markers=True,
)

ax.set_xscale("log")
ax.set_yscale("log")
# ax.set_title("compression_method = zlib")
ax.set_xlabel("# of samples")
ax.set_ylabel("Execution time (s)")

# for nsnps in df_melted["nsnps"].unique():
#     for file_type in df_melted["file_type"].unique():
#         df_tmp = df_melted[df_melted["file_type"] == file_type][
#             df_melted["nsnps"] == nsnps
#         ]
#         ax.text(
#             df_tmp["nsamples"].max(),
#             df_tmp["value"].max(),
#             "%.1f" % float(df_tmp["value"].max()),
#             fontsize=8,
#             color="black",
#             ha="center",
#             va="bottom",
#         )
fig.savefig(graph_dir + "experiment2_summarize_times.png")
plt.draw()

# %% [markdown]
# ## Experimento 2.3 - Exportação de dados para formatos originais
# - Exportação de um par mapa/amostra de um indivíduo, escolhido
# aleatóriamente, para os formatos 0125 e PLINK.

# %%
df_melted = pd.melt(
    df[df["file_type"] == "ALL"][df["nsnps"] == 1000000.0],
    id_vars=["nsamples"],
    value_vars=["export_time_0125", "export_time_plink", "export_time_bin"],
)
df_melted.loc[:, "nsamples"] = df_melted["nsamples"].astype(np.int)
sns.set(style="whitegrid", palette=sns.color_palette("muted", n_colors=6, desat=1.0))

snsplot = sns.catplot(
    x="nsamples",
    y="value",
    hue="variable",
    data=df_melted,
    # col="experiment_id",
    kind="bar",
    ci=None,
    height=6,
)

snsplot = snsplot.set_axis_labels("# of samples", "Execution time (s)")
snsplot.savefig(graph_dir + "experiment2_export_times.png")
plt.draw()

# %% [markdown]
# ## Experimento 2.4 - Busca de indivíduos, dada uma lista de SNPs
# Dada uma lista de marcadores, deseja-se buscar todos os animais que os
# possuem, nos diferentes tipos dearquivos.

# %%
df_melted = pd.melt(
    df[df["experiment_id"].str.startswith("2")],
    id_vars=["experiment_id", "compression_method", "file_type", "nsnps", "nsamples"],
    value_vars=["individuals_of_snps"],
)
sns.set(style="whitegrid", palette=sns.color_palette("muted", n_colors=6, desat=1.0))

fig = plt.figure(figsize=(6.4, 6.4))
ax = sns.lineplot(
    data=df_melted,
    x="nsamples",
    y="value",
    hue="file_type",
    hue_order=["0125", "PLINK", "VCF", "ALL"],
    palette="tab10",
    style="nsnps",
    ci=None,
    legend="auto",
    markers=True,
)

ax.set_xscale("log")
ax.set_yscale("log")
# ax.set_title("compression_method = zlib")
ax.set_xlabel("# of samples")
ax.set_ylabel("Execution time (s)")

# for nsnps in df_melted["nsnps"].unique():
#     for file_type in df_melted["file_type"].unique():
#         df_tmp = df_melted[df_melted["file_type"] == file_type][
#             df_melted["nsnps"] == nsnps
#         ]
#         ax.text(
#             df_tmp["nsamples"].max(),
#             df_tmp["value"].max(),
#             "%.1f" % float(df_tmp["value"].max()),
#             fontsize=8,
#             color="black",
#             ha="center",
#             va="bottom",
#         )
fig.savefig(graph_dir + "experiment2_indiv_search_times.png")
plt.draw()

# %% [markdown]
# ## Experimento 2.5 - Exportação de dados brutos/binários
# - Exportação de todos os arquivos armazenados em modo "bruto" de um indivíduo
# aleatoriamente selecionado.
# - Para cada experimento, é armazena um arquivo FastQ e um arquivo de mídia.

# %%
df_melted = pd.melt(
    df[df["file_type"] == "ALL"],
    id_vars=["nsamples", "nsnps"],
    value_vars=["export_time_0125", "export_time_plink", "export_time_bin"],
)
df_melted.loc[:, "nsamples"] = df_melted["nsamples"].astype(np.int)
sns.set(style="whitegrid", palette=sns.color_palette("muted", n_colors=6, desat=1.0))

snsplot = sns.catplot(
    x="nsamples",
    y="value",
    hue="variable",
    data=df_melted,
    col="nsnps",
    kind="bar",
    ci=None,
    height=6,
)

snsplot = snsplot.set_axis_labels("# of samples", "Execution time (s)")
snsplot.savefig(graph_dir + "experiment2_export_times_bin.png")
plt.draw()

# %% [markdown]
# ## Experimento 2.6 - Remoção de todos os dados de um indivíduo
# Dada uma lista de indivíduos, remover todos os dados associados

# %%
df_melted = pd.melt(
    df[df["experiment_id"].str.startswith("2")],
    id_vars=["experiment_id", "compression_method", "file_type", "nsnps", "nsamples"],
    value_vars=["delete_individual"],
)
sns.set(style="whitegrid", palette=sns.color_palette("muted", n_colors=6, desat=1.0))

fig = plt.figure(figsize=(6.4, 6.4))
ax = sns.lineplot(
    data=df_melted,
    x="nsamples",
    y="value",
    hue="file_type",
    hue_order=["0125", "PLINK", "VCF", "ALL"],
    palette="tab10",
    style="nsnps",
    ci=None,
    legend="auto",
    markers=True,
)

ax.set_xscale("log")
ax.set_yscale("log")
# ax.set_title("compression_method = zlib")
ax.set_xlabel("# of samples")
ax.set_ylabel("Execution time (s)")

# for nsnps in df_melted["nsnps"].unique():
#     for file_type in df_melted["file_type"].unique():
#         df_tmp = df_melted[df_melted["file_type"] == file_type][
#             df_melted["nsnps"] == nsnps
#         ]
#         ax.text(
#             df_tmp["nsamples"].max(),
#             df_tmp["value"].max(),
#             "%.1f" % float(df_tmp["value"].max()),
#             fontsize=8,
#             color="black",
#             ha="center",
#             va="bottom",
#         )
plt.savefig(graph_dir + "experiment2_remove_times.png")
plt.draw()

# %% [markdown]
# ## Experimento 2.7 - Comparação de importação/exportação no modo “bruto”

# %%
# %% [markdown]
# ### Tempos de execução

# %%
df["exp_27"] = "logical"
df.loc[df["experiment_id"] == "2.7", "exp_27"] = "binary"
df_melted = pd.melt(
    df[df["experiment_id"].isin(["2.A", "2.7"])][df["nsnps"] == 100000.0][
        df["nsamples"] == 10000.0
    ],
    id_vars=["exp_27"],
    value_vars=["time"],
)
sns.set(style="whitegrid", palette=sns.color_palette("muted", n_colors=6, desat=1.0))
snsplot = sns.catplot(
    x="exp_27",
    y="value",
    data=df_melted,
    kind="bar",
    height=6,
)

# snsplot._legend.set_title("Experimento")
# snsplot.set(xscale='log')
# snsplot.set(yscale="log")
snsplot = snsplot.set_axis_labels("Ingestion mode", "Execution time (s)")
# .set_titles(template="Método de compressão: {row_name} - {col_name}")
snsplot.savefig(graph_dir + "experiment2_7_times.png")
plt.draw()

# %% [markdown]
# ### Comparação de tamanhos de arquivos antes/depois das inserções (100k SNPs)

# %%
df["exp_27"] = "logical"
df.loc[df["experiment_id"] == "2.7", "exp_27"] = "binary"
df_melted = pd.melt(
    df[df["experiment_id"].isin(["2.A", "2.7"])][df["nsnps"] == 100000.0][
        df["nsamples"] == 10000.0
    ],
    id_vars=["exp_27"],
    value_vars=["fsize", "dbsize"],
)
df_melted["value"] /= 1024 ** 2
sns.set(style="whitegrid", palette=sns.color_palette("muted", n_colors=6, desat=1.0))
snsplot = sns.catplot(
    x="exp_27",
    y="value",
    hue="variable",
    data=df_melted,
    # col="experiment_id",
    kind="bar",
    height=6,
)

snsplot._legend.set_title("")
# snsplot.set(xscale='log')
# snsplot.set(yscale="log")
snsplot = snsplot.set_axis_labels("Ingestion mode", "Size (MiB)")
snsplot.savefig(graph_dir + "experiment2_7_sizes.png")
plt.draw()

# %%
