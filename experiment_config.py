import importlib
import os
import readers
import snpdb
import subprocess
import testing.random_file_generator as rfgen
import time
from typing import Tuple, Callable

# Constants and global variables
data_dir: str = "./data/"  # Data output directory
fastq_dir_1: str = "../1_FASTq_VCF/1_SRR10752699_MISSOURI/"  # data directory
fastq_dir_2: str = "../1_FASTq_VCF/2_VCF_bos_taurus/"  # data directory
phenotype_dir_1: str = "../2_DADOS_NOVELPHENOTYPES/1/"
phenotype_dir_2: str = "../2_DADOS_NOVELPHENOTYPES/2/"
phenotype_dir_3: str = "../2_DADOS_NOVELPHENOTYPES/3/02_csv/"
graph_dir: str = "./graphs/"  # Graph output directory
results_fname: str = data_dir + "experiment_results.json"

# Dictionary with experiment parameters
exps: dict = {
    "1.A": {
        "compression_methods": {"snappy", "zlib"},
        "file_type": "0125",
        "nsnps_list": [100000.0, 1000000.0, 50000000.0],
        "nsamples_list": [1.0],
        "readers": {"map": readers.Z125MapReader, "ped": readers.Z125SampleReader},
        "file_extensions": {"map": ".0125map", "ped": ".0125ped", "ids": ".0125ids"},
    },
    "1.B": {
        "compression_methods": {"snappy", "zlib"},
        "file_type": "PLINK",
        "nsnps_list": [100000.0, 1000000.0],
        "nsamples_list": [1.0],
        "readers": {"map": readers.PlinkMapReader, "ped": readers.PlinkSampleReader},
        "file_extensions": {"map": ".plmap", "ped": ".plped", "ids": ".plids"},
    },
    "1.C": {
        "compression_methods": {"snappy", "zlib"},
        "file_type": "FR",
        "nsnps_list": [100000.0, 1000000.0],
        "nsamples_list": [1.0],
        "readers": {
            "map": readers.FinalReportMapReader,
            "ped": readers.FinalReportSampleReader,
        },
        "file_extensions": {"ext": ".fr", "ids": ".frids"},
    },
    "1.D": {
        "compression_methods": {"snappy", "zlib"},
        "file_type": "VCF",
        "nsnps_list": [
            100000.0,
            1000000.0,
        ],
        "nsamples_list": [1.0],
        "readers": {"map": readers.VcfMapReader, "ped": readers.VcfSampleReader},
        "file_extensions": {"ext": ".vcf", "ids": ".vcfids"},
    },
    "1.E": {
        "compression_methods": {"snappy", "zlib"},
        "file_type": "FastQ",
        "nsnps_list": [-1],
        "nsamples_list": [-1],
        "file_extensions": {"ext": ".fastq"},
    },
    "1.F": {
        "compression_methods": {"snappy", "zlib"},
        "file_type": "Media",
        "nsnps_list": [-1],
        "nsamples_list": [-1],
        "file_extensions": {"ext": ".jpg"},
    },
    "1.G": {
        "compression_methods": {"snappy", "zlib"},
        "file_type": "ALL",
        "nsnps_list": [
            100000.0,
            1000000.0,
        ],
        "nsamples_list": [1.0],
        "readers": {
            "0125": {"map": readers.Z125MapReader, "ped": readers.Z125SampleReader},
            "PLINK": {"map": readers.PlinkMapReader, "ped": readers.PlinkSampleReader},
            "FR": {
                "map": readers.FinalReportMapReader,
                "ped": readers.FinalReportSampleReader,
            },
            "VCF": {"map": readers.VcfMapReader, "ped": readers.VcfSampleReader},
        },
        "file_extensions": {
            "0125": {"map": ".0125map", "ped": ".0125ped", "ids": ".0125ids"},
            "PLINK": {"map": ".plmap", "ped": ".plped", "ids": ".plids"},
            "FR": {"ext": ".fr", "ids": ".frids"},
            "VCF": {"ext": ".vcf", "ids": ".vcfids"},
            "FastQ": {"ext": ".fastq"},
            "Media": {"ext": ".jpg"},
        },
    },
    "2.A": {
        "compression_methods": {"snappy", "zlib"},
        "file_type": "0125",
        "nsnps_list": [
            100000.0,
            1000000.0,
            50000000.0,
        ],
        "nsamples_list": [
            1.0,
            5.0,
            10.0,
            50.0,
            100.0,
            500.0,
            1000.0,
            5000.0,
            10000.0,
            # 50000.0,
            # 100000.0,
            # 1000000.0,
        ],
        "readers": {"map": readers.Z125MapReader, "ped": readers.Z125SampleReader},
        "file_extensions": {"map": ".0125map", "ped": ".0125ped", "ids": ".0125ids"},
    },
    "2.B": {
        "compression_methods": {"snappy", "zlib"},
        "file_type": "VCF",
        "nsnps_list": [100000.0, 1000000.0],
        "nsamples_list": [
            1.0,
            5.0,
            10.0,
            50.0,
            100.0,
            500.0,
            1000.0,
            5000.0,
            10000.0,
        ],  # , 100000.0],  # , 1000000.0],
        "readers": {"map": readers.VcfMapReader, "ped": readers.VcfSampleReader},
        "file_extensions": {"ext": ".vcf", "ids": ".vcfids"},
    },
    "2.C": {
        "compression_methods": {"snappy", "zlib"},
        "file_type": "PLINK",
        "nsnps_list": [100000.0, 1000000.0],
        "nsamples_list": [
            1.0,
            5.0,
            10.0,
            50.0,
            100.0,
            500.0,
            1000.0,
            5000.0,
            10000.0,
            50000.0,
            100000.0,
        ],  # , 1000000.0],
        "readers": {"map": readers.PlinkMapReader, "ped": readers.PlinkSampleReader},
        "file_extensions": {"map": ".plmap", "ped": ".plped", "ids": ".plids"},
    },
    "2.D": {
        "compression_methods": {"snappy", "zlib"},
        "file_type": "ALL",
        "nsnps_list": [100000.0, 1000000.0],
        "nsamples_list": [
            1.0,
            5.0,
            10.0,
            50.0,
            100.0,
            500.0,
            1000.0,
            5000.0,
            10000.0,
        ],
        #   100000.0],  # , 1000000.0],
        "readers": {
            "0125": {"map": readers.Z125MapReader, "ped": readers.Z125SampleReader},
            "PLINK": {"map": readers.PlinkMapReader, "ped": readers.PlinkSampleReader},
            "FR": {
                "map": readers.FinalReportMapReader,
                "ped": readers.FinalReportSampleReader,
            },
            "VCF": {"map": readers.VcfMapReader, "ped": readers.VcfSampleReader},
        },
        "file_extensions": {
            "0125": {"map": ".0125map", "ped": ".0125ped", "ids": ".0125ids"},
            "PLINK": {"map": ".plmap", "ped": ".plped", "ids": ".plids"},
            "FR": {"ext": ".fr", "ids": ".frids"},
            "VCF": {"ext": ".vcf", "ids": ".vcfids"},
            "FastQ": {"ext": ".fastq"},
            "Media": {"ext": ".jpg"},
        },
    },
    "2.7": {
        "compression_methods": {"snappy", "zlib"},
        "file_type": "0125",
        "nsnps_list": [100000.0],
        "nsamples_list": [10000.0],
        "readers": {"map": readers.Z125MapReader, "ped": readers.Z125SampleReader},
        "file_extensions": {"map": ".0125map", "ped": ".0125ped", "ids": ".0125ids"},
    },
}
nsnps_ids: dict = {
    100000.0: "100k",
    1000000.0: "1m",
    50000000.0: "50m",
}
nsamples_ids: dict = {
    1.0: "1",
    5.0: "5",
    10.0: "10",
    50.0: "50",
    100.0: "100",
    500.0: "500",
    1000.0: "1k",
    5000.0: "5k",
    10000.0: "10k",
    50000.0: "50k",
    100000.0: "100k",
    500000.0: "500k",
    1000000.0: "1m",
    50000000.0: "50m",
}
bin_file_types: set = {"FastQ", "Media"}


def reset_db(compression_method: str) -> None:
    """Reset MongoDB database for experiments.

    Args:
        compression_method (str, optional): Compression method for database
            after reset operation. May be either 'snappy' or 'zlib'.
            Defaults to 'snappy'.
    """
    subprocess.run(
        ["/home/rodrigo/snpdb/reset_db.sh", compression_method],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.STDOUT,
        check=True,
    )
    importlib.reload(snpdb)


def generate_random_file(
    filename: str, file_type: str, verbose: bool = False, **kwargs
) -> Tuple[float, float]:
    """Generate random files for snpdb experiments.

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
        ".0125map": rfgen.random_0125_map,
        ".0125ped": rfgen.random_0125_samples,
        ".fr": rfgen.random_final_report,
        ".vcf": rfgen.random_vcf,
        ".plmap": rfgen.random_plink_map,
        ".plped": rfgen.random_plink_samples,
        ".ids": rfgen.id_mapping,
    }[file_type]
    t: float = time.time()
    with open(filename, "w") as f:
        func(outfile=f, **kwargs)
    t = time.time() - t
    if verbose:
        print(
            "Generated "
            + file_type
            + "\tTime: "
            + str(round(t, 3))
            + "s\tSize: "
            + str(round(os.stat(filename).st_size / (1024 ** 2), 2))
            + "MB"
        )
    return t, os.stat(filename).st_size
