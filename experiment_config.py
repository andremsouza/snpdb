import readers

# Constants and global variables
data_dir: str = './data/'  # Data output directory
fastq_dir_1: str = '../1_FASTq_VCF/1_SRR10752699_MISSOURI/'  # data directory
fastq_dir_2: str = '../1_FASTq_VCF/2_VCF_bos_taurus/'  # data directory
graph_dir: str = './graphs/'  # Graph output directory
results_fname: str = data_dir + 'experiment_results.json'

# Dictionary with experiment parameters
exps: dict = {
    '1.A': {
        'compression_methods': {'snappy', 'zlib'},
        'file_type': '0125',
        'nsnps_list': [100000.0, 1000000.0],
        'nsamples_list': [1.0],
        'readers': {
            'map': readers.Z125MapReader,
            'ped': readers.Z125SampleReader
        },
        'file_extensions': {
            'map': '.0125map',
            'ped': '.0125ped',
            'ids': '.0125ids'
        },
    },
    '1.B': {
        'compression_methods': {'snappy', 'zlib'},
        'file_type': 'PLINK',
        'nsnps_list': [100000.0, 1000000.0],
        'nsamples_list': [1.0],
        'readers': {
            'map': readers.PlinkMapReader,
            'ped': readers.PlinkSampleReader
        },
        'file_extensions': {
            'map': '.plmap',
            'ped': '.plped',
            'ids': '.plids'
        },
    },
    '1.C': {
        'compression_methods': {'snappy', 'zlib'},
        'file_type': 'FR',
        'nsnps_list': [100000.0, 1000000.0],
        'nsamples_list': [1.0],
        'readers': {
            'map': readers.FinalReportMapReader,
            'ped': readers.FinalReportSampleReader
        },
        'file_extensions': {
            'ext': '.fr',
            'ids': '.frids'
        },
    },
    '1.D': {
        'compression_methods': {'snappy', 'zlib'},
        'file_type': 'VCF',
        'nsnps_list': [
            100000.0,
            1000000.0,
        ],
        'nsamples_list': [1.0],
        'readers': {
            'map': readers.VcfMapReader,
            'ped': readers.VcfSampleReader
        },
        'file_extensions': {
            'ext': '.vcf',
            'ids': '.vcfids'
        },
    },
    '1.E': {
        'compression_methods': {'snappy', 'zlib'},
        'file_type': 'FastQ',
        'nsnps_list': [-1],
        'nsamples_list': [-1],
        'file_extensions': {
            'ext': '.fastq'
        },
    },
    '1.F': {
        'compression_methods': {'snappy', 'zlib'},
        'file_type': 'Media',
        'nsnps_list': [-1],
        'nsamples_list': [-1],
        'file_extensions': {
            'ext': '.jpg'
        }
    },
    '1.G': {
        'compression_methods': {'snappy', 'zlib'},
        'file_type': 'ALL',
        'nsnps_list': [
            100000.0,
            1000000.0,
        ],
        'nsamples_list': [1.0],
        'readers': {
            '0125': {
                'map': readers.Z125MapReader,
                'ped': readers.Z125SampleReader
            },
            'PLINK': {
                'map': readers.PlinkMapReader,
                'ped': readers.PlinkSampleReader
            },
            'FR': {
                'map': readers.FinalReportMapReader,
                'ped': readers.FinalReportSampleReader
            },
            'VCF': {
                'map': readers.VcfMapReader,
                'ped': readers.VcfSampleReader
            }
        },
        'file_extensions': {
            '0125': {
                'map': '.0125map',
                'ped': '.0125ped',
                'ids': '.0125ids'
            },
            'PLINK': {
                'map': '.plmap',
                'ped': '.plped',
                'ids': '.plids'
            },
            'FR': {
                'ext': '.fr',
                'ids': '.frids'
            },
            'VCF': {
                'ext': '.vcf',
                'ids': '.vcfids'
            },
            'FastQ': {
                'ext': '.fastq'
            },
            'Media': {
                'ext': '.jpg'
            }
        }
    },
    '2.A': {
        'compression_methods': {'snappy', 'zlib'},
        'file_type':
        '0125',
        'nsnps_list': [100000.0, 1000000.0],
        'nsamples_list': [
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
        'readers': {
            'map': readers.Z125MapReader,
            'ped': readers.Z125SampleReader
        },
        'file_extensions': {
            'map': '.0125map',
            'ped': '.0125ped',
            'ids': '.0125ids'
        },
    },
    '2.B': {
        'compression_methods': {'snappy', 'zlib'},
        'file_type':
        'VCF',
        'nsnps_list': [100000.0, 1000000.0],
        'nsamples_list': [
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
        'readers': {
            'map': readers.VcfMapReader,
            'ped': readers.VcfSampleReader
        },
        'file_extensions': {
            'ext': '.vcf',
            'ids': '.vcfids'
        },
    },
    '2.C': {
        'compression_methods': {'snappy', 'zlib'},
        'file_type':
        'ALL',
        'nsnps_list': [100000.0, 1000000.0],
        'nsamples_list': [
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
        'readers': {
            '0125': {
                'map': readers.Z125MapReader,
                'ped': readers.Z125SampleReader
            },
            'PLINK': {
                'map': readers.PlinkMapReader,
                'ped': readers.PlinkSampleReader
            },
            'FR': {
                'map': readers.FinalReportMapReader,
                'ped': readers.FinalReportSampleReader
            },
            'VCF': {
                'map': readers.VcfMapReader,
                'ped': readers.VcfSampleReader
            }
        },
        'file_extensions': {
            '0125': {
                'map': '.0125map',
                'ped': '.0125ped',
                'ids': '.0125ids'
            },
            'PLINK': {
                'map': '.plmap',
                'ped': '.plped',
                'ids': '.plids'
            },
            'FR': {
                'ext': '.fr',
                'ids': '.frids'
            },
            'VCF': {
                'ext': '.vcf',
                'ids': '.vcfids'
            },
            'FastQ': {
                'ext': '.fastq'
            },
            'Media': {
                'ext': '.jpg'
            }
        }
    },
}
nsnps_ids: dict = {100000.0: '100k', 1000000.0: '1m'}
nsamples_ids: dict = {
    1.0: '1',
    5.0: '5',
    10.0: '10',
    50.0: '50',
    100.0: '100',
    500.0: '500',
    1000.0: '1k',
    5000.0: '5k',
    10000.0: '10k',
    50000.0: '50k',
    100000.0: '100k',
    500000.0: '500k',
    1000000.0: '1m',
}
bin_file_types: set = {'FastQ', 'Media', 'ALL'}
