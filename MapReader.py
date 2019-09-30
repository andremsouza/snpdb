from abc import ABC, abstractmethod

class MapReader(ABC): 
    SNP_NAME = "name"
    SNP_CHROM = "chr"
    SNP_POS = "pos"
    MAP_FORMAT = "format"

    def __init__(self, map_file_path):
        self._MAP_FILE = map_file_path

    @abstractmethod
    def __iter__(self):
        pass

    @abstractmethod
    def __len__(self):
        pass

    @abstractmethod
    def map_meta(self):
        pass

