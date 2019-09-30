from abc import ABC, abstractmethod

class SampleReader(ABC):
    SAMPLE_ID = "sample"
    SAMPLE_GENOTYPE = "genotype"

    def __init__(self, ped_file_path):
        self._PED_FILE = ped_file_path

    @abstractmethod
    def __iter__(self):
        pass

    @abstractmethod
    def __len__(self):
        pass
