from abc import ABC, abstractmethod
import json

class AbstractImporter(ABC):

    def read_config(self):
        s = "{\n"
        with open("config.js", "r") as f:
            next(f)
            for line in f:
                s += line + "\n"
        return json.loads(s)

    def __init__(self):
        try:
            for (k, v) in self.read_config():
                setattr(self, k, v)
        except FileNotFoundError:
            raise Exception("config.js not found.")
        except JSONDecodeError:
            raise Exception("Syntax error on config.js.")

    @abstractmethod
    def import_map(self):
        pass

    @abstractmethod
    def import_sample(self):
        pass
