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
            for (k, v) in self.read_config().items():
                setattr(self, k, v)
        except FileNotFoundError:
            raise Exception("config.js not found.") from None
        except json.decoder.JSONDecodeError:
            raise Exception("Syntax error on config.js.") from None


    @abstractmethod
    def import_map(self, filename):
        pass


    @abstractmethod
    def import_sample(self, filename):
        pass
