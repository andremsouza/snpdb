from abc import ABC, abstractmethod
from pymongo import MongoClient, UpdateOne
from pymongo.errors import DuplicateKeyError
from pprint import pprint
import json

class AbstractImporter(ABC):
    """Subclasses of this class are able to import files (maps and samples)
    into the database.
    
    In order to configure itself, this class expects a "config.js" JSON file
    on the current working directory.
    The file is read on initialization and the attributes contained are copied
    to the class as instance attributes.

    Public methods:
        import_map: imports a map from a file into the database.
        import_samples: imports samples from a file into the database.
    
    Public abstract methods:
        read_map: extracts a map from a file.
        read_samples: extract all samples of a file.
    """

    SAMPLE_DICT_ID = "sid"
    SAMPLE_DICT_GENOTYPE = "gen"

    @abstractmethod
    def read_map(self, filename):
        """Extracts a map from file "filename".
    
        It should return a tuple (snps, meta), with `snps` being a list
        of SNPs and `meta` being a dict with map metadata to be included.
        Each SNP should be a dict which may contain the following keys:
            self.SNPS_NAME_ATTR: contains the name of the SNP, as given
                                 by the file.
            self.SNPS_CHROMOSOME_ATTR: contains the SNP's chromosome,
                                       if available.
            self.SNPS_POSITION_ATTR: contains the SNP's position,
                                     if available.
        These are the keys used by the database to index and organize itself.
        Any other key will have no special meaning and will be imported to
        the database preserving its name. You should not use "_id" as key.

        The map will be imported preserving the order given in the file.
        
        As for the `meta` dict, the only special keys is self.MAPS_FORMAT_ATTR,
        which may be used to indicate the format of the file being imported,
        and self.MAPS_SIZE_ATTR, which is automatically set to the map's
        number of SNPs.
        """
        pass

    @abstractmethod
    def read_samples(self, filename):
        """Extracts a list of samples from file "filename".
        
        It should return a list of samples.
        Each sample should be a dict containing at least the following special
        keys:
            self.SAMPLE_DICT_ID: id of the sample.
            self.SAMPLE_DICT_GENOTYPE: list of dicts with genotype data.
                                       Each dict should contain the key
                                       self.SNPBLOCKS_SNP_GENOTYPE_INSIDE_LIST,
                                       whose value is the genotype data.

        Any other key will have no special meaning and will be imported to the
        database preserving its name. You should not use "_id" as key.
        
        For each sample, the order of the SNPs in the list should be the same
        as the order on the map the samples are to be associated with.
        """
        pass
    
    #TODO: e se quisesse por metadados no mapa?
    def import_map(self, filename, mapname,
                   force_create_new=False,
                   force_use_existing=False):
        """Imports a map from file to the database. Calls read_map.
        
        Attributes:
            filename: name / path of the file to import.
            mapname: how the map will be named inside the database.
            force_create_new=False: if False, for each of the map's SNPs,
                                    looks for similar SNPs inside the database.
                                    If at least one is found, asks the user
                                    whether it wants to create a new SNP or
                                    use an existing one. If none is found,
                                    creates a new SNP on the database silently.
                                    If True, skips similarity checks entirely
                                    and always creates new SNPs.
                                    Turning this option on speeds up the
                                    importing process dramatically.
            force_use_existing=False: if True, when exactly one similar SNP is
                                      found, it'll be automatically used.
                                      If none is found, a new SNP is created.
                                      If more than one is found, the user is
                                      asked to choose which to use.
                                      If False, when one similar SNP is found,
                                      the user is asked whether to use it or
                                      create a new one. This attribute and
                                      the previous one should not be used
                                      simultaneosly.
        """
        db = self.__get_db()
        if self.__get_map(mapname) is not None:
            raise Exception("Map name already in use.")
        (snps, meta) = self.read_map(filename)
        first_new_id = self.__fill_snp_ids(snps, force_create_new, force_use_existing)
        
        if first_new_id is None:
            return          # User abort.
        
        map_doc = {"_id": mapname,
                  self.MAPS_SNP_LIST_ATTR: [snp["_id"] for snp in snps],
                  self.MAPS_SIZE_ATTR: len(snps)}
        if meta is not None:
            map_doc.update(meta)
        db[self.MAPS_COLL].insert_one(map_doc)

        new_snps = [snp for snp in snps if snp["_id"] >= first_new_id]
        if len(new_snps) > 0:
            db[self.SNPS_COLL].insert_many(new_snps)
        db[self.SNPS_COLL].bulk_write([UpdateOne({"_id": snp["_id"]},
            {"$push": {self.SNPS_MAPS_ATTR: mapname}}) for snp in snps])
        
    def import_samples(self, filename, mapname):
        """Imports samples from a file to the database. Calls read_samples.

        Attributes:
            filename: name/path of the file to import.
            mapname: name of the map to use (must exist on the database).
        
        read_samples should return samples with the genotype data in the
        same order as the map on the database.
        """
        m = self.__get_map(mapname)
        if m is None:
            raise Exception("Map not found.")
        snps = m[self.MAPS_SNP_LIST_ATTR]
        samples = self.read_samples(filename)
        db = self.__get_db()
        for sample in samples:
            genotype = sample.pop(self.SAMPLE_DICT_GENOTYPE)
            id = sample.pop(self.SAMPLE_DICT_ID)
            sample.update({
                          "_id":{self.SAMPLES_MAP_ATTR: mapname,
                          self.SAMPLES_ID_ATTR: id}})
            if len(genotype) != len(snps):
                raise Exception("Sample genotype and map size mismatch.")
            db[self.SAMPLES_COLL].insert_one(sample)
            for i in range(len(genotype)):
                genotype[i][self.SNPBLOCKS_SNP_ID_INSIDE_LIST] = snps[i]
            genotype.sort(key=lambda s : s[self.SNPBLOCKS_SNP_ID_INSIDE_LIST])
            for i in range(0, len(genotype), self.SNPBLOCKS_SNPS_PER_BLOCK):
                db[self.SNPBLOCKS_COLL].insert_one({
                    self.SNPBLOCKS_MAP_ATTR: mapname,
                    self.SNPBLOCKS_SAMPLE_ATTR: id,
                    self.SNPBLOCKS_SNP_LIST_ATTR: genotype[i:i+self.SNPBLOCKS_SNPS_PER_BLOCK]}) 

    def __read_config(self):
        s = "{\n"
        with open("config.js", "r") as f:
            next(f)         # Skip first line as it contains an attribution.
            for line in f:
                s += line + "\n"
        return json.loads(s)

    def __init__(self):
        """Reads a JSON inside the config.js file on the working directory
        and initializes instance attributes with its properties.
        
        config.js is almost a JSON file. It should follow the format:
            let config = {
                JSON OBJECT
            }
        due to technical reasons.
        """
        try:
            for (k, v) in self.__read_config().items():
                setattr(self, k, v)
            self.__client = MongoClient(self.HOST)
        except FileNotFoundError:
            raise Exception("config.js not found.") from None
        except json.decoder.JSONDecodeError:
            raise Exception("Syntax error on config.js.") from None
        
    def __get_db(self):
        return self.__client[self.DB_NAME]

    def __get_map(self, mapname):
       db = self.__get_db()
       return db[self.MAPS_COLL].find_one({"_id": mapname})
    
    def __reserve_snp_ids(self, cnt):
        db = self.__get_db()
        doc = db[self.COUNTERS_COLL].find_one_and_update(
            {"_id": self.SNPS_COLL},
            {"$inc": {self.COUNTERS_SEQ_VALUE_ATTR: cnt}})
        return doc[self.COUNTERS_SEQ_VALUE_ATTR]
   
    def __find_similar_snps(self, snp):
        db = self.__get_db()
        query = {}
        for attr in [self.SNPS_NAME_ATTR,
            self.SNPS_CHROMOSOME_ATTR, self.SNPS_POSITION_ATTR]:
            if attr in snp:
                query[attr] = snp[attr]
        matches = []
        for result in db[self.SNPS_COLL].find(query):
            matches.append(result)
        return matches

    def __fill_snp_ids(self, snps,
                       force_create_new,
                       force_use_existing):
        if force_create_new and force_use_existing:
            raise Exception("force_create_new and force_use_existing cannot"
                            " be used simultaneously.") 
        db = self.__get_db()
        next_id = self.__reserve_snp_ids(len(snps))
        first_new_id = next_id
        for snp in snps:
            if force_create_new:
                snp["_id"] = next_id
                next_id += 1
                continue
            similar = self.__find_similar_snps(snp)
            if len(similar) == 0:
                snp["_id"] = next_id
                next_id += 1
            elif len(similar) == 1 and force_use_existing:
                snp["_id"] = similar[0]["_id"]
            else:
                user_choice = self.__user_snp_choice(snp, similar,
                                                     force_use_existing)
                if user_choice is None:
                    return None
                if user_choice in {"e", "E"}:
                    snp["_id"] = similar[0]["_id"]
                    force_use_existing = True
                elif user_choice in {"0", "n", "N"}:
                    snp["_id"] = next_id
                    next_id += 1
                    force_create_new = user_choice in {"n", "N"}
                else:
                    snp["_id"] = similar[int(user_choice)]["_id"]
        return first_new_id

    def __user_snp_choice(self, snp, conflicts, force_use_existing):
        id = -1
        print(str(snp) + " is similar to the following database SNPs:")
        i = 1
        for conflict_snp in conflicts:
            print("(%d)" % i, end=' ')
            pprint(conflict_snp)
            i += 1
        options = set()
        question = "What do to? "
        if not force_use_existing:
            question += "[0] - create new SNP; "
            options.add("0")
        question += "[1] to [%d] - use i-th SNP; " % (i-1)
        for j in range(1, i):
            options.add(str(j))
        if not force_use_existing:
            question += "[n] - always create new SNP"
            options.add("n")
            options.add("N")
            if len(conflicts) == 1:
                question += ("; [e] - always use existing SNP"
                           " (you'll still be asked when there's" 
                            " more than one option)")
                options.add("e")
                options.add("E")
        question += ": "
        resp = None
        while resp not in options:
            try:
                resp = input(question)
            except KeyboardInterrupt:
                return None
            if resp not in options:
                print("Invalid response.")
        return resp
