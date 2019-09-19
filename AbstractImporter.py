from abc import ABC, abstractmethod
from pymongo import MongoClient, UpdateOne
from pymongo.errors import DuplicateKeyError
from pprint import pprint
import json

class AbstractImporter(ABC):

    @abstractmethod
    def read_map(self, filename):
        pass

    @abstractmethod
    def read_samples(self, filename):
        pass
  
    def import_map(self, filename, mapname,
                   force_create_new=False,
                   force_use_existing=False):
        db = self.__get_db()
        if db[self.MAPS_COLL].find_one({"_id": mapname}) is not None:
            raise Exception("Map name already in use.")
        snps = self.read_map(filename)
        first_new_id = self.__fill_snp_ids(snps, force_create_new, force_use_existing)
        
        if first_new_id is None:
            return          # User abort.
        
        db[self.MAPS_COLL].insert_one(
            {"_id": mapname,
            self.MAPS_SNP_LIST_ATTR: [snp["_id"] for snp in snps]})
        new_snps = [snp for snp in snps if snp["_id"] >= first_new_id]
        if len(new_snps) > 0:
            db[self.SNPS_COLL].insert_many(new_snps)
        db[self.SNPS_COLL].bulk_write([UpdateOne({"_id": snp["_id"]},
            {"$push": {self.SNPS_MAPS_ATTR: mapname}}) for snp in snps])
        
    def import_sample(self, filename, mapname):
        pass 

    def __read_config(self):
        s = "{\n"
        with open("config.js", "r") as f:
            next(f)         # Skip first line as it contains an attribution.
            for line in f:
                s += line + "\n"
        return json.loads(s)

    def __init__(self):
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
