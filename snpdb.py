#!/usr/bin/python3

from pymongo import MongoClient
from gridfs import GridFS
import json
import os

def read_config():
    s = "{\n"
    with open("config.js", "r") as f:
        next(f)         # Skip first line as it contains an attribution.
        for line in f:
            s += line + "\n"
    return json.loads(s)




def find_snp(id=None, min_chrom=None, max_chrom=None,
             min_pos=None, max_pos=None):
    chrom = config["SNPS_CHROMOSOME_ATTR"]
    pos = config["SNPS_POSITION_ATTR"]
    name = config["SNPS_NAME_ATTR"]
    query, chrom_query, pos_query = {}, {}, {}
    
    if min_chrom is not None:
        chrom_query.update({"$gte": min_chrom})
    if max_chrom is not None:
        chrom_query.update({"$lte": max_chrom})
    if min_pos is not None:
        pos_query.update({"$gte": min_pos})
    if max_pos is not None:
        pos_query.update({"$lte": max_pos})

    if id is not None:
        query.update({name: id})
    if len(chrom_query) > 0:
        query.update({chrom: chrom_query})
    if len(pos_query) > 0:
        query.update({pos: pos_query})
        
    res = []
    for snp in snpc.find(query):
        res.append(snp)
    return res




def find_maps(id=None, min_size=None, max_size=None, format=None):
    query = {}
    size_query = {}
    if id is not None:
        query.update({"_id": id})
    if format is not None:
        query.update({config["MAPS_FORMAT_ATTR"]: format})
    if min_size is not None:
        size_query.update({"$gte": min_size})
    if max_size is not None:
        size_query.update({"$lte": max_size})
    if len(size_query) > 0:
        query.update({config["MAPS_SIZE_ATTR"]: size_query})

    res = []
    for mp in mapc.find(query, {config["MAPS_FORMAT_ATTR"]: 1, config["MAPS_SIZE_ATTR"]: 1}):
        res.append(mp)
    return res




def find_individuals(id=None, tatoo=None, sample=None):
    query = {}
    if id is not None:
        query.update({"_id": id})
    if tatoo is not None:
        query.update({config[INDIVIDUALS_ID_LIST_ATTR]: tatoo})
    if sample is not None:
        query.update({config[INDIVIDUALS_SAMPLE_LIST_ATTR]: sample})
    res = []
    for individual in indc.find(query):
        res.append(individual)
    return res


def find_snp_of_sample(mapname, sample, snp_id):
    fsnp = (config["SNPBLOCKS_SNP_LIST_ATTR"] + 
           ".0." + config["SNPBLOCKS_SNP_ID_INSIDE_LIST"])
    block = snpblocksc.find_one({
                                config["SNPBLOCKS_MAP_ATTR"]: mapname,
                                config["SNPBLOCKS_SAMPLE_ATTR"]: sample,
                                fsnp: {"$lte": snp_id}},
                                sort=[(fsnp, -1)])
    if block is None:
        return None
    try:
        res = [x for x in block[config["SNPBLOCKS_SNP_LIST_ATTR"]]
              if x[config["SNPBLOCKS_SNP_ID_INSIDE_LIST"]] == snp_id][0]
    except IndexError:
        return None
    return res




def insert_file(file, individual_id=None):
    if individual_id is None:
        gfs.put(file, filename=filename)
    else:
        gfs.put(file, filename=os.path.basename(file.name),
        **{config["FILES_INDIVIDUAL_ATTR"]: individual_id})




def list_files(individual_id=None):
    res = []
    cur = None
    if individual_id is None:
        cur = db.fs.files.find({})
    else:
        cur = db.fs.files.find({config["FILES_INDIVIDUAL_ATTR"]: individual_id})
    for f in cur:
        res.append(f)
    return res




def get_files(files):
    for file_doc in files:
        grid_out = gfs.get(file_doc["_id"])
        filename = grid_out.filename
        with open(filename, "wb+") as f:
            f.write(grid_out.read())


config = read_config()
client = MongoClient(config["HOST"])
db = client[config["DB_NAME"]]
snpc = db[config["SNPS_COLL"]]
mapc = db[config["MAPS_COLL"]]
indc = db[config["INDIVIDUALS_COLL"]]
snpblocksc = db[config["SNPBLOCKS_COLL"]]
gfs = GridFS(db)
