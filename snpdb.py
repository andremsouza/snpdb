"""SNP database manipulation tools.

This module contains methods to import to and query the SNP database.
When it's loaded, reads database configuration file with a call to read_config.
It expects it to be named "config.js" and exist in the working directory.
"""

from abc import ABC, abstractmethod
from pymongo import MongoClient, UpdateOne
from gridfs import GridFS
import json
import os


class MapReader(ABC):
    """Base class for implementing a map (panel) reader for a file format.

    Classes derived from this one are used by import_map in order to import
    a map file's data to the database. The  __iter__ method should return
    an iterator which, for each iteration, returns a dict instance containing
    SNP information under the appropriate keys, whose name are stored on the
    first three attributes below. There should be no need to alter the values
    of such attributes. Other keys will be imported directly to the database
    with no special meaning (beware of conflicts, see names used for the SNPs
    collection on config.js).

    Attributes
    ----------
    SNP_NAME    The key which will contain the SNP's name.
    SNP_CHROM   The key which will contain the SNP's chromosome.
    SNP_POS     The key which will contain the SNP's position.

    MAP_FORMAT  The key which will contain the map's format identifier
                on the dict returned by map_meta().
    """
    
    SNP_NAME = "name"
    SNP_CHROM = "chr"
    SNP_POS = "pos"
    MAP_FORMAT = "format"   #TODO: make this come from config

    def __init__(self, map_file):
        """Initializes the class to in order extract SNPs from given file.
        
        Attributes
        ----------
        map_file   Path to the file containg map data to be read.
        """
        self._MAP_FILE = map_file

    @abstractmethod
    def __iter__(self):
        """Returns an iterator that reads and extracts SNPs from _MAP_FILE.
        
        Abstract method. Must return an iterator which reads the file
        located in _MAP_FILE and, for each iteration, returns a dict instance
        containing a SNP's information under the appropriate keys, whose name
        are stored on the SNP_* attrs. It's ideal, though not mandatory, that
        the SNPs are returned in the same order they appear in the file.
        
        For best performance, the iterator should be a generator that
        yields SNPs while scanning the file progressively without loading it
        entirely into memory.
        """
        pass

    @abstractmethod
    def __len__(self):
        """Returns the number of SNPs contained in _MAP_FILE.
        
        Abstract method. Must return number of readable SNPs from _MAP_FILE.
        
        For best performance, it should read _MAP_FILE once on first call
        then store the SNP count. On subsequent calls, just return the count.
        """
        pass

    @abstractmethod
    def map_meta(self):
        """Returns a dict containing map metadata for _MAP_FILE.
        
        The only key with special significance is the one whose name
        is contained in MAP_FORMAT. It should contain an identifier
        of the file format.
        """
        pass




class SampleReader(ABC):
    """Base class for implementing a sample reader for a file format.

    Classes derived from this one are used by import_samples in order to
    import a sample/ped file's data to the database. The  __iter__ method
    should return an iterator which, for each iteration, returns a dict
    instance containing sample information under the appropriate keys, 
    whose name are stored on the attributes below (other keys will be ignored).
    There should be no need to alter the values of such attributes.

    Attributes
    ----------
    SAMPLE_ID           The key which will contain the sample's identifier
                        within the file.
    SAMPLE_GENOTYPE     The key which will contain the sample's genotype data
                        (see below).
    
    The sample's genotype data itself should be a dict. It should contain
    one or more arrays with the genotype data. The length M of each array
    should be the same as the size of the associated map/panel.
    In case the data consist of only one character per SNP, a string of
    length M may be used instead of the array.
    """

    SAMPLE_ID = "sample"
    SAMPLE_GENOTYPE = "genotype"

    def __init__(self, ped_file):
        """Initializes the class to in order extract samples from given file.
        
        Attributes
        ----------
        ped_file   Path to the file containg sample data to be read.
        """
        self._PED_FILE = ped_file

    @abstractmethod
    def __iter__(self):
        """Returns an iterator that reads and extracts samples from _PED_FILE.
        
        Abstract method. Must return an iterator which reads the file
        located in _PED_FILE and, for each iteration, returns a dict instance
        containing a sample's information under the appropriate keys, whose name
        are stored on the SAMPLE_* attrs.
        
        For best performance, the iterator should be a generator that
        yields samples while scanning the file progressively without loading it
        entirely into memory.
        """
        pass

    @abstractmethod
    def __len__(self):
        """Returns the number of samples contained in _PED_FILE.
        
        Abstract method. Must return number of readable samples from
        _PED_FILE.
        
        For best performance, it should read _PED_FILE once on first call
        then store the sample count. On subsequent calls, just return the count.
        """
        pass




class MapWriter(ABC):
    """Base class for implementing a map/panel writer for a file format.
    
    Classes derived from this one are used by export_map in order to
    write a map data from the database to a file. The class
    expects data of the maps SNPs, received through __init__. The write method
    takes this data and writes to a file following the appropriate format
    rules. SNP data consists of a list of dictionaries, one for each SNP,
    following the same rules as MapReader.
    
    Attributes
    ----------
    
    SNP_NAME    The dict key which will contain the SNP's name.
    SNP_CHROM   The dict key which will contain the SNP's chromosome.
    SNP_POS     The dict key which will contain the SNP's position.
    
    Any keys imported by a MapReader will also be available here.
    """
    
    SNP_NAME = MapReader.SNP_NAME
    SNP_CHROM = MapReader.SNP_CHROM
    SNP_POS = MapReader.SNP_POS
    
    
    def __init__(self, snps):
        """Initializes the writer with the appropriate SNP data.
        
        Attributes
        ----------
        snps    list (or any iterable) of dicts, each one representing a SNP to be exported.
        """
        self._snps = snps


    @abstractmethod
    def write(self, out_file_path):
        """Write received SNP data to file, preserving order. Abstract method.
        
        Attributes
        ----------
        out_file_path   Absolute or relative path of the output file.
        """
        pass




class SampleWriter(ABC):
    """Base class for implementing a sample writer for a file format.
    
    Classes derived from this one are used by export_samples in order to
    write sample/ped data from the database to a file. The class
    expects to receive data of the samples received through __init__.
    The write method takes this data and writes to a file following the appropriate
    format rules. The data received is a list of dictionaries, each one representing
    a sample, following the same rules as SampleReader.
    
    Attributes
    ----------
    
    SAMPLE_ID           The dict key which will contain the sample's
                        within file id.
    SAMPLE_GENOTYPE     The dict key which will contain the sample's
                        genotype data (see SampleReader doc for more info).
    """
    
    SAMPLE_ID = SampleReader.SAMPLE_ID
    SAMPLE_GENOTYPE = SampleReader.SAMPLE_GENOTYPE


    def __init__(self, samples):
        """Initializes the writer with the appropriate sample data.
        
        Attributes
        ----------
        samples    list (or any iterable) of dicts, each one representing a sample to be exported.
        """
        self._samples = samples


    @abstractmethod
    def write(self, out_file_path):
        """Write received sample data to file, preserving order. Abstract method.
        
        Attributes
        ----------
        out_file_path   Absolute or relative path of the output file.
        """
        pass




def read_config(path="config.js"):
    """Load configuration file and return its contents.
    
    The configuration file a JSON file, except that it begins with an
    assignment: config = {JSON...}
    This is so it can also be used by mongo_setup.js, which runs on the
    MongoDB shell.
    Returns a dict with the contents of the JSON.
    
    Attributes
    ----------
    path="config.js"    Path to the configuration file.
    """
    s = ""
    with open(path, "r") as f:
        jsonBegin = False 
        for line in f:
            if jsonBegin:
                s += line
            else:
                open_bracket = line.find("{")
                if open_bracket != -1:
                    jsonBegin = True
                    s += line[open_bracket:]
    return json.loads(s)


# Initialization.
_config = read_config()
_client = MongoClient(_config["HOST"], w=1)
_db = _client[_config["DB_NAME"]]
_SNPS = _db[_config["SNPS_COLL"]]
_MAPS = _db[_config["MAPS_COLL"]]
_INDS = _db[_config["INDIVIDUALS_COLL"]]
_SNPBLOCKS = _db[_config["SNPBLOCKS_COLL"]]
_COUNTERS = _db[_config["COUNTERS_COLL"]]
_SAMPLES = _db[_config["SAMPLES_COLL"]]
_MAPSNPS = _db[_config["MAPSNPS_COLL"]]
_GFS = GridFS(_db)


def find_snp(id=None, min_chrom=None, max_chrom=None,
             min_pos=None, max_pos=None, map=None, iid=None,
             chr=None):
    """Search SNPs in the database.

    Returns a list of dicts, each representing a SNP.

    Parameters
    ----------
    id=None          Match only SNPs with given name.
    min_chrom=None   Match only SNPs with a chromosome that compares greater
                     or equal to the one given.
    max_chrom=None   Match only SNPs with a chromosome that compares lesser
                     or equal to the one given.
    min_pos=None     Match only SNPs with a position value that compares 
                     greater or equal than the one given.
    max_pos=None     Match only SNPs with a position value that compares
                     smaller or equal than the one given.
    map=None         Match only SNPs that are contained within the map given.
    iid=None         Match only the SNP with the given internal numeric id.
    chr=None         Match only the SNPs with the given chromosome.
    """
    chrom = _config["SNPS_CHROMOSOME_ATTR"]
    pos = _config["SNPS_POSITION_ATTR"]
    name = _config["SNPS_NAME_ATTR"]
    mp = _config["SNPS_MAPS_ATTR"]
    query, chrom_query, pos_query = {}, {}, {}
    
    if min_chrom is not None:
        chrom_query.update({"$gte": min_chrom})
    if max_chrom is not None:
        chrom_query.update({"$lte": max_chrom})
    if chr is not None:
        try:
            chrom_query.update({"$eq": int(chr)})
        except ValueError:
            chrom_query.update({"$eq": chr})
    if min_pos is not None:
        pos_query.update({"$gte": min_pos})
    if max_pos is not None:
        pos_query.update({"$lte": max_pos})
    
    if id is not None:
        query.update({name: id})
    if map is not None:
        query.update({mp: map})
    if iid is not None:
        query.update({"_id": iid})
    if len(chrom_query) > 0:
        query.update({chrom: chrom_query})
    if len(pos_query) > 0:
        query.update({pos: pos_query})
    return list(_SNPS.find(query))




def find_maps(id=None, min_size=None, max_size=None, format=None):
    """Search maps in the database.

    Returns a list of dicts, each representing a map.

    Parameters
    ----------
    id=None          Match only maps with given name.
    min_size=None    Match only maps with size at least min_size.
    max_size=None    Match only maps with size at most map_size.
    format=None      Match only maps with the specified format identifier.
    """
    query = {}
    size_query = {}
    if id is not None:
        query.update({"_id": id})
    if format is not None:
        query.update({_config["MAPS_FORMAT_ATTR"]: format})
    if min_size is not None:
        size_query.update({"$gte": min_size})
    if max_size is not None:
        size_query.update({"$lte": max_size})
    if len(size_query) > 0:
        query.update({_config["MAPS_SIZE_ATTR"]: size_query})

    return list(_MAPS.find(query, {_config["MAPS_FORMAT_ATTR"]: 1, 
                                   _config["MAPS_SIZE_ATTR"]: 1,
                                   _config["MAPS_BLOCK_SIZE_ATTR"]: 1}))





def get_map_snps(id):
    """Retrive the SNPs a map contains in the form of list of internal IDs.

    A tuple of two lists is returned, first one with the IDs sorted in the
    map's original import order, and the second sorted by id value.

    Parameters
    ----------
    id      The map's id.
    """
    cur = _MAPSNPS.find({_config["MAPSNPS_MAP_ATTR"]:id},
                        sort=[(_config["MAPSNPS_IDX_ATTR"], 1)])
    snps, ssnps = [], []
    for doc in cur:
        snps.extend(doc[_config["MAPSNPS_LIST_ATTR"]])
        ssnps.extend(doc[_config["MAPSNPS_SORTED_LIST_ATTR"]])
    return snps, ssnps




def find_individuals(id=None, tatoo=None, sample_map=None, sample_id=None):
    """Search individuals in the database.

    Returns a list of dicts, one for each individual in the result.

    Parameters
    ----------
    id=None             Match only the individual with given internal id.
    tatoo=None          Match only individuals which contain tatoo among
                        their alternate IDs.
    sample_map=None     Match only individuals that have data on the specified map.
    sample_id=None      Match only individuals that have the specified sample id
                        under some map.
    """
    query = {}
    if id is not None:
        query.update({"_id": id})
    if tatoo is not None:
        query.update({_config["INDIVIDUALS_ID_LIST_ATTR"]: tatoo})
    if sample_map is not None:
        attr = (_config["INDIVIDUALS_SAMPLE_LIST_ATTR"] + "."
                + _config["SAMPLES_MAP_ATTR"])
        query.update({attr: sample_map})
    if sample_id is not None:
        attr = (_config["INDIVIDUALS_SAMPLE_LIST_ATTR"] + "."
                + _config["SAMPLES_ID_ATTR"])
        query.update({attr: sample_id})
    return list(_INDS.find(query))




def find_snp_of_sample(mapname, sample, snp_id):
    """Fetch SNP data of a given sample from a given map.

    Returns a dict with all data available.
    Parameters
    ----------
    mapname     Map to use.
    sample      Sample id within map.
    snp_id      Internal id of the SNP to fetch.
    """
    GEN = _config["SNPBLOCKS_GENOTYPE"]
    BLOCK_SIZE = _config["MAPS_BLOCK_SIZE_ATTR"]
    SORTED_SNPS = _config["MAPSNPS_SORTED_LIST_ATTR"]
    IDX = _config["MAPSNPS_IDX_ATTR"]
    MAX_LIST_SIZE = _config["MAPSNPS_MAX_LIST_SIZE"]
    MAP = _config["MAPSNPS_MAP_ATTR"]
    try:
        map = find_maps(id=mapname)[0]
        pipeline = [{"$match": {MAP: mapname}},
                    {"$project": {"idx": {"$indexOfArray": ["$" + SORTED_SNPS,
                                                            snp_id]},
                                  IDX: 1}}]
        for part in _MAPSNPS.aggregate(pipeline):
            if part["idx"] != -1:
                index = (part["idx"] + part[IDX] * MAX_LIST_SIZE)
                blk = index // map[BLOCK_SIZE]
                pos = index % map[BLOCK_SIZE]
                break

        block = _SNPBLOCKS.find_one({_config["SNPBLOCKS_MAP_ATTR"]: mapname,
                                     _config["SNPBLOCKS_BLOCK_NUMBER"]: blk,
                                     _config["SNPBLOCKS_SAMPLE_ATTR"]: sample})
         
    except (IndexError, ValueError, UnboundLocalError):
        return None
    

    if block is None:
        return None
    res = {}
    for key in block[GEN]:
        if block[GEN][key][0] == " ":
            res[key] = block[GEN][key].split()[pos]
        else:
            res[key] = block[GEN][key][pos]
    return res




def find_sample(id=None, map=None):
    """Search samples in the database.

    Parameters
    ----------
    id=None     Match samples with the specified within-map id for some map.
    map=None    Match samples which contain data under the specified map.
    """
    query = {}
    if id is not None:
        query[_config["SAMPLES_ID_ATTR"]] = id
    if map is not None:
        query[_config["SAMPLES_MAP_ATTR"]] = map
    return list(_SAMPLES.find(query))




def get_sample_data(id, map):
    """Retrive sample data.
       
       The data is returned as a dict, following the same format produced by
       the SampleWriter used to import it.

       Parameters
       ----------
       id   The sample within-map id.
       map  The sample's associated map.
    """
    samples = find_sample(id, map)
    if len(samples) > 1:
        raise Exception("Homonymous samples within the same map.")
    if len(samples) == 0:
        return None
    sample = samples[0]

    maps = find_maps(id=map)
    if len(maps) == 0:
        raise Exception("Sample map data is missing.")
    if len(maps) > 1:
        raise Exception("Homonymous maps with the same ID.")
    m = maps[0]
    snps, sorted_snps = get_map_snps(map)

    SNPBLOCKS_MAP = _config["SNPBLOCKS_MAP_ATTR"]
    SNPBLOCKS_SAMPLE = _config["SNPBLOCKS_SAMPLE_ATTR"]
    SNPBLOCKS_NO = _config["SNPBLOCKS_BLOCK_NUMBER"]
    blocks = _SNPBLOCKS.find({SNPBLOCKS_MAP: map,
                              SNPBLOCKS_SAMPLE: id},
                              sort=[(SNPBLOCKS_NO, 1)])
    genotype = {}
    for block in blocks:
        g = block[_config["SNPBLOCKS_GENOTYPE"]]
        for key in g:
            if g[key][0] == " ":
                data = g[key].split()
            else:
                data = list(g[key])
            if key not in genotype:
                genotype[key] = []
            genotype[key].extend(data)
    
    where = {}
    for i, snp_id in enumerate(snps):
        where[snp_id] = i

    perm = [where[snp_id] for snp_id in sorted_snps]

    for key in genotype:
        if len(genotype[key]) != m[_config["MAPS_SIZE_ATTR"]]:
            raise Exception("Sample genotype and map size mismatch.")
        genotype[key] = [x for _, x in sorted(zip(perm, genotype[key]))]
    
    return genotype
        
                


def insert_file(file, individual_id=None):
    """Upload a file to the database, optionally linking it to an individual.

    Parameters
    ----------
    file                    File-like object to import
    individual_id=None      Identifier of the individual associated with the
                            the file. Currently has no special meaning and any
                            identifier (including non-existing ones) can be used.
    """
    if individual_id is None:
        _GFS.put(file, filename=os.path.basename(file.name))
    else:
        _GFS.put(file, filename=os.path.basename(file.name),
        **{_config["FILES_INDIVIDUAL_ATTR"]: individual_id})




def list_files(individual_id=None, name=None):
    """Search files in the database.
    
    Returns a list of dicts, each one containing metadata from
    a file.

    Parameters
    ----------
    individual_id=None      Match only files with the specified individual_id.
    name=None               Match only files with the specified original name.
    """
    query = {}
    if individual_id is not None:
        query.update({_config["FILES_INDIVIDUAL_ATTR"]: individual_id})
    if name is not None:
        query.update({"filename": name})

    return list(_db.fs.files.find(query, {"chunkSize": 0}))




def get_files(files):
    """Download a list of files from the database.

    The files are downloaded to the current working directory,
    maintaning their original names. Already existing files
    will be overwritten. If two files with the same name
    are specified, one will overwrite the other.

    Parameters
    ----------
    files   List of dicts describing the files to download, such as the one
            returned by list_files. Each dict must contain the "_id" field
            of the file to download.
    """
    for file_doc in files:
        grid_out = _GFS.get(file_doc["_id"])
        filename = grid_out.filename
        with open(filename, "wb+") as f:
            f.write(grid_out.read())




def import_map(map_reader, map_name,
               force_create_new=False,
               force_use_existing=False,
               report=False):
    """Import map into the database using a MapReader.

    Parameters
    ----------
    map_reader                  A MapReader instance.
    map_name                    The identifier the map will have inside the
                                database.
    force_create_new=False      If False, for each SNP being imported, tries
                                to find similar SNPs already on the database.
                                If at least one is found, the user is asked
                                to choose between creating a new SNP or using
                                one of the similars found.
                                SNPs are considered similar if their chromosomes
                                and positions match exactly. If True, new SNPs
                                are always created.
    force_use_existing=False    If True, will always try to use similar 
                                (see parameter above) SNPs instead of creating
                                new ones. If only one similar SNP is found, it
                                is used automatically. If more than one is
                                found, the user is asked to choose. If none is
                                found, a new SNP is created. This and the
                                previous parameter cannot be True at the same
                                time.
    report=False                If True, import results are printed after
                                finished. If False, nothing is printed.
    """

    if len(find_maps(id=map_name)) > 0:
        raise Exception("Map name already in use.")

    # Determine internal IDs for each SNP, possibly with user
    # interaction.
    nsnps = len(map_reader)
    snp_ids = __fill_snp_ids(map_reader, force_create_new, force_use_existing)

    if snp_ids is None:
        return

    # Insert new map into maps collection.
    map_doc = {"_id": map_name,
              _config["MAPS_SIZE_ATTR"]: nsnps,
              _config["MAPS_BLOCK_SIZE_ATTR"]: _config["SNPBLOCKS_SNPS_PER_BLOCK"]}
    map_doc.update(map_reader.map_meta())
    _MAPS.insert_one(map_doc)

    # Insert map snp list (both original order and sorted by id)
    # into map snps collection.
    BS = _config["MAPSNPS_MAX_LIST_SIZE"]
    snp_list = [x[0] for x in snp_ids]
    s_snp_list = sorted(snp_list)
    _MAPSNPS.insert_many((
        {_config["MAPSNPS_MAP_ATTR"]:map_name,
         _config["MAPSNPS_IDX_ATTR"]:i,
         _config["MAPSNPS_LIST_ATTR"]: snp_list[i*BS:i*BS+BS],
         _config["MAPSNPS_SORTED_LIST_ATTR"]: s_snp_list[i*BS:i*BS+BS]}
         for i in range(0, (nsnps-1)//BS + 1)))
    
    # Insert new SNPs into snps collection.
    new_snps = [__adjust_snp(snp, map_reader, snp_ids[i][0])
               for (i, snp) in enumerate(map_reader) if snp_ids[i][1]]
    if len(new_snps) > 0:
        _SNPS.insert_many(new_snps)

    # For each SNP (old or new), add the new map to the SNP's map list.
    _SNPS.bulk_write([UpdateOne({"_id": snp_ids[j][0]},
        {"$push": {_config["SNPS_MAPS_ATTR"]: map_name}}) for j in range(nsnps)])
    
    if report:
        print(f"Added map {map_name} with {nsnps} SNPs, " +
              f"{len(new_snps)} new SNPs created.")




def import_samples(sample_reader, map_name, id_map={}, report=False):
    """Import samples into the database using a SampleReader.

    Parameters
    ----------
    sample_reader   A SampleReader instance.
    map_name        The map id associated with the samples being imported.
                    Must exist on the database.
    id_map={}       A dict associating sample IDs (key) with individual
                    tatoos (value). Sample IDs which have an associated
                    individual tatoo will be associated to the individual
                    in question. If it doesn't exist, it'll be created.
                    If more than one exists with the same tatoo, the user
                    will be asked to choose between them.
    report=False    If True, import results are printed after finished.
                    IF False, nothing is printed.
    """
    try:
        m = find_maps(id=map_name)[0]
    except IndexError:
        raise Exception("Map not found.") from None

    snps, sorted_snps = get_map_snps(map_name)
    bsize = m[_config["MAPS_BLOCK_SIZE_ATTR"]]

    new_samples = 0
    new_blocks = 0
    new_individuals = 0
    old_individuals = 0

    for sample in sample_reader:
        genotype = sample.pop(sample_reader.SAMPLE_GENOTYPE)
        id = sample.pop(sample_reader.SAMPLE_ID)
       
        for key in genotype:
            if len(genotype[key]) != len(snps):
                raise Exception("Sample genotype and map size mismatch.")

        # Prepare sample object to be inserted.
        sample_key = {_config["SAMPLES_MAP_ATTR"]: map_name,
                      _config["SAMPLES_ID_ATTR"]: id}
        sample.update(sample_key)

        _SAMPLES.insert_one(sample)
        new_samples += 1

        # Sort genotype lists using snp id as key.
        for key in genotype:
            if isinstance(genotype[key], str):
                genotype[key] = "".join([x for _, x in sorted(zip(snps, genotype[key]))])
            else:
                genotype[key] = [str(x) for _, x in
                                 sorted(zip(snps, genotype[key]))]

        # Break genotype into blocks and insert into SNP blocks collection.
        current_block = 0
        for i in range(0, len(snps), bsize):
            b_genotype = {}
            for key in genotype:
                if isinstance(genotype[key], str):
                    b_genotype[key] = genotype[key][i:i+bsize]
                else:
                    b_genotype[key] = " " + " ".join(genotype[key][i:i+bsize])
            _SNPBLOCKS.insert_one({
                _config["SNPBLOCKS_MAP_ATTR"]: map_name,
                _config["SNPBLOCKS_SAMPLE_ATTR"]: id,
                _config["SNPBLOCKS_BLOCK_NUMBER"]: current_block,
                _config["SNPBLOCKS_GENOTYPE"]: b_genotype})
            new_blocks += 1
            current_block += 1

        # Try to associate the sample with an individual, possibly
        # interacting with the user.
        if id in id_map:
            individuals = find_individuals(tatoo=id_map[id])
            option = 0
            if len(individuals) > 1:
                option = __user_individual_choice(id_map[id], individuals)
            elif len(individuals) == 1:
                option = 1
            if option == 0:
                _INDS.insert_one({
                    "_id": __next_individual_id(),
                     _config["INDIVIDUALS_ID_LIST_ATTR"]: [id_map[id]],
                     _config["INDIVIDUALS_SAMPLE_LIST_ATTR"]: [sample_key]})
                new_individuals += 1
            else:
                _INDS.update_one(
                    {"_id": individuals[option-1]["_id"]},
                    {"$push": {_config["INDIVIDUALS_SAMPLE_LIST_ATTR"]: sample_key}})
                old_individuals += 1
        
    if report:
        print(f"{new_samples} samples added, {new_blocks} blocks, " +
              f"{new_individuals} individuals created, " +
              f"{old_individuals} pre-existing individuals updated.")




def export_map(id, map_writer, out_file_path):
    """Export map from database to file using a MapWriter.

    The SNPs of the map to be exported should contain
    all the fields required by the specific MapWriter.

    Parameters
    ----------
    id              Name of the map in the database to export.
    map_writer      A MapWriter instance.
    out_file_path   Path of file to export to.
    """
    snps, _ = get_map_snps(id)
    if len(snps) == 0:
        raise Exception("Map not found.")
    where = dict()
    for i in range(len(snps)):
        where[snps[i]] = i
    wsnps = [None] * len(snps)
    for snp in _SNPS.find({_config["SNPS_MAPS_ATTR"]:id}):
        wsnps[where[snp["_id"]]] = snp
    for snp in wsnps:
        __rev_adjust_snp(snp, map_writer)
    writer = map_writer(wsnps)
    writer.write(out_file_path)




def export_samples(samples, map, sample_writer, out_file_path):
    """Export samples to file using a SampleWriter.

    The samples to be exported should contain all the fields required by
    the specific SampleWriter.
    
    Parameters
    ----------
    samples         List containing within-map IDs of the samples to be exported.
                    If empty, all the samples associated with the map will be
                    exported.
    map             ID of map with which the samples to be exported are associated.
    sample_writer   A SampleWriter instance.
    out_file_path   Path of file to export to.
    """
    #TODO: optimize performance by reducing number of calls to find_sample.
    wsamples = []
    if len(samples) == 0:
        samples = [sample[_config["SAMPLES_ID_ATTR"]] for sample in find_sample(map=map)]
    for id in samples:
        current = {sample_writer.SAMPLE_ID: id,
                   sample_writer.SAMPLE_GENOTYPE: get_sample_data(id, map)}
        sample_info = find_sample(id, map)[0]
        sample_info.pop("_id")
        current.update(sample_info)
        wsamples.append(current)
    writer = sample_writer(wsamples)
    writer.write(out_file_path)



# Used for testing purposes.
def get_db_stats(scale=1):
    colls = [_config[key] for key in _config if key[-4:] == "COLL"]
    for coll in colls:
        _db.command("validate", coll, full=True)
    return _db.command("dbstats", scale=scale)



# Used for testing purposes.
#def create_individuals(individuals):
#    _INDS.insert_many(individuals)




def __reserve_snp_ids(cnt):
    doc = _COUNTERS.find_one_and_update({"_id": _config["SNPS_COLL"]}, 
        {"$inc": {_config["COUNTERS_SEQ_VALUE_ATTR"]: cnt}})
    return doc[_config["COUNTERS_SEQ_VALUE_ATTR"]]




def __adjust_snp(snp, map_reader, id=None):
    if id is not None:
        snp["_id"] = id
    if map_reader.SNP_NAME in snp:
        snp[_config["SNPS_NAME_ATTR"]] = snp.pop(map_reader.SNP_NAME)
    if map_reader.SNP_CHROM in snp:
        snp[_config["SNPS_CHROMOSOME_ATTR"]] = snp.pop(map_reader.SNP_CHROM)
    if map_reader.SNP_POS in snp:
        snp[_config["SNPS_POSITION_ATTR"]] = snp.pop(map_reader.SNP_POS)
    return snp




def __rev_adjust_snp(snp, map_writer):
    snp.pop("_id")
    if _config["SNPS_NAME_ATTR"] in snp:
        snp[map_writer.SNP_NAME] = snp.pop(_config["SNPS_NAME_ATTR"])
    if _config["SNPS_CHROMOSOME_ATTR"] in snp:
        snp[map_writer.SNP_CHROM] = snp.pop(_config["SNPS_CHROMOSOME_ATTR"])
    if _config["SNPS_POSITION_ATTR"] in snp:
        snp[map_writer.SNP_POS] = snp.pop(_config["SNPS_POSITION_ATTR"])




def __find_similar_snps(snp):
    query = {}
    for attr in [_config["SNPS_CHROMOSOME_ATTR"],
                 _config["SNPS_POSITION_ATTR"]]:
        if attr not in snp:
            return []
        query[attr] = snp[attr]
    return list(_SNPS.find(query))




def __user_snp_choice(snp, conflicts, force_use_existing):
        id = -1
        print(str(snp) + " is similar to the following database SNPs:")
        i = 1
        for conflict_snp in conflicts:
            print("(%d)" % i, end=' ')
            print(conflict_snp)
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




def __fill_snp_ids(map_reader, force_create_new, force_use_existing):
    if force_create_new and force_use_existing:
        raise Exception("force_create_new and force_use_existing cannot"
                        " be used simultaneously.") 
    next_id = __reserve_snp_ids(len(map_reader))
    snp_ids = []
    for snp in map_reader:
        __adjust_snp(snp, map_reader)
        if force_create_new:
            snp_ids.append((next_id, True))
            next_id += 1
            continue
        similar = __find_similar_snps(snp)
        if len(similar) == 0:
            snp_ids.append((next_id, True))
            next_id += 1
        elif len(similar) == 1 and force_use_existing:
            snp_ids.append((similar[0]["_id"], False))
        else:
            user_choice = __user_snp_choice(snp, similar,
                                            force_use_existing)
            if user_choice is None:
                return None
            if user_choice in {"e", "E"}:
                snp_ids.append((similar[0]["_id"], False))
                force_use_existing = True
            elif user_choice in {"0", "n", "N"}:
                snp_ids.append((next_id, True))
                next_id += 1
                force_create_new = user_choice in {"n", "N"}
            else:
                snp_ids.append((similar[int(user_choice)-1]["_id"], False))
    return snp_ids



    
def __user_individual_choice(tatoo, individuals):
    print("Ambigous match for individual %s:" % tatoo)
    i = 1
    for individual in individuals:
        print("(%d)" % i, end=' ')
        pprint(individual)
        i += 1
    question = "What do to? [0] - create new individual; "
    question += "[1] to [%d] - use i-th individual: " % (i-1)
    resp = None
    while resp == None:
        try:
            resp = int(input(question))
            if resp < 0 or resp > len(individuals):
                resp = None
        except KeyboardInterrupt:
            return None
        if resp is None:
            print("Invalid response.")
    return resp




def __next_individual_id():
    doc = _COUNTERS.find_one_and_update(
        {"_id": _config["INDIVIDUALS_COLL"]},
        {"$inc": {_config["COUNTERS_SEQ_VALUE_ATTR"]: 1}})
    return doc[_config["COUNTERS_SEQ_VALUE_ATTR"]]

