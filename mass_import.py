import import_map
import import_samples
import os
import time
import snpdb
import random


def _stopwatch(f, *args, **kwargs):
    start = time.time()
    f(*args, **kwargs)
    return time.time() - start

_EXT_FORMAT = {".0125map": import_map.Z125,
               ".0125ped": import_map.Z125,
               ".plmap": import_map.PLINK,
               ".plped": import_map.PLINK,
               ".vcf": import_map.VCF,
               ".fr": import_map.ILMFR}

_MAP_EXTS = {".0125map", ".plmap"}
_PED_EXTS = {".0125ped", ".plped", ".vcf", ".fr"}
_IDS_EXTS = {".ids"}





def mass_import(directory, maps_only=False, **kwargs):
    exts = {}
    for filename in os.listdir(directory):
        name, ext = os.path.splitext(filename)
        if name not in exts:
            exts[name] = {ext}
        else:
            exts[name].add(ext)
    for name in exts:
        imap = exts[name].intersection(_MAP_EXTS)
        iped = exts[name].intersection(_PED_EXTS)
        iids = exts[name].intersection(_IDS_EXTS)

        if len(imap) > 1:
            raise Exception("Multiple map files with the same name.")
        if len(iped) > 1:
            raise Exception("Multiple ped files with the same name.")
        if len(iids) > 1:
            raise Exception("Multiple id files with the same name.")
        
        mapfileext, pedfileext, idsfileext = None, None, None
        if len(imap) == 1:
            mapfileext = imap.pop()
        if len(iped) == 1:
            pedfileext = iped.pop()
        if len(iids) == 1:
            idsfileext = iids.pop()
        
        if mapfileext is None:
            mapfileext = pedfileext
        
        if (mapfileext is not None and pedfileext is not None and
            _EXT_FORMAT[mapfileext] != _EXT_FORMAT[pedfileext]):
            raise Exception("Incompatible map and ped files.")

        fmt = _EXT_FORMAT[mapfileext]

        t_m = _stopwatch(import_map.import_map,
                         os.path.join(directory, name + mapfileext),
                         fmt,
                         name,
                         report=True,
                         **kwargs)
        
        if not maps_only and pedfileext is not None:
            idfilename = None
            if idsfileext is not None:
                idfilename = os.path.join(directory, name + idsfileext)
            t_p = _stopwatch(import_samples.import_samples,
                            os.path.join(directory, name + pedfileext),
                            fmt,
                            name,
                            idfilename=idfilename,
                            report=True)
        else:
            t_p = 0.0            

        stats = snpdb.get_db_stats(2 ** 20)
        data = stats["dataSize"]
        storage = stats["storageSize"]
        print(f"{name}: "+
              f"map: {t_m:.3f} s, ped: {t_p:.3f} s, " +
              f"db size (raw): {data:.1f} MiB, compressed: {storage:.1f} MiB")
    



def create_individuals(n):
    inds = [{"_id": f"A{str(i+1)}"} for i in range(n)]
    snpdb.create_individuals(inds)




def import_media_randomly(directory):
    ids = [ind["_id"] for ind in snpdb.find_individuals()]
    total_t =  0.
    for filename in os.listdir(directory):
        id = random.choice(ids)
        f = open(os.path.join(directory, filename), "rb")
        t = _stopwatch(snpdb.insert_file, f, id)
        total_t += t
        f.close()
        print(f"Added file {filename} to individual {id} in {t:.3f} s.")
    
    time.sleep(61)
    stats = snpdb.get_db_stats(2 ** 20)
    data = stats["dataSize"]
    storage = stats["storageSize"]

    print(f"Total time: {total_t:.3f} s, " +
           f"storage: {data:.1f} MiB(raw), {storage:.1f} MiB (compressed)")
