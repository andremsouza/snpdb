import import_map
import import_samples
import os
import time
import snpdb

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
    print(exts) 
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
        
        if _EXT_FORMAT[mapfileext] != _EXT_FORMAT[pedfileext]:
            raise Exception("Incompatible map and ped files.")

        fmt = _EXT_FORMAT[mapfileext]

        t_m = _stopwatch(import_map.import_map,
                         os.path.join(directory, name + mapfileext),
                         fmt,
                         name,
                         report=True,
                         **kwargs)
        
        if not maps_only:
            idfilename = None
            if idsfileext is not None:
                idfilename = name + idsfileext
            t_p = _stopwatch(import_samples.import_samples,
                            os.path.join(directory, name + pedfileext),
                            fmt,
                            name,
                            idfilename=os.path.join(directory, idfilename),
                            report=True)
        else:
            t_p = 0.0            

        stats = snpdb.get_db_stats(2 ** 20)
        data = stats["dataSize"]
        storage = stats["storageSize"]
        print(f"{name}: "+
              f"{t_m:.3f} s / {t_p:.3f} s / " +
              f"{data:.1f} MiB / {storage:.1f} MiB")
