import import_map
import import_samples
import os
import time
import snpdb

def _stopwatch(f, *args, **kwargs):
    start = time.time()
    f(*args, **kwargs)
    return time.time() - start

def mass_import(directory, fmt, **kwargs):
    for filename in sorted(os.listdir(directory)):
        mfilename = filename
        pfilename = filename
        if fmt in (import_map.PLINK, import_map.Z125):
            if filename[0:3] == "map":
                pfilename = "ped" + filename[3:]
            else:
                continue
        t_m = _stopwatch(import_map.import_map,
                        os.path.join(directory, mfilename),
                        fmt,
                        mfilename,
                        report=True,
                        **kwargs)
        t_p = _stopwatch(import_samples.import_samples,
                        os.path.join(directory, pfilename),
                        fmt,
                        mfilename,
                        report=True)
        stats = snpdb.get_db_stats(2 ** 20)
        data = stats["dataSize"]
        storage = stats["storageSize"]
        print(f"{os.path.join(directory, filename)}: "+
              f"{t_m:.3f} s / {t_p:.3f} s / " +
              f"{data:.1f} MiB / {storage:.1f} MiB")
