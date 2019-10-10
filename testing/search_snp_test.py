import snpdb
import random
import time
import sys

CHROM = snpdb.config["SNPS_CHROMOSOME_ATTR"]
POS = snpdb.config["SNPS_POSITION_ATTR"]

snps = sorted([(snp[CHROM], snp[POS]) for snp in snpdb.find_snp()
               if (CHROM in snp) and (POS in snp)])

queries = random.choices(snps, k=int(sys.argv[1]))

t = 0.
for query in queries:
    start = time.time()
    res = snpdb.find_snp(min_chrom=query[0],
                         max_chrom=query[0],
                         min_pos=query[1],
                         max_pos=query[1])
    t += (time.time() - start)
    assert(len(res) > 0)


print(f"{len(queries)} in {t:.3f} s, " + 
      f"avg: {t/len(queries):.9f}")
