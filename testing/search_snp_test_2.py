import snpdb
import random
import time
import sys

CHROM = snpdb.config["SNPS_CHROMOSOME_ATTR"]
POS = snpdb.config["SNPS_POSITION_ATTR"]
MAPS = snpdb.config["SNPS_MAPS_ATTR"]
SAMPLE_ID = snpdb.config["SAMPLES_ID_ATTR"]

snps = [snp for snp in snpdb.find_snp() if (CHROM in snp) and (POS in snp)]
query_snps = random.choices(snps, k=int(sys.argv[1]))

t = 0.
multiple_snps = 0
miss = 0
for qsnp in query_snps:
    
    map = random.choice(qsnp[MAPS])
    sample = random.choice(snpdb.find_sample(map=map))

    start = time.time()
    res = snpdb.find_snp(min_chrom=qsnp[CHROM],
                         max_chrom=qsnp[CHROM],
                         min_pos=qsnp[POS],
                         max_pos=qsnp[POS])
    res2 = snpdb.find_snp_of_sample(map, 
                                    sample[SAMPLE_ID],
                                    res[0]["_id"])
    t += (time.time() - start)
    if len(res) > 1:
        multiple_snps += 1
    if res2 is None:
        miss += 1


print(f"{len(query_snps)} in {t:.3f} s, " + 
      f"avg: {t/len(query_snps):.9f}, " + 
      f"multiple: {multiple_snps}, " + 
      f"miss: {miss}")
