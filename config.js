const config = {
	"HOST": "localhost",
	"DB_NAME": "snpdb",

	"SNPS_COLL": "snps",
	"SNPS_NAME_ATTR": "i",
	"SNPS_CHROMOSOME_ATTR": "c",
	"SNPS_POSITION_ATTR": "p",
	"SNPS_MAPS_ATTR": "m",

	"INDIVIDUALS_COLL": "individuals",
	"INDIVIDUALS_ID_LIST_ATTR": "tatoos",
	"INDIVIDUALS_SAMPLE_LIST_ATTR": "samples",

	"MAPS_COLL": "maps",
	"MAPS_BLOCK_SIZE_ATTR": "block",
	"MAPS_FORMAT_ATTR": "format",
	"MAPS_SIZE_ATTR": "size",

	"MAPSNPS_COLL": "mapsnps",
	"MAPSNPS_MAP_ATTR": "map",
	"MAPSNPS_IDX_ATTR": "i",
	"MAPSNPS_LIST_ATTR": "snps",
	"MAPSNPS_SORTED_LIST_ATTR": "ssnps",
	"MAPSNPS_MAX_LIST_SIZE": 100000,

	"SAMPLES_COLL": "samples",
	"SAMPLES_MAP_ATTR": "map",
	"SAMPLES_ID_ATTR": "id",

	"SNPBLOCKS_COLL": "snpblocks",
	"SNPBLOCKS_MAP_ATTR": "m",
	"SNPBLOCKS_BLOCK_NUMBER": "no",
	"SNPBLOCKS_SAMPLE_ATTR": "s",
	"SNPBLOCKS_GENOTYPE": "g",
	"SNPBLOCKS_SNPS_PER_BLOCK": 10000,

	"COUNTERS_COLL": "counters",
	"COUNTERS_SEQ_VALUE_ATTR": "next",

	"FILES_INDIVIDUAL_ATTR": "individual"
}

