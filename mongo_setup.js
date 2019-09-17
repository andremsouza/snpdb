
const HOST = "localhost/snpdb";

/* Names for collections and their indexed attributes. 
Single-letter attribute names save space on dense collections. */

const SNPS_COLL = "snps";
const SNPS_ID_LIST_ATTR = "i";
const SNPS_CHROMOSOME_ATTR = "c";
const SNPS_POSITION_ATTR = "p";
const SNPS_MAPS_ATTR = "m";

const INDIVIDUALS_COLL = "individuals";
const INDIVIDUALS_ID_LIST_ATTR = "tatoos";
const INDIVIDUALS_SAMPLE_LIST_ATTR = "samples";

const SNPBLOCKS_COLL = "snpblocks";
const SNPBLOCKS_SAMPLE_ATTR = "s";
const SNPBLOCKS_SNP_LIST_ATTR = "a";
const SNPBLOCKS_SNP_ID_INSIDE_LIST = "i";

const MAPS_COLL = "maps";
const MAPS_SNP_LIST_ATTR = "snps";

const COUNTERS_COLL = "counters";
const COUNTERS_SEQ_VALUE_ATTR = "next";

/* Creates an object {key1: value1, key2: value2...}
to be passed to MongoDB functions. */
function keyValueObject(key1, value1, ...args)
{
	if (args.length%2 == 1)
		throw new Error("Even number of arguments required.");
	obj = {};
	obj[key1] = value1;
	for (let i = 0; i < args.length; i += 2)
		obj[args[i]] = args[i+1];
	return obj;
}


let db = connect(HOST);

db.createCollection(SNPS_COLL);
db[SNPS_COLL].createIndex(keyValueObject(SNPS_ID_LIST_ATTR, 1));
db[SNPS_COLL].createIndex(keyValueObject(SNPS_CHROMOSOME_ATTR, 1, SNPS_POSITION_ATTR, 1));
db[SNPS_COLL].createIndex(keyValueObject(SNPS_MAPS_ATTR, 1));

db.createCollection(INDIVIDUALS_COLL);
db[INDIVIDUALS_COLL].createIndex(keyValueObject(INDIVIDUALS_ID_LIST_ATTR, 1));
db[INDIVIDUALS_COLL].createIndex(keyValueObject(INDIVIDUALS_SAMPLE_LIST_ATTR, 1), unique=true);

db.createCollection(SNPBLOCKS_COLL);
db[SNPBLOCKS_COLL].createIndex(keyValueObject(SNPBLOCKS_SAMPLE_ATTR, 1, SNPBLOCKS_SNP_LIST_ATTR + ".0." + SNPBLOCKS_SNP_ID_INSIDE_LIST, 1));

db.createCollection(MAPS_COLL);
db[MAPS_COLL].createIndex(keyValueObject(MAPS_SNP_LIST_ATTR, 1));

/* Collection for generating sequential ids. */
db.createCollection(COUNTERS_COLL);
db[COUNTERS_COLL].insert(keyValueObject("_id", SNPS_COLL, COUNTERS_SEQ_VALUE_ATTR, NumberInt(0)));
