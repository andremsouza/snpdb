
load("config.js");

/* Creates an object {key1: value1, key2: value2...}
to be passed to MongoDB functions. */
function keyValueObject(key1, value1, ...args) {
	if (args.length % 2 == 1)
		throw new Error("Even number of arguments required.");
	obj = {};
	obj[key1] = value1;
	for (let i = 0; i < args.length; i += 2)
		obj[args[i]] = args[i + 1];
	return obj;
}


db = connect(config.HOST + "/" + config.DB_NAME);
db.dropDatabase()

db.createCollection(config.SNPS_COLL);
db[config.SNPS_COLL].createIndex(keyValueObject(config.SNPS_NAME_ATTR, 1));
db[config.SNPS_COLL].createIndex(keyValueObject(config.SNPS_CHROMOSOME_ATTR, 1, config.SNPS_POSITION_ATTR, 1));
db[config.SNPS_COLL].createIndex(keyValueObject(config.SNPS_MAPS_ATTR, 1));

db.createCollection(config.INDIVIDUALS_COLL);
db[config.INDIVIDUALS_COLL].createIndex(keyValueObject(config.INDIVIDUALS_ID_LIST_ATTR, 1));
db[config.INDIVIDUALS_COLL].createIndex(keyValueObject(config.INDIVIDUALS_SAMPLE_LIST_ATTR + "." + config.SAMPLES_MAP_ATTR, 1));
db[config.INDIVIDUALS_COLL].createIndex(keyValueObject(config.INDIVIDUALS_SAMPLE_LIST_ATTR + "." + config.SAMPLES_ID_ATTR, 1));

db.createCollection(config.MAPS_COLL);

db.createCollection(config.MAPSNPS_COLL);
db[config.MAPSNPS_COLL].createIndex(keyValueObject(config.MAPSNPS_MAP_ATTR, 1, config.MAPSNPS_IDX_ATTR, 1), unique = true)

db.createCollection(config.SAMPLES_COLL)
db[config.SAMPLES_COLL].createIndex(keyValueObject(config.SAMPLES_MAP_ATTR, 1, config.SAMPLES_ID_ATTR, 1), unique = true)
db[config.SAMPLES_COLL].createIndex(keyValueObject(config.SAMPLES_ID_ATTR, 1))

db.createCollection(config.SNPBLOCKS_COLL);
db[config.SNPBLOCKS_COLL].createIndex(keyValueObject(config.SNPBLOCKS_MAP_ATTR, 1, config.SNPBLOCKS_SAMPLE_ATTR, 1, config.SNPBLOCKS_BLOCK_NUMBER, 1));

/* Collection for generating sequential ids. */
db.createCollection(config.COUNTERS_COLL);
db[config.COUNTERS_COLL].insert(keyValueObject("_id", config.SNPS_COLL, config.COUNTERS_SEQ_VALUE_ATTR, NumberInt(0)));
db[config.COUNTERS_COLL].insert(keyValueObject("_id", config.INDIVIDUALS_COLL, config.COUNTERS_SEQ_VALUE_ATTR, NumberInt(0)));

db.fs.data.createIndex(keyValueObject(config.FILES_INDIVIDUAL_ATTR, 1))
