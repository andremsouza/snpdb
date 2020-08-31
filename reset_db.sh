#!/bin/bash
if [ -n "$1" ]; then # start mongd as a deamon
    echo "Shutting down MongoDB"
    /home/rodrigo/anaconda3/envs/snpdb/bin/mongo -eval "db.adminCommand({shutdown: 1})" ||
        echo "Warning: MongoDB server not being executed"
    echo "Deleting data directory"
    rm -rf /home/rodrigo/snpdb/mongodb_data/ --verbose
    mkdir /home/rodrigo/snpdb/mongodb_data --verbose
    echo "Starting database and loading schema"
    /home/rodrigo/anaconda3/pkgs/numactl-2.0.11-0/bin/numactl --interleave=all /home/rodrigo/anaconda3/envs/snpdb/bin/mongod --dbpath /home/rodrigo/snpdb/mongodb_data/ --wiredTigerCollectionBlockCompressor $1 --fork --logpath /home/rodrigo/snpdb/mongod.log --logappend &&
        /home/rodrigo/anaconda3/envs/snpdb/bin/mongo --eval "load('/home/rodrigo/snpdb/mongo_setup.js')" &&
        exit 0
else
    echo "usage: start_mongo.sh {snappy|zlib}" &&
        exit 1
fi
