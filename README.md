This is a Python application for managing single-nucleotide polymorphism data using MongoDB as backend. It seeks to be highly flexible and extensible. Currently it works with the `0125`, `PLINK PED/MAP`, `Illumina Final Report` and `VCF` formats, although exporting is only implemented for the first two. Also, it may not work for some variations of those formats.

Please refer to the documentation in `snpdb.py` for more info.


## Installation

Requires Python 3.6 or higher, MongoDB 4.0 or higher.

Run the `mongo_setup.js` script inside MongoDB (e.g. `load(mongo_setup.js)` inside `mongo` shell) to setup the database.
Then, make sure `HOST` and `DB_NAME` parameters are correct on `config.js`.
