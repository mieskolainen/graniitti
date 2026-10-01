# HEPData downloads

Run `python HEPData/download.py` from the repository root to download every JSON table for every `HEPData-ins<ID>-v<VERSION>-json` directory under `HEPData/`. The directory names select the records and fix their versions. New records can be included by creating a directory with this naming convention in the appropriate category.

Use `python HEPData/download.py PHOTOPROD` for one category, `python HEPData/download.py --list` to list the selected records, or `python HEPData/download.py --output tmp/hepdata` to download into a separate directory. Requires Python and `requests` from the `graniitti` environment.

The script discovers all tables from the [HEPData JSON service](https://www.hepdata.net/formats), including correlation tables, and verifies each downloaded table DOI before writing the original response bytes. It preserves existing filenames by DOI and uses HEPData filenames for new tables. Changed files are retained as `filename._old`, with further `._old` suffixes when needed. Identical files are left untouched. Failed downloads produce a nonzero exit status. Non-JSON ancillary files and Rivet files are outside this download.
