Search the two diaPASEF runs of the 1k Pepyrus pool (runs 4960 and 4962 of MSV000096555) together, with the search
settings of the DDA module (non-specific, 8 to 12 residues, 700 to 2,000 Da, 15 ppm precursor and 20 ppm fragment
tolerance, the same four variable modifications, at most two).

A spectral-library search can only report library entries, so its FDP depends on how much entrapment the library
contains. **Upload the library as well** (TSV or parquet; binary `.speclib` files cannot be read): *r* is then
computed from the library instead of the FASTA. Without a library, the FASTA's *r* is used, which is only correct
for a library-free search of the benchmark FASTA. A library built from the sample itself contains almost no
entrapment; its FDP upper bound is then reported as not estimable.
