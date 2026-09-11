S2-S10 Code and source data files

S2 contains all original members of the former S2 archive except two large, internally XZ-compressed RDS objects. Those two objects are distributed losslessly across S3-S6 and S7-S10, respectively, so every upload file remains below the PLOS 20 MB limit.

Extract S2, place S3-S10 in one directory, and run:
python reassemble_large_rds.py --archives-dir <directory-containing-S3-S10>

The script verifies every part and each reconstructed RDS object against SHA-256 values in large_rds_manifest.json. No scientific data were changed.
