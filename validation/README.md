# Reproducible public-data validation

## Public 10x 5-prime Oxford Nanopore dataset

The validation uses the first 10,000 complete FASTQ records from the public
10x Genomics lung carcinoma DTC Oxford Nanopore dataset. The following command
downloads only the first 64 MiB of the source archive, extracts the fixed subset,
runs AnchorScope and checks the report against the values in
`public_10x_lung_ont/expected_results.json`:

```bash
python -B validation/public_10x_lung_ont/reproduce_public_example.py
```

The validation checks the FASTQ SHA-256, selected summary fields, count invariants
and the presence of all report outputs. To recheck an existing local run without
downloading or rerunning AnchorScope:

```bash
python -B validation/public_10x_lung_ont/validate_public_results.py
```

The public subset demonstrates report behaviour and does not contain curated
read-level structural truth. It must not be used to estimate sensitivity or
specificity.
