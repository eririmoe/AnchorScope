# Benchmark design and interpretation

## Deterministic synthetic structure benchmark

`benchmarks/benchmark_synthetic.py` creates eight equally sized scenarios using two adapter-scale anchors and an 80-base payload:

1. exact full structure;
2. one anchor substitution;
3. one anchor insertion;
4. one anchor deletion;
5. missing 5-prime anchor;
6. reversed anchor order;
7. two complete concatemer cycles;
8. reverse-complement full structure.

Substitution, insertion, deletion, and reverse scenarios all have the expected label `full_structure`. Seeded payloads are rejected if they contain any sequence within one Levenshtein edit of an anchor or its reverse complement. This prevents accidental, unlabelled approximate anchors from invalidating the truth set.

The benchmark evaluates identical reads with `max_edits=1` and legacy `max_mismatches=1`, recording the confusion matrix, overall structure accuracy, throughput, Python/platform details, matcher backend, Git commit, and dirty-tree state.

Run:

```bash
python -B benchmarks/benchmark_synthetic.py \
  --reads-per-class 100 \
  --seed 20260710 \
  --output benchmarks/results/synthetic_benchmark.json
```

## Interpretation of the bundled result

On the bundled 800-read truth set, the indel-aware mode classified all structures correctly, whereas the Hamming mode missed the insertion and deletion cases. Throughput is environment- and backend-specific and must not be generalized from the bundled Windows/Python fallback run.

This benchmark establishes deterministic algorithmic correctness for the represented error model. It does not establish biological sensitivity or specificity. Broader performance claims require independent aligner comparisons, curated real-read truth, multiple protocols and chemistries, and resource measurements on a documented compute environment.
