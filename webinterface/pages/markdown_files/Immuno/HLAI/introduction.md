This module benchmarks **HLA class I peptide identification** on a sample of known content: a pool of
synthetic Pepyrus peptides (Manakongtreecheep, Ctortecka *et al.*, *Nat. Biotechnol.* 2026; MassIVE
MSV000096555). Because the content of the sample is known, two numbers are measured rather than assumed:

- **Pool peptides identified**: the number of designed pool peptides that a workflow identifies per run.
- **False discovery proportion (FDP)**: the fraction of accepted peptides that cannot be in the sample,
  estimated with the combined entrapment method (Wen *et al.*, *Nat. Methods* 2025).

Every workflow is held to the same rule: rank-1 target identifications with a peptide-level q-value of at most 1%.
