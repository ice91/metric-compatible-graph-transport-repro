# Metric-compatible graph transport: reproducibility artifact

Reproducibility artifact for the exact finite-graph rank witnesses in
"Gauge Orbits, Cycle Completeness, and Reference Design in
Metric-Compatible Graph Transport."

Author: Chien-Chih Chen  
ORCID: [0000-0002-0019-1164](https://orcid.org/0000-0002-0019-1164)

## Scope

Version 1.0.0 reproduces one deterministic witness: the twelve-vertex,
seventeen-edge grid at sector `r = 2`, with the normalized reference
`g_v = I_2` and `L_e = I_2`. It computes the cycle-space ranks, the
internal Jacobian, the gauge tangent, and the restrictions of Designs A
and B.

The analytical theorems in the paper do not depend on this software.
This repository reproduces the deterministic finite-graph rank witnesses
reported in the article.

No external dataset is used. The graph is written in
`src/mcgt_repro/graph_fixture.py`.

## Python

Python 3.11 or Python 3.12. The standard library is sufficient. No
third-party numerical package is required.

## One-command reproduction

From the repository root:

```bash
./reproduce.sh
```

The command prints the witness and writes `reproduced_results.json` in
the current directory. Arithmetic is exact rational elimination.
No floating-point tolerance is used.

## Expected output

```text
Twelve-vertex witness
---------------------

|V|                     = 12
|E_free|                = 17
r                       = 2
q_max                   = 6
q                       = 4

rank(J_int)             = 84
dim ker(J_int)          = 100
rank(T_gauge)           = 92
dim H_miss              = 8

rank(A_A | T_gauge)     = 92
rank(A_B | T_gauge)     = 92
rank(A_A | ker J_int)   = 92
rank(A_B | ker J_int)   = 92

J_int(T_gauge)          = exact zero

arithmetic              = exact rational
floating-point tolerance = none

RESULT                  = PASS
```

## Where these numbers appear in the article

The twelve-vertex table in the article records `q_max = 6`, `q = 4`,
kernel dimension 100, gauge-tangent dimension 92, missing-cycle
dimension 8, and Design A/B rank 92 on the gauge tangent. The same
section, and the appendix that lists the edges and declared cycles,
records Design A/B rank 92 on the full internal kernel. In the real
coordinate space of that calculation, the internal Jacobian has rank 84
and kernel dimension 100.

## Tests

```bash
PYTHONPATH=src python3 -m unittest discover -s tests -v
```

## License

A license has not been selected. See `LICENSE`. This candidate is not
offered for public reuse until a license is added.
