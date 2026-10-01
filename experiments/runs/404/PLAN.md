# Run404 refrozen independent preflight

Board107 reserved/linked before execution. Corrects Torch module version from
zeroevaluation403, preserving fit402 and scientific inputs. Frozen config SHA
24f6cbf7d81b98e2ce6c4a9b5e64fbc2d92d9e9f2c649a3f46c432c62e41cd24;
checker SHA9422a752053361ede634b2dde88933def820c7525957724412684e33975cca92.
Budget3fixedpoints,6Torch/NumPyfamily+6gatelist+2prefix/directproducts,
14top-level builds,one invocation,one thread,30sec. No targetloss/AD/optimizer.

```bash
timeout --signal=TERM --kill-after=5s 30s env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 python3 experiments/runs/404/preflight.py > experiments/runs/404/stdout.txt 2> experiments/runs/404/stderr.txt
```
