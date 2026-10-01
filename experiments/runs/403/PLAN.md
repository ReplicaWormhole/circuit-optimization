# Run403 independent preflight

Reserved/linked board103 before execution. Three fixed parameter points;
6 full-family builds (Torch/independent NumPy), 6 native/compiled list builds,
2 prefix/direct decoder products. No objectives, AD, optimizers, or searches.
One invocation, one thread, external30sec. Frozen checker/config hashes
9422a752053361ede634b2dde88933def820c7525957724412684e33975cca92 /
4cb144f1b4046a5e2138460127baf2c4a1d1bc011ae21174988ca82de6dc9432.
Bind correctedfit402; gate chronology and prefix orientation independently checked.

```bash
timeout --signal=TERM --kill-after=5s 30s env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 python3 experiments/runs/403/preflight.py > experiments/runs/403/stdout.txt 2> experiments/runs/403/stderr.txt
```
