# Run405 exact primitive comparator

Board116 reserved/linked before computation. Deterministic seed0/no randomness.
Explicit19CX word: Bellprefix2+R2+CSdg2+CCX6+CH1+CCH6. One numerical checker
invocation30sec and one exact cyclotomic checker60sec; external100sec/one thread.
SourceSHA a6ae1b4cbd754e6f35f308ed667034011e9fb5eefc9ab89999d5edc433247bb9,
configSHA f42ac287926163c2f61058f9aa46af6ea3bccdf06eace6a812126d8fc29c0d93.
Run once only. Ledgerstatus/code guard prevents post-close repetition. Preserve
candidate/logs/results even on failure. This is overbudget comparator, not12CX
success or incumbent improvement. Actual CH chronology implements
Ry(+pi/4) Z Ry(-pi/4)=H; the frozen source comment reverses those signs, but
its gate list and independent by-hand review have the correct chronology.

```bash
timeout --signal=TERM --kill-after=5s 100s env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 python3 experiments/runs/405/check_baseline.py > experiments/runs/405/stdout.txt 2> experiments/runs/405/stderr.txt
```
