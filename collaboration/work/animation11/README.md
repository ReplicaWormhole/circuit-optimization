# Circuit matrix video

`video/circuit_matrices.mp4` shows the accepted 11-CNOT, 30-gate circuit.
Every state shows both complete 16 × 16 matrices: the cumulative circuit
`U_k = G_k ... G_1`, and the conjugated shift `U_k V4 U_k†`.
The initial identity and each post-gate state appear for two seconds, with
three additional seconds on the final state: 65 seconds in 4K at 10 fps.
Repeated video frames hold each matrix still; no matrix interpolation is used.

The circuit is chronological from left to right, with the gate just applied
highlighted in gold. Basis labels run from 0000 to 1111, q0 most significant.
Entry labels are rounded to three decimals; only magnitudes below 1e-10 are
classified as zero. All global phases are retained. Color records phase
(hue) and magnitude (brightness). Numerical visualization is not a new proof.

Reproduce from the repository root, choosing a fresh output folder:

```bash
python3 collaboration/work/animation11/render_video.py --output /tmp/animation11-new
ffmpeg -hide_banner -loglevel error -n -framerate 1/2 -i /tmp/animation11-new/frames/step_%03d.png -vf fps=10,tpad=stop_mode=clone:stop_duration=3 -c:v libx264 -preset veryfast -crf 18 -threads 2 -pix_fmt yuv420p -movflags +faststart /tmp/animation11-new/circuit_matrices.mp4
python3 collaboration/work/animation11/verify_matrices.py /tmp/animation11-new
```

Uses existing NumPy, Pillow, FFmpeg and system DejaVu fonts; no installation.
The renderer preserves existing output folders rather than overwriting them.
`video/matrices.npz` retains all 31 unrounded matrix pairs, and
`video/frames/` retains each full-resolution state. `video/manifest.json`
records the input and renderer hashes and the shared numerical check.
This is a presentation of an accepted circuit, not a new search attempt.

Verification: all 30 prefixes independently recomposed by tensor and bit
operations agree exactly in floating-point arithmetic; all 31 unitarity
checks pass at 1e-12; the final diagonal agrees with the exact acceptance
record within 2.11e-15. FFprobe confirms H.264, 3840 × 2160, 10 fps.
Visual inspection covered the initial state, a complex intermediate state,
and the final diagonal. Adversarial review found no outstanding defect.
Residual limitation: numerical entry labels are rounded and the MP4 is lossy;
the PNGs and numerical data are retained for closer inspection.

## Faster version

`video/circuit_matrices_3.5x.mp4` runs 3.5 times faster (about 18.57 seconds).
It keeps all 650 encoded frames at 35 fps, so all 31 matrix states remain
visible. Each ordinary state lasts about 0.571 seconds.

```bash
ffmpeg -hide_banner -loglevel error -n -i collaboration/work/animation11/video/circuit_matrices.mp4 -vf settb=AVTB,setpts=PTS/3.5 -r 35 -frames:v 650 -c:v libx264 -preset veryfast -crf 18 -threads 2 -pix_fmt yuv420p -movflags +faststart collaboration/work/animation11/video/circuit_matrices_3.5x.mp4
```
