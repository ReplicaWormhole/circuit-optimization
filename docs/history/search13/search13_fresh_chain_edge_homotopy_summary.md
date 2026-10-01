# Smooth relocation of one interaction in the fresh chain lead

Run344 follows two specified paths: move slot6 from its source interaction to03, or slot8 to02. The powered interaction is C(a)=(I+C)/2+exp(i*pi*a)(I-C)/2. Since C is a Hermitian involution, this is unitary; C(0)=I and C(1)=C. At each varied slot the matrix is C_new(1-a) C_old(a). Thus a=1 gives the source CNOT and a=0 the proposed replacement CNOT. Interior points use non-native powered entanglers and are outside the counted gate model.

Each path visits a=1,.875,.75,.5,.25,.125,0 with at most150 local iterations per stage, then at most300 native-endpoint polish iterations. All168 local parameters vary. Both final native candidates and both polished candidates have13 ordinary CNOTs and independently fail the cycle checker. Polished losses are0.1750182 and0.4970417, with maximum entries0.3582416 and0.5317426. Both are worse than the source. No interior point is offered as a13CNOT circuit or upper bound.

Checks: source loss roundtrip, endpoint identities, unitarity of every powered-product stage, native endpoint serialization and exact CNOT count, independent NumPy endpoint losses within1e-11, original cycle validity reevaluation. Adversarial review found no confirmed issue; non-native stages and local-search limits are explicit. This is bounded search evidence, not a topology exclusion.

```bash
python3 search13_fresh_chain_edge_homotopy.py --maxiter 150 --polish-maxiter 300
python3 experiment_log.py show 344
```

Reproduce in a separate copy without the existing result file. Script, stage records, four native candidate gate lists and this summary are kept artifacts.
