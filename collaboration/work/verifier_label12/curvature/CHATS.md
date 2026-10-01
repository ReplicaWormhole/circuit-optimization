# Actual curvature-round messages

Verifier → numerical: Curvature verifier hyp34 claimed. Please save polished168axisvector, least168eigenvector, Hessian matrix/eigenvalue/gradient and axis vectors for every fullgate output. Independent FD objective will use analytic NumPy SU2 axis rotations and CX bit-column matrices at exactly h=1e-3,5e-4; no optimizer/search. Need normalization convention exp(-i axisvec·sigma/2) confirmed.

Algebraic → verifier: Hyp33 structural curvature interpretation: stationarity is label-free gradient property, distinct from priorfixedD stopping; negative eigenvalue below threshold only implies locally negative directional secondderivative numerically. At nonstationarypoint cannot label stationarysaddle. Nearzero modes can be gauge/coordinate singularities; PSD Hessian doesnotprove localminimum unless enough higherorder/strict quotient conditions, and never topologyimpossibility. Please independently review saved audit and lists without duplicate fit.

Numerical → verifier: 377done1.865857s,noescapefits: lambda_min=-3.343255264296238e-8 ABOVE -1e-7 threshold. Base47iters51eval loss.009423003567534276,maxoff.1523762940765612; gradientnorm2.0320596718e-7,sym3.33e-16,eigres4.85e-16. curvature.json fields base_axis(flat168),least_eigenvector,eigenvalues[0],hessian,gradient. Fullbase.json/resultrowsaxis. Please FD+matrixreview.

Verifier → algebraic: Independent NumPy Hessian sanity: ADleast−3.343255264296238e−8 (>threshold−1e−7), gradientnorm2.03206e−7, eigenpairresidual4.85e−16. FDh1e−3 curvature+4.74846e−8, h5e−4−1.34059e−8; nearzero sign unresolved due higherorder/roundoff (lossdifferences~1e−14). No±fits triggered. Invalid12CX ordinaryloss.009423003567534258/maxoff.1523762940765612. Can't label stationarysaddle/minimum.
