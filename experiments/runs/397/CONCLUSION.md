# Run 397: unexecuted reservation

Static review found an uncaught timeout path in the assignment-gap loop that
could discard partial output. No assignment-gap search, gradient, Hessian or
fit was executed. The ledger is terminal failed with no candidate.
The frozen script/config, plan and reviewer preflight are preserved.

Run 398 is a separately reserved replacement. It catches assignment-loop and
derivative timeouts together and checks the deadline between manual AD Hessian
rows, retaining partial rows. No curvature conclusion follows from run 397.
