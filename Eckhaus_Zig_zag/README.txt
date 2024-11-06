README file

MATLAB 2023 scripts for chapter 
“Numerical continuation in PDE -- wavetrains and invasion in the Swift-Hohenberg equation” 
by 
Ryan Goh (Boston University, USA), rgoh@bu.edu
David Lloyd (University of Surrey), d.lloyd@surrey.ac.uk
Jens D.M. Rademacher (Universit\"at Hamburg), jens.rademacher@uni-hamburg.de


To compute and plot zig-zag and Eckhaus boundaries and more, 
and thus generate the figures do the following:

Run each of these cell-by-cell!

1. Computation of upper Eckhaus boundary: continuation_EK_boundary_upper.m

2. Computation of zigzag boundary: continuation_EK_boundary_lower.m

3. Computation of lower Eckhaus boundary (slower): continuation_EK_boundary_lower.m

4. Plot results and reproduce figures: plot_EK_ZZ_boundaries.m

Continuation uses code developped by Danielle Avitabile (SecantContinuation.m) 
with small modification for stopping criteria.
