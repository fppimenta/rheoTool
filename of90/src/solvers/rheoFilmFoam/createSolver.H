autoPtr<sparseSolver<vector> > spSolverU
(
  sparseSolver<vector>::New(U, mesh, mesh.solutionDict()).ptr()
);


