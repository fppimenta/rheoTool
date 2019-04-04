//- Check type of solver to be used
const dictionary* cSDict = mesh.solutionDict().subDictPtr("coupledSolvers");
 
bool solveCoupled(false);
if (cSDict != NULL)
{
 const dictionary* cSDictUpTau = cSDict->subDictPtr("Uptau");
 if (cSDictUpTau != NULL)
   solveCoupled = readBool(cSDictUpTau->lookup("solveCoupledUp"));
} 
 
//- Create coupled solver if needed
autoPtr<coupledSolver>  cps
(
  solveCoupled
  ? 
  new coupledSolver("Uptau", mesh)
  :
  NULL
);

//- Add fields to coupled solver
if (solveCoupled)
{
 cps->insertField(p);
 cps->insertField(U);  
}

//- Create segregated solvers for p|U if needed
autoPtr<sparseSolver<scalar> > spSolverp
(
  !solveCoupled
  ? 
  sparseSolver<scalar>::New(p, mesh, mesh.solutionDict()).ptr()
  :
  NULL
);

autoPtr<sparseSolver<vector> > spSolverU
(
  !solveCoupled
  ? 
  sparseSolver<vector>::New(U, mesh, mesh.solutionDict()).ptr()
  :
  NULL
);


