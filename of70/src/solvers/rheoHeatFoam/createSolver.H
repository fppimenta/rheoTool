//- Check type of solver to be used
const dictionary* cSDict = mesh.solutionDict().subDictPtr("coupledSolvers");
 
bool solveCoupled(false);
if (cSDict != NULL)
{
 const dictionary* cSDictUpTau = cSDict->subDictPtr("Uptau");
 if (cSDictUpTau != NULL)
   solveCoupled = readBool(cSDictUpTau->lookup("solveCoupledUp"));
} 

// We cannot solve p_rgh-U coupled if there is gravity. Check for it.
if (solveCoupled && thermo().hasGravity())
{
  FatalErrorIn("rheoHeatFoam\n")
  << "\nIt is not possible to solve p-U coupled when the thermo model enables gravity. Either set Uptau as false or disable gravity.\n"
  << abort(FatalError);
}

// Dummy null autoptr
autoPtr<coupledSolver> cpsT(NULL);
 
//- Create coupled solver if needed
autoPtr<coupledSolver>  cps
(
  solveCoupled
  ? 
  new coupledSolver("Uptau",mesh.name(),runTime,mesh.solutionDict())
  :
  NULL
);

//- Add fields to coupled solver
if (solveCoupled)
{
 cps->insertMesh(mesh);
 cps->insertField(p_rgh);
 cps->insertField(U);  
}

//- Create segregated solvers for p|U if needed
autoPtr<sparseSolver<scalar> > spSolverp_rgh
(
  !solveCoupled
  ? 
  sparseSolver<scalar>::New(p_rgh, mesh, mesh.solutionDict()).ptr()
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


