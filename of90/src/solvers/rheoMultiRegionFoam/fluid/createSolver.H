//- Check type of solver to be used
const dictionary* cSDict = fluidRegions[i].solutionDict().subDictPtr("coupledSolvers");
 
solveCoupledFluid[i]=false;
if (cSDict != NULL)
{
 const dictionary* cSDictUpTau = cSDict->subDictPtr("Uptau");
 if (cSDictUpTau != NULL)
   solveCoupledFluid[i] = readBool(cSDictUpTau->lookup("solveCoupledUp"));
} 

// We cannot solve p_rgh-U coupled if there is gravity. Check for it.
if (solveCoupledFluid[i] && thermoFluid[i].hasGravity())
{
  FatalErrorIn("rheoMultiRegionFoam\n")
  << "\nIt is not possible to solve p-U coupled when the thermo model for region: "
  << fluidRegions[i].name() << " enables gravity. Either set Uptau as false or disable gravity.\n"
  << abort(FatalError);
}
 
//- Create coupled solver if needed
cpsFluid.set
( 
  i,
  solveCoupledFluid[i]
  ? 
  new coupledSolver("Uptau",fluidRegions[i].name(),runTime,fluidRegions[i].solutionDict())
  :
  NULL
);

//- Add fields to coupled solver
if (solveCoupledFluid[i])
{
 cpsFluid[i].insertMesh(fluidRegions[i]);
 cpsFluid[i].insertField(p_rghFluid[i]);
 cpsFluid[i].insertField(UFluid[i]);  
}

//- Create segregated solvers for p|U if needed
spSolverp_rghFluid.set
(  
  i,
  !solveCoupledFluid[i]
  ? 
  sparseSolver<scalar>::New(p_rghFluid[i], fluidRegions[i], fluidRegions[i].solutionDict()).ptr()
  :
  NULL
);

spSolverUFluid.set
(  
  i,
  !solveCoupledFluid[i]
  ? 
  sparseSolver<vector>::New(UFluid[i], fluidRegions[i], fluidRegions[i].solutionDict()).ptr()
  :
  NULL
);
 

