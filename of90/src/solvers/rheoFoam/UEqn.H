// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(U) 
  + fvm::div(phi, U)
  + MRF.DDt(U)
 ==
    fvModels.source(U)
  + constEq.divTau(U)       
);

fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvConstraints.constrain(UEqn);
 
spSolverU->solve(UEqn == -fvc::grad(p));  

fvConstraints.constrain(U);
 
