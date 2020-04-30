// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(U) 
  + fvm::div(phi, U)
  + MRF.DDt(U)
 ==
    fvOptions(U)
  + constEq.divTauThermo(U)  
  + elecM.Fe()/constEq.rho()      
);

fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

spSolverUFluid[i].solve
(
  UEqn
  ==
  ( thermo.hasGravity() 
    ?    
    fvc::reconstruct
    (
      (
       - thermo.ghf()*fvc::snGrad(thermo.rhok())
       - fvc::snGrad(p_rgh)
      ) * mesh.magSf()
    )
    :
    -fvc::grad(p_rgh)
  )
);

fvOptions.correct(U);
