    MRF.correctBoundaryVelocity(U);

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
     ==
       phaseChange.SU(rho, rhoPhi, U)
     + fvModels.source(rho, U)
     + mixture.divTauMF(U)   
    );

    UEqn.relax();

    fvConstraints.constrain(UEqn);

    if (pimple.momentumPredictor())
    {        
      spSolverU->solve
      (
       UEqn
       ==
       fvc::reconstruct
       (
          (
            mixture.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
          - fvc::snGrad(p_rgh)
         ) * mesh.magSf()
       )
      );

      fvConstraints.constrain(U);
    }
