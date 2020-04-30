//------- U -----------//
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

volScalarField rAU(1.0/UEqn.A());

// UEqn
cps->insertEquation
(
 U.name(), 
 U.name(),
 UEqn
);

// -grad(p)
cps->insertEquation
(
 U.name(),
 p_rgh.name(),
 fvmb::grad(1., p_rgh)
);

// Insert div(tau) if tau is to be solved coupled,
// otherwise do no nothing (div(tau) becomes inserted explicitly
// in UEqn) 
constEq.divTauImplCoupled();

// Insert -rhoE*grad(psi_) if ci-psi are to be solved coupled,
// otherwise do no nothing (-rhoE*grad(psi_) becomes inserted
// explicitly in UEqn). 
elecM.FeImplCoupled(constEq.rho());

//------- p -----------//

// Pressure and flux ddt correction
surfaceScalarField pRC
(
  "pRC",
  fvc::interpolate(rAU)*(fvc::interpolate(fvc::grad(p_rgh)) & mesh.Sf())
 +MRF.zeroFilter(fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)) 
);

fvScalarMatrix p_rghEqn
(
  fvm::laplacian(rAU, p_rgh)
 ==
  fvc::div(pRC)
);

p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

// pEqn
cps->insertEquation
(
 p_rgh.name(),
 p_rgh.name(),
 p_rghEqn
);

// div(U)
cps->insertEquation
(
 p_rgh.name(),
 U.name(),
 fvmb::div(-1., U)
);
