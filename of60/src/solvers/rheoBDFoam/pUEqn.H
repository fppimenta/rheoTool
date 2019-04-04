//------- U -----------//
tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(U) 
  + fvm::div(phi, U)
  + MRF.DDt(U)
 ==
    fvOptions(U)
  + constEq.divTau(U) 
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
 p.name(),
 fvmb::grad(1., p)
);

// Insert div(tau) if tau is to be solved coupled,
// otherwise do no nothing (div(tau) becomes inserted explicitly
// in UEqn) 
constEq.divTauImplCoupled();

// Insert -rhoE*grad(psi_) if ci-psi are to be solved coupled,
// otherwise do no nothing (-rhoE*grad(psi_) becomes inserted
// explicitly in UEqn). 
if (solveElecM)
{
  elecM.FeImplCoupled(constEq.rho());
}

//------- p -----------//

surfaceScalarField pRC
(
  "pRC",
  fvc::interpolate(rAU)*(fvc::interpolate(fvc::grad(p)) & mesh.Sf())
);

fvScalarMatrix pEqn
(
  fvm::laplacian(rAU, p)
 ==
  fvc::div(pRC)
);

pEqn.setReference(pRefCell, pRefValue);

// pEqn
cps->insertEquation
(
 p.name(),
 p.name(),
 pEqn
);

// div(U)
cps->insertEquation
(
 p.name(),
 U.name(),
 fvmb::div(-1., U)
);


