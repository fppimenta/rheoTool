volScalarField::Boundary& bndPsiCorr = psiCorr.boundaryFieldRef();

// Reset to match psi_
if (i != 0)
 bndPsiCorr = psi_.boundaryField();

forAll(ci.boundaryField(), pI)
{
  if (isType<zeroIonicFluxFvPatchScalarField>(ci.boundaryField()[pI])) 
    bndPsiCorr.set(pI, new zeroGradientFvPatchField<scalar>(psi_.mesh().boundary()[pI], psiCorr.internalField()));
}
  
psiCorr.correctBoundaryConditions();
