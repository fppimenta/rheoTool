 volTensorField L(gradU == nullptr ? fvc::grad(U())() : *gradU);
  
 volTensorField B
 (
   IOobject
   (
     "B",
     U().time().timeName(),
     U().mesh(),
     IOobject::NO_READ,
     IOobject::NO_WRITE,
     false
   ),
   U().mesh(),
   dimensionedTensor
   (
     "zero",
     dimless/dimTime,
     pTraits<tensor>::zero
   ),
   extrapolatedCalculatedFvPatchField<tensor>::typeName
 );
 
 volTensorField   omega(B);
 
 #if defined(PTTLog_H) || defined(SaramitoLog_H)
 // volTensorField   M(eigVecs_.T() & ( L.T() - zeta_*symm(L) ) & eigVecs_);
  volTensorField M(innerP(eigVecs_,( L.T() - zeta_*symm(L) ) ,true));
 #else
 // volTensorField   M(eigVecs_.T() & L.T() & eigVecs_);
  volTensorField M(innerP(eigVecs_,L.T(),true));
 #endif
 
 decomposeGradU(M, eigVals_, eigVecs_, omega, B);
 
 dimensionedTensor Itensor("I",dimless,tensor::I);
