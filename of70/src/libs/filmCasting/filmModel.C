/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "filmModel.H"
#include "freeSurfaceDisplacementPointPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(filmModel, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

filmModel::filmModel
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    IOdictionary
    (
        IOobject
        (
            "constitutiveProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    U_(U),
    isThermo_(readBool(subDict("filmProperties").lookup("isThermo"))),
    ht_(isThermo_ ? subDict("filmProperties").lookup("h") : dimensionedScalar("-", dimless, 0)),
    Cp_(isThermo_ ? subDict("filmProperties").lookup("Cp") : dimensionedScalar("-", dimless, 0)),
    Tair_(isThermo_ ? subDict("filmProperties").lookup("Tair") : dimensionedScalar("-", dimless, 0)),
    eqPtr_(NULL),
    TPtr_(NULL),
    h_
    (
        IOobject
        (
            "h",
            U.time().timeName(),
            U.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh()
    ),
    p_
    (
        IOobject
        (
            "p",
            U.time().timeName(),
            U.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U.mesh(),
        dimensionedScalar
        (
            "0",
            dimPressure,
            0
        ),
        extrapolatedCalculatedFvPatchField<scalar>::typeName
         
    ),
    isZeroTractionBC_(readBool(subDict("filmProperties").lookup("isTanZeroTraction"))),
    ppEnabled_(readBool(U.mesh().solutionDict().subDict("filmPostProcessing").lookup("enabled"))),
    ppFreq_(ppEnabled_ ? readInt(U.mesh().solutionDict().subDict("filmPostProcessing").lookup("writeFrequency")) : 0),
    ppCnt_(0),
    outW_(nullptr),
    outH_(nullptr)
{

 // Check which method to be used in fs (assuming the same method for all free-surfaces)
 pointVectorField::Boundary& pmXB = const_cast<pointVectorField&>
 (
   U.mesh().lookupObject<pointVectorField>("pointDisplacement")
 ).boundaryFieldRef();
    
 absFluxNeeded_ = false;
 forAll(pmXB, i)
 {
  if (isType<freeSurfaceDisplacementPointPatchVectorField>(pmXB[i]))
  {
    if (refCast<freeSurfaceDisplacementPointPatchVectorField>(pmXB[i]).method() == "streamline")
    {
       absFluxNeeded_ = true;                 
    }     
  }
 }

 // Thermo
 if (isThermo_)
 {
   TPtr_ = new  
   volScalarField
   (
     IOobject
     (
       "T",
       U.time().timeName(),
       U.mesh(),
       IOobject::MUST_READ,
       IOobject::AUTO_WRITE
     ),
     U.mesh()
   );
 }

 // Create constitutive equation only after T field exists
 eqPtr_ = constitutiveEq::New(word::null, U, phi, subDict("parameters")).ptr();
 
 // This is a hotfix for the multimode model, has this model doesn't call checkForStab on its own (only for each sub-model).
 // As such multimode doesn't have stabOption_ member defined, but we need it here in divTau(). Thus, force the function call.
 // For any model other than multimode, this call is just a second call and will do nothing.
 eqPtr_->checkForStab(subDict("parameters"));
 
 // Post-processing
 if (ppEnabled_)
 {
  ppDir_ =
  (
     Pstream::parRun() 
   ? U.mesh().time().path()/".."/"filmCastingPP"/U.mesh().time().timeName() 
   : U.mesh().time().path()/"filmCastingPP"/U.mesh().time().timeName()  
  );
   
  if (Pstream::master()) 
   mkDir(ppDir_);
  
  outW_.reset(new OFstream(ppDir_/"width.txt"));   
 }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void filmModel::checkFlux
(
  const volVectorField& U,
  const surfaceScalarField& phiRel
) const
{     
   pointVectorField::Boundary& pmXB = const_cast<pointVectorField&>
   (
     U.mesh().lookupObject<pointVectorField>("pointDisplacement")
   ).boundaryFieldRef();
    
   forAll(pmXB, i)
   {
    if (isType<freeSurfaceDisplacementPointPatchVectorField>(pmXB[i]))
    {           
        scalar absFlux = gSum(mag(U.boundaryField()[i] & U.mesh().Sf().boundaryField()[i]));
        scalar relFlux = gSum(mag(phiRel.boundaryField()[i]));
        Info << "|Flux| at " << U.mesh().boundary()[i].name() << " (abs/rel/%): " << absFlux 
             << tab <<  relFlux << tab << relFlux*100/(absFlux+1e-16) << endl;            
    }     
  }
  
  label outID = U.mesh().boundaryMesh().findPatchID("outlet");
  label inID = U.mesh().boundaryMesh().findPatchID("inlet");
  if (outID > -1)
  {
    scalar absFlux = gSum(mag( (U.boundaryField()[outID] * h_.boundaryField()[outID]) & U.mesh().Sf().boundaryField()[outID]));
    Info << "|Flux| at " << U.mesh().boundary()[outID].name() << " (abs): " << absFlux << endl; 
  }
  
  if (inID > -1)
  {
    scalar absFlux = gSum(mag( (U.boundaryField()[inID] * h_.boundaryField()[inID]) & U.mesh().Sf().boundaryField()[inID]));
    Info << "|Flux| at " << U.mesh().boundary()[inID].name() << " (abs): " << absFlux << endl; 
  }
}

void filmModel::updateFreeSurface
(
  bool isFinalIter,
  const volVectorField& U,
  const surfaceScalarField& phi
) const
{  
   pointVectorField::Boundary& pmXB = const_cast<pointVectorField&>
   (
     U.mesh().lookupObject<pointVectorField>("pointDisplacement")
   ).boundaryFieldRef();
    
   forAll(pmXB, i)
   {
    if (isType<freeSurfaceDisplacementPointPatchVectorField>(pmXB[i]))
    {
      if (absFluxNeeded_)
      {
        // phi should be absolute
        refCast<freeSurfaceDisplacementPointPatchVectorField>(pmXB[i]).updatePoints
        (
         isFinalIter, 
         phi.boundaryField()[i]
        );  
      }
      else
      {
        // phi should be relative
        refCast<freeSurfaceDisplacementPointPatchVectorField>(pmXB[i]).updatePoints
        (
         phi.boundaryField()[i]
        );  
      }                    
    }     
  }
}
 
void filmModel::forcedZeroFluxFreeSurface
(
  bool forced,
  surfaceScalarField& phiRel
) const
{
 if (forced)
 {
   pointVectorField::Boundary& pmXB = const_cast<pointVectorField&>
   (
     phiRel.mesh().lookupObject<pointVectorField>("pointDisplacement")
   ).boundaryFieldRef();
    
   forAll(pmXB, i)
   {
    if (isType<freeSurfaceDisplacementPointPatchVectorField>(pmXB[i]))
    {                     
      phiRel.boundaryFieldRef()[i] *= 0.;                      
    }     
   }
 }
}

void filmModel::updateHeight
(
  const surfaceScalarField& phi
) 
{
  // Solve Height equation
  fvScalarMatrix hEqn
  (
     fvm::ddt(h_)
   + fvm::div(phi, h_)
  );

  hEqn.relax(); 
  hEqn.solve();    
}


tmp<fvVectorMatrix> filmModel::divTau
(   
   const volVectorField& U
) 
{
 volTensorField gradU = fvc::grad(U);
 
 // Must use uncorrected gradU to update pressure
 updatePressure(gradU);
 
 // Make zz component of grad(U) 
 correctGradU(gradU);  
 
 // Note: in div(p*h*I), I is 3D, but no problem since u_z is not solved for 
 if (isGNF())
 { 
   return
   (
      fvm::laplacian(eqPtr_->eta()*h_, U, "laplacian(eta,U)")
    + fvc::div(eqPtr_->eta()*h_*dev2(T(gradU)), "div(eta*e*dev2(T(gradU)))") 
    - fvc::div(p_*h_*symmTensor::I, "div(pI)")
   ); 
 }
 else
 {  
   volSymmTensorField tauC(tau());
   
   // Here we ensure that tauXY contribution to Cauchy stress is null at the outlet 
   // if zero traction is selected. The null contribution from Newtonian stresses
   // is achieved by setting dv/dx = 0 and du/dy = 0.
   if (isZeroTractionBC_) 
   {
     label outID = U.mesh().boundaryMesh().findPatchID("outlet");
     volSymmTensorField::Boundary& tauCB = tauC.boundaryFieldRef();
     symmTensorField& tauCBI(tauCB[outID]);
     forAll(tauCBI, i)
     {
      tauCBI[i].xy() = 0.;
      tauCBI[i].xz() = 0.;
     }
   }
   
   volScalarField etaS(eqPtr_->etaSThermo());  
   switch (eqPtr_->stabOption_) 
   {      
     case constitutiveEq::stabOptions::soNone : // none   
     return
     (
          fvc::div(tauC*h_, "div(tau)")
        + fvm::laplacian(etaS*h_, U, "laplacian(eta,U)")
        + fvc::div(etaS*h_*dev2(T(gradU)), "div(eta*e*dev2(T(gradU)))")
        - fvc::div(p_*h_*symmTensor::I, "div(pI)")  
     );    
    
     // Note: for multimode modeling, etaP will be the sum of etaP from all modes.
     case constitutiveEq::stabOptions::soBSD : // BSD 
     { 
      volScalarField etaP(eqPtr_->etaPThermo());    
      return
      (
         fvc::div(tauC*h_, "div(tau)")
       - fvc::laplacian(etaP*h_, U, "laplacian(etaP,U)")
       + fvm::laplacian((etaP + etaS)*h_, U, "laplacian(eta,U)")
       + fvc::div(etaS*h_*dev2(T(gradU)), "div(eta*e*dev2(T(gradU)))")
       - fvc::div(p_*h_*symmTensor::I, "div(pI)")  
      );
     }
     
     case constitutiveEq::stabOptions::soCoupling : // coupling   
     {
      volScalarField etaP(eqPtr_->etaPThermo());      
      return
      (
         fvc::div(tauC*h_, "div(tau)")
       - fvc::div(etaP*h_*gradU, "div(etaP,grad(U))" )
       + fvm::laplacian((etaP + etaS)*h_, U, "laplacian(eta,U)")
       + fvc::div(etaS*h_*dev2(T(gradU)), "div(eta*e*dev2(T(gradU)))")
       - fvc::div(p_*h_*symmTensor::I, "div(pI)")  
      );     
     }
     
     default: // This will never happen
     return fvVectorMatrix(U, U.dimensions()); 
   }     
 }    
}


void filmModel::updatePressure
(
 const volTensorField& gradU
) 
{
 // Note: the pressure is made from xx and yy components of the newtonian
 // stress tensor + the zz component of the polymeric stress tensor. Therefore
 // grad(U_) must not be corrected, since the tr() would retrieve a wrong value.
 if (isGNF())
 {      
   p_ = -2.*eqPtr_->eta()*tr(gradU); 
 }
 else
 {
   p_ = -2.*eqPtr_->etaSThermo()*tr(gradU) + tau()().component(5); 
 }
}

void filmModel::correctGradU
(
  volTensorField& gradU
) const
{
  // Correct gradU: make zz component from continuity principle.
  // Both internal and boundaryField.
  gradU.replace(tensor::ZZ, -(gradU.component(tensor::XX) + gradU.component(tensor::YY)));
}


void filmModel::correctStresses
(
  const volVectorField& U
)
{
  volTensorField L = fvc::grad(U);
  correctGradU(L);
  eqPtr_->correct(nullptr,&L);
}

void filmModel::correctThermo
(
  const volVectorField& U,
  const surfaceScalarField& phi
)
{
  if (isThermo_)
  {
    volScalarField& T_ = TPtr_();
    
    fvScalarMatrix TEqn
    (
        fvm::ddt(T_) 
      + fvm::div(phi, T_) // Must use bounded scheme, since phi is not conservative
      ==
        2.*ht_*Tair_/(rho()*Cp_*h_) 
      - fvm::Sp(2.*ht_/(rho()*Cp_*h_), T_)
    );
      
    TEqn.relax();
    TEqn.solve(); 
  }
}


tmp<volSymmTensorField> filmModel::tauTotal
(
  const volVectorField& U
) const
{
 volTensorField L = fvc::grad(U);
 correctGradU(L);
   
 // Last term is only meaningfull for two-phase flows
 if (isGNF())
 {      
   return eqPtr_->eta()*(symm(L+L.T()) - (2./3.)*tr(L)*symmTensor::I); 
 }
 else
 {
   return tau() + eqPtr_->etaSThermo()*(symm(L+L.T()) - (2./3.)*tr(L)*symmTensor::I); 
 }
}


void filmModel::postProcess
(
  const volVectorField& U
) 
{
  if (!ppEnabled_)
   return;
  
  ppCnt_++;
  
  if (ppCnt_ > ppFreq_)
  {      
   //- Width at outlet
   pointVectorField::Boundary& pmXB = const_cast<pointVectorField&>
   (
     U.mesh().lookupObject<pointVectorField>("pointDisplacement")
   ).boundaryFieldRef();
    
   const pointField& pts = U.mesh().points();
   
   if (Pstream::master())
     outW_() << U.mesh().time().value() << tab;
   
   forAll(pmXB, i)
   { 
    scalar gW(-1e6);
    if (isType<freeSurfaceDisplacementPointPatchVectorField>(pmXB[i]))
    {                     
      const labelList& patchpts = pmXB[i].patch().meshPoints();
      scalar maxX(-1);
      forAll(patchpts, i)
      {
        if (pts[patchpts[i]].x() > maxX)
        {
           gW = pts[patchpts[i]].y();
           maxX = pts[patchpts[i]].x();
        }
      } 
      reduce(gW, maxOp<scalar>() ); 
          
      if (Pstream::master()) 
        outW_() << gW << tab;                       
    }  
   }
   
   if (Pstream::master())
     outW_() << endl;
  
   //- Thickness along the outlet
   if (Pstream::master())
     outH_.reset(new OFstream(ppDir_/"thickness.txt")); 
   
   label outID = U.mesh().boundaryMesh().findPatchID("outlet");
   
   if (Pstream::parRun())
   {
    // Share H
    List< List<scalar>  > allH(Pstream::nProcs());     
    allH[Pstream::myProcNo()] = h_.boundaryField()[outID];
    Pstream::gatherList(allH);           
    Pstream::scatterList(allH);
     
    List<scalar> sH = ListListOps::combine< List<scalar> >
    (
      allH,
      accessOp< List<scalar> >()
    );
    
    // Share Y
    List< List<scalar>  > allY(Pstream::nProcs());     
    allY[Pstream::myProcNo()] = h_.boundaryField()[outID].patch().Cf().component(1)();  
    Pstream::gatherList(allY);           
    Pstream::scatterList(allY);
     
    List<scalar> sY = ListListOps::combine< List<scalar> >
    (
      allY,
      accessOp< List<scalar> >()
    );
    
    if (Pstream::master())
    {
      forAll(sH, i)
      {
        outH_() << sY[i] << tab << sH[i] << endl; 
      }
    }
   }
   else
   {
      const fvPatchField<scalar>& hb = h_.boundaryField()[outID];
      const vectorField& Cf = h_.boundaryField()[outID].patch().Cf();
      forAll(hb, i)
      {
        outH_() << Cf[i].y() << tab << hb[i] << endl; 
      }
   } 
   
   //- Reset counter
   ppCnt_ = 0;
  } 
}


tmp<volSymmTensorField> filmModel::tau() const
{
    return eqPtr_->tau();
}
 

const dimensionedScalar filmModel::rho() const
{
    return eqPtr_->rho();
}

bool filmModel::isGNF() const
{
    return eqPtr_->isGNF();
}


bool filmModel::read()
{
    if (regIOobject::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
