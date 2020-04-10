/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "gaussDefCmpwConvectionScheme.H"
#include "fvcSurfaceIntegrate.H"
#include "fvMatrices.H"
#include "fvCFD.H"
#include "coupledFvPatchFields.H"
#include "limiters.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
gaussDefCmpwConvectionScheme<Type>::interpolate
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{

    return phifDefC(faceFlux, Foam::pos(faceFlux), vf, false); 
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
gaussDefCmpwConvectionScheme<Type>::flux
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    return faceFlux*interpolate(faceFlux, vf);
}


template<class Type>
tmp<fvMatrix<Type> >
gaussDefCmpwConvectionScheme<Type>::fvmDiv
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{

    const fvMesh& mesh = vf.mesh();
 
    tmp<fvMatrix<Type> > tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            faceFlux.dimensions()*vf.dimensions()
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    if (scheme_ == "none") {fvm.diag() = 0; return tfvm;}
  
    surfaceScalarField upw = Foam::pos(faceFlux);

    const labelUList& own = mesh.owner();
    const labelUList& neig = mesh.neighbour();

    fvm.lower() = -upw.primitiveField()*faceFlux.primitiveField(); // Negative because the contribution from neigb is (-)
    fvm.upper() = (1.0-upw.primitiveField() )* faceFlux.primitiveField();
    fvm.negSumDiag();

    forAll(vf.boundaryField(), patchI)
    {
        const fvPatchField<Type>& psf = vf.boundaryField()[patchI];
        const fvsPatchScalarField& patchPhi = faceFlux.boundaryField()[patchI];
         
        scalarField plim = scalarField( psf.size(), pTraits<scalar>::one );  
   
        if (vf.boundaryField()[patchI].coupled())
          {
              plim = upw.boundaryField()[patchI];     // This will give the same upwind weigthing as if it was an internalFace           
          } 

        fvm.internalCoeffs()[patchI] = patchPhi*psf.valueInternalCoeffs(plim);
        fvm.boundaryCoeffs()[patchI] = -patchPhi*psf.valueBoundaryCoeffs(plim);
    }

    if (scheme_ == "upwind") {return tfvm;}

    // Explicit correction
  
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >  tphifDC( phifDefC(faceFlux, upw, vf, true));
         
    GeometricField<Type, fvsPatchField, surfaceMesh>& phifDC = tphifDC.ref();
 
    GeometricField<Type, fvPatchField, volMesh> souT
    (
    IOobject
    (
       "souT",
        vf.instance(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensioned<Type>
    (
        "0",
        vf.dimensions()*faceFlux.dimensions(),
        pTraits<Type>::zero
    )
    );

    
    // Internal Field 

    forAll(neig, f)
    {            
     souT[own[f]] += phifDC[f] * faceFlux[f];
     souT[neig[f]] -= phifDC[f] * faceFlux[f];
    }

    // Boundary Field

    forAll(vf.boundaryField(), patchI)
    {
   
      if (vf.boundaryField()[patchI].coupled())
        {               
          const fvsPatchField<Type>& phifDCp = phifDC.boundaryField()[patchI];
          const fvsPatchScalarField& patchPhi = faceFlux.boundaryField()[patchI];
          const labelUList& fC= mesh.boundaryMesh()[patchI].faceCells();
       
          forAll(patchPhi, faceI)
           {
	       souT[fC[faceI]] += phifDCp[faceI] * patchPhi[faceI]; // Only contributes once (to owner cell)
       
           } //ForAll BC faces
              
        }// If patch is coupled        
 
    } //For all patches
 
   fvm.source() = -souT;
   
   return tfvm;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
gaussDefCmpwConvectionScheme<Type>::fvcDiv
(
    const surfaceScalarField& faceFlux,
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    tmp<GeometricField<Type, fvPatchField, volMesh> > tConvection
    (
        fvc::surfaceIntegrate(flux(faceFlux, vf))
    );

    tConvection.ref().rename
    (
        "convection(" + faceFlux.name() + ',' + vf.name() + ')'
    );

    return tConvection;
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
gaussDefCmpwConvectionScheme<Type>::phifDefC
(
    const surfaceScalarField& faceFlux,
    const surfaceScalarField& upw,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const bool& onlyDCphi
) const
{
 
    const fvMesh& mesh = vf.mesh();

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsf
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "interpolateDef("+vf.name()+')',
                vf.instance(),
                vf.db()
            ),
            mesh,
            vf.dimensions()
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& sf = tsf.ref();

    if (scheme_ == "none") {sf = sf*0.; return tsf;} 

    scalar swit(0.); // 0 (false) to return full phi_HRS; 1 (true) to return only the (phi_HRS - phi_upw) difference/correction
    if (onlyDCphi) {swit = 1.;} 

    Field<Type>& sfi = sf.primitiveFieldRef();

    const labelUList& own = mesh.owner();
    const labelUList& neig = mesh.neighbour();

    scalarList alphaL; scalarList betaL; scalarList bL;
    lims(alphaL,betaL,bL,faceFlux,vf);
    
    scalar phitc(0.0); scalar alpha(1.0);  scalar beta(0.0); 
  
    const vectorField& C = mesh.C();


    for (int cmp=0; cmp < pTraits<Type>::nComponents; cmp++)
     {

   // Internal field
   
       tmp<volScalarField> tvfc = vf.component(cmp);
       
       // The component of a volScalarField is the volScalarField itself and
       // it is not returned as true tmp. Simple hacking to render it so (of4? exclusive). 
       if ( pTraits<Type>::nComponents == 1)      
           tvfc = tmp<volScalarField>(vf.component(cmp)*1.);
                      
       volScalarField& vfc = tvfc.ref();
       
       volVectorField gradcmp = fvc::grad(vfc); 
       Field<scalar>  tf(sfi.size(), pTraits<scalar>::zero);
      
       forAll (sfi, f)
        {
           phitc = 1.0 -
           (
           ( vfc[neig[f]] - vfc[own[f]] )/
           (  2.0 * ( ( gradcmp[own[f]]*upw[f] + (1.0 - upw[f]) * gradcmp[neig[f]] ) & (C[neig[f]]-C[own[f]])) + 1e-18 )  		
           );

           if (phitc <= 0. || phitc >= 1. ) { alpha = 1.; beta = 0.;}
           else if (phitc < bL[0] ) { alpha = alphaL[0];  beta = betaL[0]; }
           else if (phitc < bL[1] ) { alpha = alphaL[1];  beta = betaL[1]; }
           else { alpha = alphaL[2];  beta = betaL[2]; }
 
           tf[f] =
              (1.0 - alpha - beta) * ( vfc[neig[f]]  - 2.0 * (  gradcmp[own[f]]  & (C[neig[f]]-C[own[f]])) ) * upw[f] 
            + (1.0 - alpha - beta) * ( vfc[own[f]]  + 2.0 * (  gradcmp[neig[f]]  & (C[neig[f]]-C[own[f]])) ) * (1.0-upw[f])
            + ( (alpha - 1.0*swit) * upw[f] + beta * ( 1.0 - upw[f] ) ) * vfc[own[f]] 
            + ( beta * upw[f] + (alpha - 1.0*swit) * ( 1.0 - upw[f] ) ) * vfc[neig[f]];
           
        }
        
       sfi.replace(cmp, tf);
    
   // Boundaries

       typename GeometricField<Type, fvsPatchField, surfaceMesh>::
       Boundary& sfb = sf.boundaryFieldRef();

        forAll(vf.boundaryField(), patchI)
         {
         
           Field<Type>& sfp = sfb[patchI];
           Field<scalar> sfpt(sfp.size(), pTraits<scalar>::zero);
           
           if (vf.boundaryField()[patchI].coupled())
             {

              const Field<scalar> pPF = vfc.boundaryField()[patchI].patchInternalField();
              const Field<scalar> pNF = vfc.boundaryField()[patchI].patchNeighbourField();
              const vectorField gradcmpNF = gradcmp.boundaryField()[patchI].patchNeighbourField(); 
              const vectorField gradcmpPF = gradcmp.boundaryField()[patchI].patchInternalField(); 
              const scalarField& upwP = upw.boundaryField()[patchI];
              vectorField deltasP(vf.boundaryField()[patchI].patch().delta());
              
              forAll(sfp, faceI)
                {
 
                  phitc = 1.0 -
                   (
                     ( pNF[faceI] - pPF[faceI] )/
                     (  2.0 * ( ( gradcmpPF [faceI]*upwP[faceI] + (1.0 - upwP[faceI]) * gradcmpNF[faceI]) & deltasP[faceI] ) + 1e-18 )  		
                   );

                  if (phitc <= 0. || phitc >= 1. ) { alpha = 1.; beta = 0.;}
          	  else if (phitc < bL[0] ) { alpha = alphaL[0];  beta = betaL[0]; }
          	  else if (phitc < bL[1] ) { alpha = alphaL[1];  beta = betaL[1]; }
           	  else { alpha = alphaL[2];  beta = betaL[2]; }
    
                  sfpt[faceI] = 
                    (1.0 - alpha - beta) * (pNF[faceI] - 2.0 * (  gradcmpPF[faceI]  & deltasP[faceI] ) ) * upwP[faceI]
            	  + (1.0 - alpha - beta) * (pPF[faceI] + 2.0 * (  gradcmpNF[faceI]  & deltasP[faceI] ) ) * (1.0-upwP[faceI])
            	  + ( (alpha - 1.0*swit) * upwP[faceI] + beta * ( 1.0 - upwP[faceI] ) ) * pPF[faceI] 
            	  + ( beta * upwP[faceI] + (alpha - 1.0*swit) * ( 1.0 - upwP[faceI] ) ) * pNF[faceI];

                } //ForAll BC faces
                
                sfp.replace(cmp, sfpt);
          
              } // If patch is coupled
            else 
              {             
                 sfp.replace(cmp, vf.boundaryField()[patchI].component(cmp));
              }        
         
           } //For all patches

 
    } // For all components

    return tsf; 
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
