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


template<>
tmp<GeometricField<scalar, fvsPatchField, surfaceMesh> >
gaussDefCmpwConvectionScheme<scalar>::interpolate
(
    const surfaceScalarField& faceFlux,
    const GeometricField<scalar, fvPatchField, volMesh>& vf
) const
{    

   return phifDefCs(faceFlux, Foam::pos(faceFlux), vf, false);

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
    fvMatrix<Type>& fvm = tfvm();

    if (scheme_ == "none") {tfvm().diag() = 0; return tfvm;}
  
    surfaceScalarField upw = Foam::pos(faceFlux);

    const labelUList& own = mesh.owner();
    const labelUList& neig = mesh.neighbour();

    fvm.lower() = -upw.internalField()*faceFlux.internalField(); // Negative because the contribution from neigb is (-)
    fvm.upper() = (1.0-upw.internalField() )* faceFlux.internalField();
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
    
    GeometricField<Type, fvsPatchField, surfaceMesh>& phifDC = tphifDC();
 
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

    for (int cmp=0; cmp<Type::nComponents; cmp++)
     {
        // Internal Field
      
        forAll(neig, f)
         {   

            souT[own[f]].component(cmp) += phifDC[f].component(cmp) * faceFlux[f];
            souT[neig[f]].component(cmp) -= phifDC[f].component(cmp) * faceFlux[f];
     
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

		   souT[fC[faceI]].component(cmp) += phifDCp[faceI].component(cmp) * patchPhi[faceI]; // Only contributes once (to owner cell)
                
                 } //ForAll BC faces
                      
             }// If patch is coupled        
         
         } //For all patches

     } // For all components

  fvm.source() = -souT;  
   
  return tfvm;
}

template<>
tmp<fvMatrix<scalar> >
gaussDefCmpwConvectionScheme<scalar>::fvmDiv
(
    const surfaceScalarField& faceFlux,
    const GeometricField<scalar, fvPatchField, volMesh>& vf
) const
{

    const fvMesh& mesh = vf.mesh();
 
    tmp<fvMatrix<scalar> > tfvm
    (
        new fvMatrix<scalar>
        (
            vf,
            faceFlux.dimensions()*vf.dimensions()
        )
    );
    fvMatrix<scalar>& fvm = tfvm();

    if (scheme_ == "none") {tfvm().diag() = 0; return tfvm;}
  
    surfaceScalarField upw = Foam::pos(faceFlux);

    const labelUList& own = mesh.owner();
    const labelUList& neig = mesh.neighbour();

    fvm.lower() = -upw.internalField()*faceFlux.internalField(); // Negative because the contribution from neigb is (-)
    fvm.upper() = (1.0-upw.internalField() )* faceFlux.internalField();
    fvm.negSumDiag();

    forAll(vf.boundaryField(), patchI)
    {
        const fvPatchField<scalar>& psf = vf.boundaryField()[patchI];
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
  
    tmp<GeometricField<scalar, fvsPatchField, surfaceMesh> >  tphifDC( phifDefCs(faceFlux, upw, vf, true));

    GeometricField<scalar, fvsPatchField, surfaceMesh>& phifDC = tphifDC();
 
    GeometricField<scalar, fvPatchField, volMesh> souT
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
            dimensioned<scalar>
            (
                "0",
                vf.dimensions()*faceFlux.dimensions(),
                pTraits<scalar>::zero
            )
        );

        // Internal Field
      
        forAll(neig, f)
         {   

            souT[own[f]]  += phifDC[f] * faceFlux[f];
            souT[neig[f]] -= phifDC[f] * faceFlux[f];
     
         }

        // Boundary Field

        forAll(vf.boundaryField(), patchI)
         {
           
           if (vf.boundaryField()[patchI].coupled())
             {               
              const fvsPatchField<scalar>& phifDCp = phifDC.boundaryField()[patchI];
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

    tConvection().rename
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
    GeometricField<Type, fvsPatchField, surfaceMesh>& sf = tsf();

    if (scheme_ == "none") {sf = sf*0; return tsf;} // Force to be 0

    scalar swit(0.); // 0 (false) to return full phi_HRS; 1 (true) to return only the (phi_HRS - phi_upw) difference/correction
    if (onlyDCphi) {swit = 1.;} 

    Field<Type>& sfi = sf.internalField();

    const labelUList& own = mesh.owner();
    const labelUList& neig = mesh.neighbour();

    scalarList alphaL; scalarList betaL; scalarList bL;
    lims(alphaL,betaL,bL,faceFlux,vf);
    
    scalar phitc(0.0); scalar alpha(1.0);  scalar beta(0.0); 
    
    const vectorField& C = mesh.C(); 

    for (int cmp=0; cmp<Type::nComponents; cmp++)
     {
       volVectorField gradcmp = fvc::grad(vf.component(cmp));
      
       forAll(sfi, f)
        {
           phitc = 1.0 -
           (
           ( vf[neig[f]].component(cmp) - vf[own[f]].component(cmp) )/
           (  2.0 * ( ( gradcmp[own[f]]*upw[f] + (1.0 - upw[f]) * gradcmp[neig[f]] ) & (C[neig[f]]-C[own[f]]) ) + 1e-18 )  		
           );

           if (phitc <= 0 || phitc >= 1 ) { alpha = 1.; beta = 0.;}
           else if (phitc < bL[0] ) { alpha = alphaL[0];  beta = betaL[0]; }
           else if (phitc < bL[1] ) { alpha = alphaL[1];  beta = betaL[1]; }
           else { alpha = alphaL[2];  beta = betaL[2]; }
 
           sfi[f].component(cmp) =
              (1.0 - alpha - beta) * ( vf[neig[f]].component(cmp) - 2.0 * (  gradcmp[own[f]]  & (C[neig[f]]-C[own[f]]) ) ) * upw[f] 
            + (1.0 - alpha - beta) * ( vf[own[f]].component(cmp) + 2.0 * (  gradcmp[neig[f]]  & (C[neig[f]]-C[own[f]]) ) ) * (1.0-upw[f])
            + ( (alpha - 1.0*swit) * upw[f] + beta * ( 1.0 - upw[f] ) ) * vf[own[f]].component(cmp) 
            + ( beta * upw[f] + (alpha - 1.0*swit) * ( 1.0 - upw[f] ) ) * vf[neig[f]].component(cmp);
           
        }

   // Boundaries

        forAll(vf.boundaryField(), patchI)
         {
         
           fvsPatchField<Type>& sfp = sf.boundaryField()[patchI];
           
           if (vf.boundaryField()[patchI].coupled())
             {
              
              const Field<Type> pPF = vf.boundaryField()[patchI].patchInternalField();
              const Field<Type> pNF = vf.boundaryField()[patchI].patchNeighbourField();
              const vectorField gradcmpNF = gradcmp.boundaryField()[patchI].patchNeighbourField(); 
              const vectorField gradcmpPF = gradcmp.boundaryField()[patchI].patchInternalField(); 
              const scalarField& upwP = upw.boundaryField()[patchI];
              const vectorField& deltasP = vf.boundaryField()[patchI].patch().delta();
              
              forAll(sfp, faceI)
                {
 
                  phitc = 1.0 -
                   (
                     ( pNF[faceI].component(cmp) - pPF[faceI].component(cmp) )/
                     (  2.0 * ( ( gradcmpPF [faceI]*upwP[faceI] + (1.0 - upwP[faceI]) * gradcmpNF[faceI]) & deltasP[faceI] ) + 1e-18 )  		
                   );

                  if (phitc <= 0 || phitc >= 1 ) { alpha = 1.; beta = 0.;}
          	  else if (phitc < bL[0] ) { alpha = alphaL[0];  beta = betaL[0]; }
          	  else if (phitc < bL[1] ) { alpha = alphaL[1];  beta = betaL[1]; }
           	  else { alpha = alphaL[2];  beta = betaL[2]; }
    
                 sfp[faceI].component(cmp) = 
                   (1.0 - alpha - beta) * (pNF[faceI].component(cmp) - 2.0 * (  gradcmpPF[faceI]  & deltasP[faceI] ) ) * upwP[faceI]
            	 + (1.0 - alpha - beta) * (pPF[faceI].component(cmp) + 2.0 * (  gradcmpNF[faceI]  & deltasP[faceI] ) ) * (1.0-upwP[faceI])
            	 + ( (alpha - 1.0*swit) * upwP[faceI] + beta * ( 1.0 - upwP[faceI] ) ) * pPF[faceI].component(cmp) 
            	 + ( beta * upwP[faceI] + (alpha - 1.0*swit) * ( 1.0 - upwP[faceI] ) ) * pNF[faceI].component(cmp);

                } //ForAll BC faces
                      
              } // If patch is coupled
              else 
              {             
                 sfp.component(cmp) = vf.boundaryField()[patchI].component(cmp);
              }        
         
           } //For all patches

 
    } // For all components

    return tsf; 
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
gaussDefCmpwConvectionScheme<Type>::phifDefCs
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
    GeometricField<Type, fvsPatchField, surfaceMesh>& sf = tsf();

    if (scheme_ == "none") {sf = sf*0; return tsf;} // Force to be 0 

    scalar swit(0.); // 0 (false) to return full phi_HRS; 1 (true) to return only the (phi_HRS - phi_upw) difference/correction
    if (onlyDCphi) {swit = 1.;}   

    Field<Type>& sfi = sf.internalField();

    const labelUList& own = mesh.owner();
    const labelUList& neig = mesh.neighbour();

    scalarList alphaL; scalarList betaL; scalarList bL;
    lims(alphaL,betaL,bL,faceFlux,vf);
    
    scalar phitc(0.0); scalar alpha(1.0);  scalar beta(0.0); 
    
    const vectorField& C = mesh.C(); 
  
    volVectorField gradcmp = fvc::grad(vf);
      
    forAll(sfi, f)
        {
      
           phitc = 1.0 -
           (
           ( vf[neig[f]] - vf[own[f]] )/
           (  2.0 * ( ( gradcmp[own[f]]*upw[f] + (1.0 - upw[f]) * gradcmp[neig[f]] ) & (C[neig[f]]-C[own[f]]) ) + 1e-18 )  		
           );

           if (phitc <= 0 || phitc >= 1 ) { alpha = 1.; beta = 0.;}
           else if (phitc < bL[0] ) { alpha = alphaL[0];  beta = betaL[0]; }
           else if (phitc < bL[1] ) { alpha = alphaL[1];  beta = betaL[1]; }
           else { alpha = alphaL[2];  beta = betaL[2]; }
 
           sfi[f] =
              (1.0 - alpha - beta) * ( vf[neig[f]]  - 2.0 * (  gradcmp[own[f]]  & (C[neig[f]]-C[own[f]]) ) ) * upw[f] 
            + (1.0 - alpha - beta) * ( vf[own[f]]  + 2.0 * (  gradcmp[neig[f]]  & (C[neig[f]]-C[own[f]]) ) ) * (1.0-upw[f])
            + ( (alpha - 1.0*swit) * upw[f] + beta * ( 1.0 - upw[f] ) ) * vf[own[f]]  
            + ( beta * upw[f] + (alpha - 1.0*swit) * ( 1.0 - upw[f] ) ) * vf[neig[f]] ;
        }

   // Boundaries

        forAll(vf.boundaryField(), patchI)
         {
         
           fvsPatchField<Type>& sfp = sf.boundaryField()[patchI];
           
           if (vf.boundaryField()[patchI].coupled())
             {
             
              const Field<Type> pPF = vf.boundaryField()[patchI].patchInternalField();
              const Field<Type> pNF = vf.boundaryField()[patchI].patchNeighbourField();
              const vectorField gradcmpNF = gradcmp.boundaryField()[patchI].patchNeighbourField(); 
              const vectorField gradcmpPF = gradcmp.boundaryField()[patchI].patchInternalField(); 
              const scalarField& upwP = upw.boundaryField()[patchI];
              const vectorField& deltasP = vf.boundaryField()[patchI].patch().delta();

              forAll(sfp, faceI)
                {
 
                  phitc = 1.0 -
                   (
                     ( pNF[faceI] - pPF[faceI])/
                     (  2.0 * ( ( gradcmpPF [faceI]*upwP[faceI] + (1.0 - upwP[faceI]) * gradcmpNF[faceI]) & deltasP[faceI] ) + 1e-18 )  		
                   );

                  if (phitc <= 0 || phitc >= 1 ) { alpha = 1.; beta = 0.;}
          	  else if (phitc < bL[0] ) { alpha = alphaL[0]; beta = betaL[0]; }
          	  else if (phitc < bL[1] ) { alpha = alphaL[1]; beta = betaL[1]; }
           	  else { alpha = alphaL[2];  beta = betaL[2]; }
    
                 sfp[faceI]  =
                   (1.0 - alpha - beta) * (pNF[faceI]  - 2.0 * (  gradcmpPF[faceI]  & deltasP[faceI] ) ) * upwP[faceI]
            	 + (1.0 - alpha - beta) * (pPF[faceI]  + 2.0 * (  gradcmpNF[faceI]  & deltasP[faceI] ) ) * (1.0-upwP[faceI])
            	 + ( (alpha - 1.0*swit) * upwP[faceI] + beta * ( 1.0 - upwP[faceI] ) ) * pPF[faceI]  
            	 + ( beta * upwP[faceI] + (alpha - 1.0*swit) * ( 1.0 - upwP[faceI] ) ) * pNF[faceI];
         
                } //ForAll BC faces
                      
              } // If patch is coupled
              else 
              {             
                 sfp = vf.boundaryField()[patchI];
              }        
         
           } //For all patches

    return tsf; 
}




template<class Type>
void gaussDefCmpwConvectionScheme<Type>::lims
(
   scalarList& alpha,
   scalarList& beta,
   scalarList& bounds,
   const surfaceScalarField& faceFlux,
   const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    if (scheme_ == "upwind")
     {
       alpha.append(1.); 
       beta.append(0.); 
       bounds.append(1.); 
     } 
     else if (scheme_ == "cubista")
     {
       alpha.append(7./4.); alpha.append(3./4.); alpha.append(1./4.); 
       beta.append(0.); beta.append(3./8.); beta.append(3./4.); 
       bounds.append(3./8.); bounds.append(3./4.);
     } 
     else if (scheme_ == "minmod")
     {
       alpha.append(1.5); alpha.append(.5); alpha.append(.5); 
       beta.append(0.); beta.append(.5); beta.append(.5);   
       bounds.append(.5); bounds.append(1.); 
     }  
     else if (scheme_ == "smart")
     {
       alpha.append(3.); alpha.append(3./4.); alpha.append(0.); 
       beta.append(0.); beta.append(3./8.); beta.append(1.); 
       bounds.append(1./6.); bounds.append(5./6.);
     } 
     else if (scheme_ == "waceb")
     {
       alpha.append(2.); alpha.append(3./4.); alpha.append(0.); 
       beta.append(0.); beta.append(3./8.); beta.append(1.); 
       bounds.append(3./10.); bounds.append(5./6.);
     } 
     else if (scheme_ == "superbee")
     {
       alpha.append(0.5); alpha.append(1.5); alpha.append(0.); 
       beta.append(0.5); beta.append(0.); beta.append(1.); 
       bounds.append(1./2.); bounds.append(2./3.);
     } 
     else if (scheme_ == "none")
     {
       // Do nothing
     } 
     else
     {
         FatalErrorIn("gaussDefCmpwConvectionScheme<Type>::lims()\n")
            << "\nError in div("<< faceFlux.name() << ", " << vf.name() << ")\n"
            << "\nThe deferred limited scheme is not specified or does not exist.\n"
            << "\nValid schemes are:\n"
            << "\n. upwind" <<"\n. cubista" << "\n. minmod" << "\n. smart"
            << "\n. waceb" << "\n. superbee" << "\n. none (no convection)" 
            << abort(FatalError);
 
     } 

}
    



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
