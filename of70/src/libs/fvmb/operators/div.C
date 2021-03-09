/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "div.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
/*
Note: List LMatrix.rel contains the row-col-cmp indices. Ie, LMatrix.field(cmp) 
is to be inserted in (row, col) of the block matrix. Both row/col assume that the
first valid component of the fields involved is 0 (fields ordering in the block 
matrix is accounted for later inside the coupled solver). The coupled solver will
only insert the equations specified in LMatrix.rel, even though the operators 
are built for all the 3 geometric dimensions.
*/

namespace fvmb
{
 
// symmTensor-vector
// Operator: c*[grad(U) + T(grad(U))]
 

// vector-symmTensor

tmp< LMatrix<tensor> >  div
(
   const dimensionedScalar c,
   const GeometricField<symmTensor, fvPatchField, volMesh>& vf
)
{
   return div(c.value(), vf);
}

tmp< LMatrix<tensor> >  div
(
   const scalar c,
   const GeometricField<symmTensor, fvPatchField, volMesh>& vf
)
{
  const fvMesh& mesh = vf.mesh();
  
  tmp< LMatrix<tensor> > tfvm
  (
     LMatrix<tensor>::New(mesh)
  );
  LMatrix<tensor>& fvm = tfvm.ref();
 
  tensorField& upper = fvm.upper();
  tensorField& lower = fvm.lower();
  tensorField& diag = fvm.diag();
  List<labelList>& rel = fvm.rowColList();
   
  const labelUList& own = mesh.owner();
  const labelUList& neig = mesh.neighbour();
  const vectorField Sf(mesh.faceAreas()*c); // Multiply constant
  
  surfaceScalarField weights = vf.mesh().surfaceInterpolation::weights();

  forAll(own, facei)
  {    
    scalar w = weights[facei];
    
    tensor wo = w*tensor(Sf[facei],Sf[facei],Sf[facei]);
    tensor wn = (1.-w)*tensor(Sf[facei],Sf[facei],Sf[facei]);
    
    // For owner
    diag[own[facei]] += wo;
    upper[facei] = wn;
    
    // For neig
    diag[neig[facei]] -= wn;   
    lower[facei] = -wo; 
  } 
  
  forAll(vf.boundaryField(), patchi)
  {
    const fvPatchField<symmTensor>& pvf = vf.boundaryField()[patchi];
    const fvsPatchField<scalar>& bw = weights.boundaryField()[patchi];
    const vectorField Sfb(pvf.patch().Sf()*c); // Multiply constant
    tensorField& internalCoeffsFvm(fvm.internalCoeffs()[patchi]);
    tensorField& boundaryCoeffsFvm(fvm.boundaryCoeffs()[patchi]);
        
    // bw refer to the f-C fractional distance on the other processor. The
    // cell on each processor is always the owner, so bw still corresponds 
    // to the fractional distance of the neighbour.
    if (pvf.coupled())
    {         
      forAll(pvf, facei)
      {
        tensor SfbTensor = tensor(Sfb[facei],Sfb[facei],Sfb[facei]);
        
        // To procowner
        internalCoeffsFvm[facei] = bw[facei] * SfbTensor; 
        
        // To procneigh
        boundaryCoeffsFvm[facei] = -(1. - bw[facei]) * SfbTensor; // the (-) is because it would go to source
      }
    }
    else
    {
      const symmTensorField internalCoeffs(pvf.valueInternalCoeffs(bw));
      const symmTensorField boundaryCoeffs(pvf.valueBoundaryCoeffs(bw));
     
      forAll(pvf, facei)
      {
        tensor SfbTensor = tensor(Sfb[facei],Sfb[facei],Sfb[facei]);
        internalCoeffsFvm[facei] = cmptMultiply(SfbTensor, tensor(internalCoeffs[facei])); 
        boundaryCoeffsFvm[facei] = cmptMultiply(-SfbTensor, tensor(boundaryCoeffs[facei]));
      }          
    }
  }
  
  // Fill the list
  
  //-- Tensor
  typename pTraits<symmTensor>::labelType validComponentsT
  (
    mesh.template validComponents<symmTensor>()
  );
  
  int i = 0;
  labelList cT(6,-1);
  forAll(cT, cmpt)
  {
    if (component(validComponentsT, cmpt) == -1) continue; 
     
    cT[cmpt]= i;
    i++;
  }
  
  //-- Vector
  typename pTraits<vector>::labelType validComponentsV
  (
    mesh.template validComponents<vector>()
  );
  
  i = 0;
  labelList cV(3,-1);
  forAll(cV, cmpt)
  {
    if (component(validComponentsV, cmpt) == -1) continue; 
    
    cV[cmpt] = i;     
    i++;
  }
 
  // Fill the row-col-cmpt list
  List<labelList> VT({
  {cV[0],cT[0]}, {cV[0],cT[1]}, {cV[0],cT[2]},
  {cV[1],cT[1]}, {cV[1],cT[3]}, {cV[1],cT[4]},
  {cV[2],cT[2]}, {cV[2],cT[4]}, {cV[2],cT[5]}
  });
   
  forAll(VT, i)
  {
    if (VT[i][0] == -1 || VT[i][1] == -1 ) continue;
    
    rel.append( {VT[i][0], VT[i][1], i} );
  }
   
  return tfvm;    
} 

// scalar-vector
tmp< LMatrix<vector> >  div
(
   const dimensionedScalar c,
   const GeometricField<vector, fvPatchField, volMesh>& vf
)
{
  return div(c.value(), vf);
}

tmp< LMatrix<vector> >  div
(
   const scalar c,
   const GeometricField<vector, fvPatchField, volMesh>& vf
)
{
  const fvMesh& mesh = vf.mesh();
  
  tmp<LMatrix<vector> > tfvm
  (
    LMatrix<vector>::New(mesh)
  );
  LMatrix<vector>& fvm = tfvm.ref();
 
  vectorField& upper = fvm.upper();
  vectorField& lower = fvm.lower();
  vectorField& diag = fvm.diag();
  List<labelList>& rel = fvm.rowColList();
   
  const labelUList& own = mesh.owner();
  const labelUList& neig = mesh.neighbour();
  const vectorField Sf(mesh.faceAreas()*c); // Multiply constant
  
  surfaceScalarField weights = vf.mesh().surfaceInterpolation::weights();

  forAll(own, facei)
  {    
    scalar w = weights[facei];
    
    // For the owner[facei] cell 
    diag[own[facei]] += w*Sf[facei];
    upper[facei] = (1. - w)*Sf[facei];
    
    // For the neig[facei] cell 
    diag[neig[facei]] -= (1. - w)*Sf[facei];
    lower[facei] = -w*Sf[facei];    
  } 
  
  forAll(vf.boundaryField(), patchi)
  {
    const fvPatchField<vector>& pvf = vf.boundaryField()[patchi];
    const fvsPatchField<scalar>& bw = weights.boundaryField()[patchi];
    const vectorField Sfb(pvf.patch().Sf()*c); // Multiply constant
    
    // bw refer to the f-C fractional distance on the other processor. The
    // cell on each processor is always the owner, so bw still corresponds 
    // to the fractional distance of the neighbour.
    if (pvf.coupled())
    {     
       // To procowner
       fvm.internalCoeffs()[patchi] =  bw * Sfb;
      
       // To procneigh
       fvm.boundaryCoeffs()[patchi] = -(1. - bw) * Sfb; // the (-) is because it would go to source
    }
    else
    {
      const vectorField internalCoeffs(pvf.valueInternalCoeffs(bw));
      const vectorField boundaryCoeffs(pvf.valueBoundaryCoeffs(bw));
      
      fvm.internalCoeffs()[patchi] = cmptMultiply(internalCoeffs, Sfb);
        
      fvm.boundaryCoeffs()[patchi] = cmptMultiply(-boundaryCoeffs, Sfb);
    }
  }
  
  // Fill the row-col-cmpt list
  typename pTraits<vector>::labelType validComponents
  (
    mesh.template validComponents<vector>()
  );
 
  int i = 0;
  for (direction cmpt=0; cmpt<pTraits<vector>::nComponents; cmpt++)
  {
    if (component(validComponents, cmpt) == -1) continue; 
    
    // 0 is for the scalar var
    rel.append({0,i,cmpt});
    
    i++;
  }

  return tfvm;
    
}
 

} // namespace fvmb
} // namespace Foam
