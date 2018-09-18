/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "linearExtrapolationFvPatchField.H"
#include "volFields.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
linearExtrapolationFvPatchField<Type>::linearExtrapolationFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    useReg_(false)
{}


template<class Type>
linearExtrapolationFvPatchField<Type>::linearExtrapolationFvPatchField
(
    const linearExtrapolationFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    useReg_(ptf.useReg_)
{}


template<class Type>
linearExtrapolationFvPatchField<Type>::linearExtrapolationFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict),
    useReg_(dict.lookupOrDefault<bool>("useRegression", false))
{}


template<class Type>
linearExtrapolationFvPatchField<Type>::linearExtrapolationFvPatchField
(
    const linearExtrapolationFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    useReg_(ptf.useReg_)
{}


template<class Type>
linearExtrapolationFvPatchField<Type>::linearExtrapolationFvPatchField
(
    const linearExtrapolationFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    useReg_(ptf.useReg_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void linearExtrapolationFvPatchField<Type>::updateCoeffs()
{
   if (this->updated())
   {
       return;
   }

   word varName( this->internalField().name() );

   const fvMesh& mesh = this->patch().boundaryMesh().mesh();  
 
   const GeometricField<Type, fvPatchField, volMesh>&
   var = this->db().objectRegistry::lookupObject< GeometricField<Type, fvPatchField, volMesh> >(varName);  

   Field<Type> varp = Field<Type>( this->patch().size(), pTraits<Type>::zero );

   if (!useReg_)
   {
     Field<scalar> varpI = Field<scalar>( this->patch().size(), pTraits<scalar>::zero );
     
     // Note: we are computing the gradient in all the mesh, but then we only use 
     // the values in the cell next to the patch. This is wasting time, but it can 
     // ensure a wide stencil. Tested a version we local computation of gaussGrad
     // but speedup was not signficative due to lookup for faces on patches.
     
     for (direction cmp = 0; cmp < pTraits<Type>::nComponents; cmp++)
      {          
        tmp<volScalarField> tvarI = var.component(cmp);
        if ( pTraits<Type>::nComponents == 1)      
           tvarI = tmp<volScalarField>(var.component(cmp)*1.);
           
        volScalarField& varI = tvarI.ref();
        volVectorField gradT = fvc::grad(varI, "linExtrapGrad");
    
        forAll(this->patch(), facei )        
          {
             vector r_face = this->patch().Cf()[facei];  

             label cellA = this->patch().faceCells()[facei];

             vector r_cellA = mesh.cellCentres()[cellA];  

             vector CtoF = (r_face - r_cellA);

             varpI[facei] = varI[cellA] + (gradT[cellA] & CtoF );
          }
       
        varp.replace(cmp, varpI);          
      }
   }
   else
   {       
     const labelUList& owner = mesh.owner();
     const labelUList& neighbour = mesh.neighbour();
     const vectorField& Cf = mesh.faceCentres();
     const vectorField& C = mesh.cellCentres();
     const vectorField& Sf = mesh.faceAreas();        
    
     forAll(this->patch(), facepi )        
      {
          label cellA = this->patch().faceCells()[facepi];
          const cell& faces = mesh.cells()[cellA];
          vector n = this->patch().Sf()[facepi]/mag(this->patch().Sf()[facepi]);
          vector fx = this->patch().Cf()[facepi];
        
          List<scalar> x(faces.size(), 0.); 
          List<Type> y(faces.size(), pTraits<Type>::zero ); 
    
          // Interpolate values to faces and compute normal distances       
          int id(0);
          scalar xav(0.);
          Type yav(pTraits<Type>::zero);
          forAll(faces, fi)     
           {
             label  facei = faces[fi];  
             
             // Only use internal cells (pitfall: faces on coupled patches
             // will not be used).
             if (mesh.isInternalFace(facei))
              {
                scalar SfdOwn = mag(Sf[facei] & (Cf[facei] - C[owner[facei]]));
                scalar SfdNei = mag(Sf[facei] & (C[neighbour[facei]] - Cf[facei]));
                scalar w = SfdOwn/(SfdOwn + SfdNei);
              
                y[id] = w*var[neighbour[facei]] + (1.-w)*var[owner[facei]];
                yav += y[id];
                
                x[id] = mag(n&(fx-Cf[facei]));
                xav += x[id]; 
                
                id++;            
              }                   
           }  
         
          // Last pair x-y is for cell P (there is always space for it in the lists
          // because at least one of the cell's faces is on the boundary)
        
          y[id] = var[cellA];
          yav += y[id]; 
          x[id] = mag(n&(fx-C[cellA]));
          xav += x[id]; 
          id++;
          
          yav /= id;
          xav /= id;
          
          // Compute regression
          Type num(pTraits<Type>::zero);
          scalar den(0.);
          for (int i=0; i<id; i++)
           {
             num += (x[i]-xav)*(y[i]-yav);
             den += (x[i]-xav)*(x[i]-xav); 
           }
          
          varp[facepi] = yav - xav*num/den;    
      }  
    
   }
   
   this->operator==(varp);      
   fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void linearExtrapolationFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeKeyword("useRegression") << useReg_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
