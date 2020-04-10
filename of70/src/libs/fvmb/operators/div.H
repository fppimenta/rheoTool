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

InNamespace
    Foam::fvmb

Description
    Implicit operators arising in coupled systems. Only a few are defined,
    which are not templated (fixed input/output). 

SourceFiles
    blockOperators.C

\*---------------------------------------------------------------------------*/

#ifndef div_H
#define div_H

#include "volFields.H"
#include "surfaceFields.H"
#include "tmp.H"
#include "dimensionedTypes.H"
#include "zero.H"

#include "LMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fvmb functions Declaration
\*---------------------------------------------------------------------------*/

namespace fvmb
{
   
   
   // Implicit operator: div(tau)
   // Input: c - multiplicative constant
   //        tau - volSymmTensorField
   // Output: LMatrix<tensor>
   tmp< LMatrix<tensor> >  div
   (
     const dimensionedScalar c,
     const GeometricField<symmTensor, fvPatchField, volMesh>& vf
   );
   
   tmp< LMatrix<tensor> >  div
   (
     const scalar c,
     const GeometricField<symmTensor, fvPatchField, volMesh>& vf
   );
   
   // Implicit operator: div(U)
   // Input: c - multiplicative constant
   //        U - volVectorField
   // Output: LMatrix<vector>
   tmp< LMatrix<vector> >  div
   (
     const dimensionedScalar c,
     const GeometricField<vector, fvPatchField, volMesh>& vf
   );
   
   tmp< LMatrix<vector> >  div
   (
     const scalar c,
     const GeometricField<vector, fvPatchField, volMesh>& vf
   );
   
   

} // End namespace fvmb

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
