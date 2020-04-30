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
 
namespace Foam
{
namespace fvmb
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
LMatrix<Type>::LMatrix
(
  const fvMesh& mesh
)
:
tmp<LMatrix<Type>>::refCount(),
LduMatrix<Type,Type,Type>::LduMatrix(mesh),
internalCoeffs_(mesh.boundary().size()),
boundaryCoeffs_(mesh.boundary().size()),
mesh_(mesh)
{
  // Initialize BC contribs
  forAll(mesh.boundary(), patchi)
  {
    internalCoeffs_.set
    (
      patchi,
      new Field<Type>
      (
        mesh.boundary()[patchi].size(),
        Zero
      )
    );

    boundaryCoeffs_.set
    (
      patchi,
      new Field<Type>
      (
        mesh.boundary()[patchi].size(),
        Zero
      )
    );
  }

}
 

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
  
 
} // End namespace fvmb
} // End namespace Foam
// ************************************************************************* //
