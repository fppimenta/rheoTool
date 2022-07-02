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
    
InNamespace
    Foam::fvmb
    
Class
    Foam::LMatrix
    
Description
    This is simply an LduMatrix to which we add internalCoeffs_ and 
    boundaryCoeffs_. It is sort of an fvMatrix, but because it derives
    from an LduMatrix, the matrix fields (diag, source, upper/lower) may
    be all field<Type>. This class could be even more general if we made use
    of the 3 different templates available in LduMatrix<Type1,Type2,Type3> (we
    are using a single template parameter). 
    Created to handle special operators arising in coupled systems.

\*---------------------------------------------------------------------------*/

#ifndef LMatrix_H
#define LMatrix_H

#include "tmp.H"
#include "LduMatrix.H"
#include "volFields.H"
 
namespace Foam
{
namespace fvmb
{

template<class Type>
class LMatrix
:
 public tmp<LMatrix<Type>>::refCount,
 public LduMatrix<Type,Type,Type> 
{
    // Private data
    
        // Contribution to diag
        FieldField<Field, Type> internalCoeffs_;
        
        // Contribution to source or upper/lower for coupled patches
        FieldField<Field, Type> boundaryCoeffs_;
        
        // Reference to mesh
        const fvMesh& mesh_;
        
        // Referencing list to be set by operators and used by coupledSolvers
        List<labelList>  rel_;
          

    // Private Member Functions

        //- Disallow default bitwise copy construct
        LMatrix(const LMatrix&);

        //- Disallow default bitwise assignment
        void operator=(const LMatrix&);
         
       
public:

    //- Runtime type information
    TypeName("LMatrix");

    // Constructors

        //- Construct from components
        LMatrix
        (
           const fvMesh& mesh
        );

        //- Hack to return tmp on demand 
        static tmp<LMatrix<Type>> New
        (
           const fvMesh& mesh
        )
        {
          return tmp<LMatrix<Type>> ( new LMatrix<Type>(mesh) ); 
        }

    // Destructor

        virtual ~LMatrix()
        {
        }


    // Member Functions (accessors)
        FieldField<Field, Type>& internalCoeffs()
        {
          return internalCoeffs_;
        }

        FieldField<Field, Type>& boundaryCoeffs()
        {
          return boundaryCoeffs_;
        }
        
        List<labelList>& rowColList() 
        {
          return rel_;
        }
        
        const fvMesh& mesh() const
        {
          return mesh_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvmb
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LMatrix.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
