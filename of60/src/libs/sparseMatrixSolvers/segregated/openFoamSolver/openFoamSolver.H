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

Class
    Foam::openFoamSolver

Description
    Simple wrapper to the default fvmatrix.solve() of OpenFOAM. 

\*---------------------------------------------------------------------------*/

#ifndef openFoamSolver_H
#define openFoamSolver_H

#include "sparseSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
 
/*---------------------------------------------------------------------------*\
                       Class openFoamSolver Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class openFoamSolver
:
    public sparseSolver<Type>
{
    // Private data
 

    // Private Member Functions

        //- Disallow default bitwise copy construct
        openFoamSolver(const openFoamSolver&);

        //- Disallow default bitwise assignment
        void operator=(const openFoamSolver&);


public:

    //- Runtime type information
    TypeName("openFoamSolver");


    // Constructors

        //- Construct from mesh
        openFoamSolver
        (
            const GeometricField<Type, fvPatchField, volMesh>& T,
            const fvMesh& mesh,
            const dictionary& fvSolution
        );
        
    // Destructor

        virtual ~openFoamSolver()
        {}
         

    // Member Functions
    
        // Does the sparse solver relies on any external package? 
        virtual bool isExternalSolver() const
        {
          return false;
        }
        
        // Interface to applications
        virtual void solve
        (
          fvMatrix<Type>&
        );
        
        virtual void solve
        (
          fvMatrix<Type>&,
          const dictionary&
        );
        
        virtual void solve
        (
          const tmp<fvMatrix<Type>>&,
          const dictionary&
        );

};

  
} // End namespace Foam

#ifdef NoRepository
    #include "openFoamSolver.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
