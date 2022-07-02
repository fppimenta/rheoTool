/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::spinesFvMotionSolver

Description
    Mesh motion solver for an fvMesh. Based on spines. Case-specific, 
    cannot be generalized. Points only move in y-dir, both x and z keep
    fixed.

SourceFiles
    spinesFvMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef spinesFvMotionSolver_H
#define spinesFvMotionSolver_H

#include "displacementMotionSolver.H"
#include "fvMotionSolver.H"

#include <vector>
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{ 

/*---------------------------------------------------------------------------*\
             Class spinesFvMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class spinesFvMotionSolver
:
    public displacementMotionSolver,
    public fvMotionSolver
{
    // Private Data
 
        //- Optionally read point-position field. Used only for position
        //  boundary conditions.
        mutable autoPtr<pointVectorField> pointLocation_;
 
        //- Frozen points (that are not on patches). -1 or points that are
        //  fixed to be at points0_ location
        label frozenPointsZone_;
        
        //- Retrives the current mesh points based on spines method
        tmp<pointField> curPointsSequential() const;
        
        //- Same as previous, but for parallel runs
        tmp<pointField> curPointsParallel() const;
        
        //- Give to each point of the mesh (both internal and boundary) an
        // index indentifying if it is interior or if it belongs to a given
        // boundary. This is just to split between points that can be moved
        // freely (according to spines) and the ones that have constrained motion.
        // Note that some boundary points also move with spines, as e.g. those at the outlet.
        void markPatches();
        
        //- Number of mesh points along each Cartesian direction
        label npx_;
        label npy_;
        label npz_;
        
        //- Contains the ID of each mesh point. See markPatches().
        mutable List<label> pID_;
        
        //- Just the global adressing, used in parallel runs.
        mutable autoPtr<labelIOList> globalAdr_;
        
        //- Local index of the mesh points belonging to patch freeSurface0,
        // order by increasing x-coord. Uses global addressing for identification
        // and the resulting order is the same as if the mesh was not decomposed 
        // (structured-like algorithm).
        // (Only used for parallel runs) 
        List<label> fsOrder0_;
        
        //- Same as previous, for freeSurface1, when it exists.
        // (Only used for parallel runs) 
        List<label> fsOrder1_;
        
        // Only one free-surface? 
        // (Only used for parallel runs) 
        bool oneFS_;
 
public:

    //- Runtime type information
    TypeName("spines");


    // Constructors

        //- Construct from polyMesh and dictionary
        spinesFvMotionSolver
        (
            const polyMesh&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        spinesFvMotionSolver
        (
            const spinesFvMotionSolver&
        );


    //- Destructor
    ~spinesFvMotionSolver();


    // Member Functions

       
         
        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Update topology
        virtual void updateMesh(const mapPolyMesh&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const spinesFvMotionSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
