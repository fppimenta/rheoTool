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

Class
    Foam::IOPositionMolec

Description
    Helper IO class to read and write particle positions

SourceFiles
    IOPositionMolec.C

\*---------------------------------------------------------------------------*/

#ifndef IOPositionMolec_H
#define IOPositionMolec_H

#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class IOPositionMolec Declaration
\*---------------------------------------------------------------------------*/


class IOPositionMolec
:
    public regIOobject
{

    // Private data

        //- Reference to the positions
        const List<vector>& positions_;
        
        //- Dummy
        word typeN_; 
        
        const polyMesh& mesh_;


public:

    // Static data
    
        //- Runtime type name information.
        // Dummy.
        virtual const word& type() const
        {
              return typeN_;
        }


    // Constructors

        //- Construct from positions (vector) field 
        IOPositionMolec(const  List<vector>& p, const Time& rt, const polyMesh& mesh);


    // Member functions

        virtual bool readData(Istream&);

        virtual bool write(const bool write = true) const;

        virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "IOPositionMolec.C" 
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
