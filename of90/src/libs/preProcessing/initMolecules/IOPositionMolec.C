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

#include "volFields.H"
#include "IOPositionMolec.H"
#include "particle.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::IOPositionMolec::IOPositionMolec(const List<vector>& p, const Time& rt, const polyMesh& mesh)
:        
    regIOobject
    (
        IOobject
        (
           "positions",
            rt.timePath(),
           "lagrangian/molecules",        
            rt,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
    positions_(p),
    typeN_("positionMolecule"),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::IOPositionMolec::write(const bool wr) const
{
    if (positions_.size())
    {
        return regIOobject::write(wr);
    }
    else
    {
        return true;
    }
}

// We create a single particle and we track it to the XYZ positions
// of p. Write the barycentric coordinates after each movement. This 
// is a fast method, but accumulates some error.

bool Foam::IOPositionMolec::writeData(Ostream& os) const
{    
    os  << positions_.size() << nl << token::BEGIN_LIST << nl;
    
    particle* pt;
    pt = new particle(mesh_, positions_[0], -1);
    pt->writePosition(os);
    scalar f = 0;
    
    for (int i = 1; i<positions_.size(); i++)
     {      
       if ( pt->track(positions_[i]-positions_[i-1], f) > 0.01)
        {
          FatalErrorIn("Foam::IOPositionMolec::writeData(Ostream& os)")
          << "\nBead with position: " << positions_[i] <<" is probably outside the mesh." 
          << " Check constant/initMoleculesDict for p0 and p1 and/or modify the mesh."
          << exit(FatalError);
        }    
       pt->writePosition(os);
           
       os  << nl;       
     }
    
    delete pt; 
    pt = nullptr;
    
    os  << token::END_LIST << endl;
    
    return os.good();
}

bool Foam::IOPositionMolec::readData(Istream& is)
{
   // Should not be needed 
   return true; 
}


// ************************************************************************* //
