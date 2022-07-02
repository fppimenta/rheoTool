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

#include "solidParticle.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const std::size_t Foam::solidParticle::sizeofFields_
(
    sizeof(solidParticle) - sizeof(particle)
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidParticle::solidParticle
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    particle(mesh, is, readFields)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            is >> ids_;
            is >> molcID_;
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&ids_), 
                sizeof(ids_) + sizeof(molcID_) 
            );
        }
    }
    
    // Check state of Istream
    is.check("solidParticle::solidParticle(Istream&)");
  
}


void Foam::solidParticle::readFields(Cloud<solidParticle>& c)
{
 
    if (!c.size())
    {
        return;
    }

    particle::readFields(c);
  
    IOField<Field<label> > ids(c.fieldIOobject("indices", IOobject::MUST_READ));
    c.checkFieldIOobject(c, ids); 
    
    IOField<label> molcID(c.fieldIOobject("molcID", IOobject::MUST_READ));
    c.checkFieldIOobject(c, molcID);  
    
    // Initialize U as (0,0,0), without reading it.
    // U is NO_READ and NO_WRITE.
    vectorField U(ids.size(), vector::zero);  

    label i = 0;
    forAllIter(Cloud<solidParticle>, c, iter)
    {
        solidParticle& p = iter();

        p.ids_ = ids[i];
        p.molcID_ = molcID[i];
        p.U_ = U[i];
        
        i++;
    }
 
}


void Foam::solidParticle::writeFields(const Cloud<solidParticle>& c)
{
 
    particle::writeFields(c);

    label np = c.size();

    IOField<Field<label> > ids(c.fieldIOobject("indices", IOobject::NO_READ), np);
    IOField<label> molcID(c.fieldIOobject("molcID", IOobject::NO_READ), np);
    
    label i = 0;
    forAllConstIter(Cloud<solidParticle>, c, iter)
    {
        const solidParticle& p = iter();

        ids[i] = p.ids_;
        molcID[i] = p.molcID_;
       
        i++;
    }

    ids.write();
    molcID.write(); 
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const solidParticle& p)
{
 
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const particle&>(p)
            << token::SPACE << p.ids_
            << token::SPACE << p.molcID_;
    }
    else
    {
        os  << static_cast<const particle&>(p);
        os.write
        (
              reinterpret_cast<const char*>(&p.ids_),
              sizeof(p.ids_) + sizeof(p.molcID_)
        );
    }

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const solidParticle&)");

    return os;
}


// ************************************************************************* //
