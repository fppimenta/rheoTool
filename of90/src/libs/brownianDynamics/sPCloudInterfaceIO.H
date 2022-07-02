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

#include "sPCloudInterface.H"
#include "particle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Functions  * * * * * * * * * * * * * * * //

bool sPCloudInterface::checkIfExist()
{      
    IOobject IONamesList
    (
      "molcID",
      U().time().timeName(),
      "lagrangian/molecules",      
      U().db(),
      IOobject::MUST_READ
    );
   
    if (!IONamesList.typeHeaderOk<labelField>(false))
    {         
        FatalErrorIn("sPCloudInterface::checkIfExist()")
        << "\nFile " << U().time().timeName() << "/lagrangian/molecules" << " could not be found.\n\n"
        << "Run initMolecules to create a valid set of molecules.\n"
        << exit(FatalError);
    }        
    
    return true;
}

void sPCloudInterface::writeRunTimeInfoDict()
{    
  // We keep only the active molecules inside MoleculesInfo,
  // but we keep the info on all groups for the properties.

    IOdictionary uD
    (
        IOobject
        (
            "MoleculesInfo",
            U().time().constant(),
            "runTimeInfo"/U().time().timeName(),
            U().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );
    
    List<List<label> > actMolc_(nMolc_, List<label>(3,0));
    
    // We have to update the active molecules
    int j(0); 
    forAll(linkM_, i)
     {
       if (linkM_[i] != -1)
        {
          label mi(linkM_[i]);         
          actMolc_[j] = mAct_[mi]; 
          j++;
        }
     }

    word actM("ActiveMolecules");
    word gP("GroupProperties");
    uD.add(actM, dictionary());
    uD.add(gP, dictionary());
    uD.subDict(actM).add("nActive", nMolc_);
    uD.subDict(actM).add("idList", actMolc_);
    uD.subDict(gP).add("names", names_);
    uD.subDict(gP).add("Nks", Nks_);
    uD.subDict(gP).add("nuEV", nuEV_);
    uD.subDict(gP).add("D", D_);
    uD.subDict(gP).add("a", a_);
    uD.subDict(gP).add("Ls", Ls_);
        
    uD.writeObject( IOstream::ASCII, IOstream::currentVersion, U().time().writeCompression(), j>0 );
 
    // Write the springs (connections)
    List<List<label> > sprTmp;
    forAll(mSpr_, mi)
     {
       if (mAct_[mi][0] != -1)
        {
          sprTmp.append(mSpr_[mi]);
        }
     }
     
    IOList<List<label> > spr
    (
      IOobject
       (
         "springs",
         U().time().constant(),
         "runTimeInfo"/U().time().timeName(),
         U().mesh(),
         IOobject::NO_READ,
         IOobject::NO_WRITE,
         false
       ),
      sprTmp
    );  
  
    spr.write();
}

void sPCloudInterface::writeM()
{ 
   if (U().time().outputTime())
    {    
      writeRunTimeInfoDict();
        
      if (writeVTK_)
        writeVTK();
      
    // if writeContFields_ is false the main solver is not writing anything, thus do it here
    // for lagrangian fields and for the time    
      if (!writeContFields_)
      {
        spc_.write();
        
        // Write uniform/time
        U().time().writeTimeDict();       
      }      
    }    
}

void sPCloudInterface::writeStatistics()
{
 ppCount_++;
 
 if (ppCount_>outputStatsN_)
 { 
  
  //- Time (write to all files except id file, the last one)
  forAll(outS_, gi)
   {
     for (int i = 0; i<outS_[gi].size()-1; i++)
      {
        outS_[gi][i]() << U().mesh().time().value() << tab;
      }
   }
      
  forAll(mx_, mi)
   {
      label gi = mAct_[mi][2];
      
      if (mAct_[mi][0] != -1)
       {                    
         // Center of mass
         vector cM(0.,0.,0.);
         
         forAll(mx_[mi], bi)
         {
           cM += mx_[mi][bi]; 
         }
         cM /= mx_[mi].size(); 

         //  Max inter-bead distance (apparent molecule length)
         scalar maxd(0.);
         for (int bi=0; bi<mx_[mi].size()-1; bi++)
         {
           for (int bj=bi+1; bj<mx_[mi].size(); bj++)
            {
              scalar di (mag(mx_[mi][bi] - mx_[mi][bj]));
              if ( di > maxd)
               {
                 maxd = di; 
               }
            }
         }
         
         //  Write
         outS_[gi][0]() << cM.x() << tab << cM.y() << tab << cM.z() << tab;
         outS_[gi][1]() << maxd << tab;     
       }
      else
       {
         // Lost molecules get 0. To know if a given 0 is real or means lost,
         // look to the stretch value: 0 will always mean lost. A
         // token::SPACE would be annoying while reading.
         outS_[gi][0]() << 0 << tab << 0 << tab << 0 << tab;
         outS_[gi][1]() << 0 << tab; 
       }
   }
   
  forAll(outS_, gi)
   {
     for (int i = 0; i<outS_[gi].size()-1; i++)
      {
        outS_[gi][i]() << endl;
      }
   }
   
  // Reset counter
  ppCount_ = 1;
 }
 
}

void sPCloudInterface::writeVTK()
{  
   // Index in the name of the file
   fileName   VTFDir(VTKDir_/word("molecules_00" + Foam::name(U().mesh().time().timeIndex()) + ".vtk"));
   OFstream   VTF(VTFDir);
   
   //- Header
   VTF << "# vtk DataFile Version 2.0" << nl;
   VTF << "Molecules data" << nl;
   VTF << "ASCII" << nl;
   VTF << "DATASET UNSTRUCTURED_GRID" << nl;
 
   //- Points coordinates
   int np(0);
   forAll(mU_, mi)
   {
      if (mAct_[mi][0] != -1)
       {
         np += mU_[mi].size();
       }
   }   
   VTF << "POINTS " << np << " double" << nl;
   
   // miTmp contains the global VTK ID of the points, for each
   // spring, in its first label of the list (remaining are ignored). We
   // cannot use the global ID of the beads because it will not match the
   // global VTK ID once the total number of beads changes (this is when
   // the global VTK ID and the global bead ID start to unmatch)
     
   PtrList<List<List<label > > > miTmp = mIds_;
   label cnt(0); 
   label nSpr(0);
   forAll(mx_, mi)
   {
      if (mAct_[mi][0] != -1)
       {
         forAll(mx_[mi], bi)
          {
            VTF << mx_[mi][bi].x() << token::SPACE << mx_[mi][bi].y() << token::SPACE << mx_[mi][bi].z() << token::SPACE;
            miTmp[mi][bi][0] = cnt; 
            cnt++;
          }
         nSpr += mSpr_[mi].size();
       }
   }   
   VTF << nl;
  
   //- Cells (springs). Molecules would be more economic in writting,
   // but it would not allow to differentiate between branches one 
   // using polyline cell type. No problem regarding the data, because
   // it is associated with points (beads), not cells
   VTF << "CELLS " << nSpr << token::SPACE << nSpr*3 << nl;  
   
   // Loop over each spring
   forAll(mSpr_, mi)
   {
      if (mAct_[mi][0] != -1)
       {
         forAll(mSpr_[mi], bi)
          {
            VTF << 2 << token::SPACE <<  miTmp[mi][mSpr_[mi][bi][0]][0] << token::SPACE << miTmp[mi][mSpr_[mi][bi][1]][0] << nl;
          }
       }
   }  
    
   VTF << "CELL_TYPES " << nSpr << nl; 
   forAll(mSpr_, mi)
   {
      if (mAct_[mi][0] != -1)
       {
         forAll(mSpr_[mi], si)
          {
            VTF << 3 << nl; // type of cell = line
          }
       }
   } 
     
   //- Point data  
   VTF << "POINT_DATA " << np << nl; 
   VTF << "SCALARS globalID int 1" << nl; 
   VTF << "LOOKUP_TABLE default" << nl; 
   forAll(mx_, mi)
   {
      if (mAct_[mi][0] != -1)
       {
         forAll(mx_[mi], bi)
          {
            VTF << mIds_[mi][bi][0] << token::SPACE;
          }
       }
   } 
   VTF << nl;
   VTF << "SCALARS localID int 1" << nl;
   VTF << "LOOKUP_TABLE default" << nl; 
   forAll(mx_, mi)
   {
      if (mAct_[mi][0] != -1)
       {
         forAll(mx_[mi], bi)
          {
            VTF << mIds_[mi][bi][1] << token::SPACE;
          }
       }
   } 
   VTF << nl;
   VTF << "SCALARS groupID int 1" << nl;
   VTF << "LOOKUP_TABLE default" << nl; 
   forAll(mx_, mi)
   {
      if (mAct_[mi][0] != -1)
       {
         forAll(mx_[mi], bi)
          {
            VTF << mIds_[mi][bi][2] << token::SPACE;
          }
       }
   }
   VTF << nl;   
   VTF << "SCALARS molcID int 1" << nl;
   VTF << "LOOKUP_TABLE default" << nl; 
   forAll(mx_, mi)
   {
      if (mAct_[mi][0] != -1)
       {
         forAll(mx_[mi], bi)
          {
            VTF << mAct_[mi][0] << token::SPACE;
          }
       }
   }

   VTF << endl;
   
   //- Cell data 
       
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
