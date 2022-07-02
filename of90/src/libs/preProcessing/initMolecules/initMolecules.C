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

#include "fvCFD.H"
#include "IOPositionMolec.H"

#include "Random.H"
#include <string> 
#include "IFstream.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class userInput
{
 public:
 
   List<int> startBead;
   List<int> parentBranchStart;
   List<int> nbeads;
   
   List<int> endBead;
   List<int> parentBranchEnd;
        
   List<vector> growV;
        
   List<bool> isRand;
   List<bool> isOpen;
 
};

int createMoleculesGI
(
  const int gi,
  const userInput& uI,
  const Field<vector>& p0,
  const scalar Ls,
  const int nMinc,
  List<vector>& xL,
  Field<Field<label> >& idsIn,
  List<List<label> >& sprIn
);
 
int main(int argc, char *argv[])
{
   argList::noParallel();
    
   #include "setRootCase.H"
   #include "createTime.H"
   #include "createMesh.H"
   
   // Dict with the info for init 
   IOdictionary  molcDict_
   (
     IOobject
     (
       "initMoleculesDict",
       runTime.constant(),
       runTime,
       IOobject::MUST_READ_IF_MODIFIED,
       IOobject::NO_WRITE,
       false
     )
   );   
       
   PtrList<entry> allG(molcDict_.lookup("groups"));
   
   // Find the number of molecules  
   int nAllMolc(0);  int nM(0); 
   int nAllbeads(0);  
   forAll(allG, i)
    {
      allG[i].dict().lookup("nMolecules") >> nM;
      nAllMolc += nM;
    }
   if (nAllMolc<1)
    {
      FatalErrorIn("initMolecules")
      << "\n0 molecules defined.\n"
      << "Please check constant/initMoleculesDict and ensure that at least one valid molecule is defined."
      << exit(FatalError);
   }

   // Create positions and ids
   List<vector> xL;
   Field<Field<label> > idsTmp;
   List<List<label> > sprTmp;
   vector p1(0,0,0), p2(0,0,0);
   scalar Lsi(0.);
   labelField nbeadsgi(allG.size(), 0); // No. of beads for molecules in each group
    
   int nMinc(0);  // Increments the numb of molecules after each group iterate (helper to build springs)
   forAll(allG, i)
    {    
        int nM(0);
        allG[i].dict().lookup("nMolecules") >> nM;
        allG[i].dict().subDict("spatialDistibutionCoeffs").lookup("p0") >> p1;
        allG[i].dict().subDict("spatialDistibutionCoeffs").lookup("p1") >> p2;
        allG[i].dict().lookup("Ls") >> Lsi;
        
        // 1- Create the line over which the molecules' first bead will lie 
        vectorField p0(nM, vector::zero);
        if (mag(p1-p2)>1e-20)
         {       
          p0[0] = p1;
        
          if (nM>1)
           {
            vector dX((p2-p1)/(nM-1.));
      
            for (int j=1; j<nM; j++)
             {
               p0[j] = p0[j-1] + dX;
             }
           }         
         }
        else
         {
          p0 = p1;
         }
        
        
        // 2- Read and create each branch. Input has fixed format.
        
        ITstream& is(allG[i].dict().subDict("spatialDistibutionCoeffs").lookup("branches"));
         
        int ti;
        vector tv;
        bool tb;
        
        userInput uI;
        
        token lastToken(is); //(
        while(is.nRemainingTokens() > 1)
         {            
           is >> lastToken;  //( 
           is >> ti; uI.parentBranchStart.append(ti);         
           is >> ti; uI.startBead.append(ti);            
           is >> ti; uI.nbeads.append(ti);  
           is >> lastToken;  //)
           
           is >> lastToken;  //(
           is >> tv.x();  
           is >> tv.y();  
           is >> tv.z(); uI.growV.append(tv);   
           is >> lastToken;  //)
           
           is >> lastToken;  //(
           is >> tb; uI.isRand.append(tb);  
           is >> tb; uI.isOpen.append(tb);  
           is >> lastToken;  //) 
           if (!uI.isOpen[uI.isOpen.size()-1])
            {
               is >> lastToken;  //(            
               is >> ti; uI.endBead.append(ti);  
               is >> ti; uI.parentBranchEnd.append(ti);  
               is >> lastToken;  //)
            }  
           else
            {
             // Fill with arbitrary values just to avoid uninit errors in case
             // of not all closed
             uI.endBead.append(0); 
             uI.parentBranchEnd.append(0);
            }        
         }  
         
        
         nbeadsgi[i]
          =
         createMoleculesGI
           (
             i,
             uI,
             p0,
             Lsi,
             nMinc,
             xL,
             idsTmp,
             sprTmp           
           );
           
         nAllbeads += nbeadsgi[i] * nM; 
         nMinc += nM;
    }

    IOField<Field<label> > ids
    (
      IOobject
       (
        "indices",   
        runTime.timePath(),
        "lagrangian/molecules",             
        runTime,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
       ),
      idsTmp
    );
    idsTmp.clear();
    
    IOList<List<label> > spr
    (
      IOobject
       (
         "springs",
         runTime.constant(),
         "runTimeInfo"/runTime.timeName(),
         runTime,
         IOobject::NO_READ,
         IOobject::NO_WRITE,
         false
       ),
      sprTmp
    );
    sprTmp.clear();

   //- Create positions (special object; has its own writeData) 
   IOPositionMolec xLIO(xL, runTime, mesh); 
   
   //- origProcId (only needed for parallel  run; required as default by of)
   
   labelField origProcIdtmp(nAllbeads, 0); 
   IOField<label> origProcId
    (
      IOobject
       (
        "origProcId",   
        runTime.timePath(),
        "lagrangian/molecules",               
        runTime,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
       ),
      origProcIdtmp
    );
   origProcIdtmp.clear();
   
   //- origId (on construction this an alias of ids[:][0]; required as default by of)
   
   labelField origIdtmp(nAllbeads, 0); 
   forAll(origIdtmp, bi)
   {
     origIdtmp[bi] = ids[bi][0];
   }
   
   IOField<label> origId
    (
      IOobject
       (
        "origId",   
        runTime.timePath(),
        "lagrangian/molecules",               
        runTime,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
       ),
      origIdtmp
    );
   origIdtmp.clear();
   
   //- molcID
    
   IOField<label > molcID_
   (
     IOobject
      (
       "molcID",
       runTime.timePath(),
       "lagrangian/molecules",               
       runTime,
       IOobject::NO_READ,
       IOobject::NO_WRITE,
       false
      ),
     nAllbeads
   );
        
   {
   label ri(0);
   label mi(0);
   forAll(allG, i)
    {
      int nM(0);
      allG[i].dict().lookup("nMolecules") >> nM;
      int nB = nbeadsgi[i];
      
      for (int j=0; j<nM; j++)
       {
         for (int k=0; k<nB; k++)
          {
            // MolcID
            molcID_[ri] = mi;
            ri++;
          }
          
         mi++; 
       }        
    } 
   } 
    
   //- MoleculesInfo (dict, not a field)
       
   IOdictionary MIDict
    (
        IOobject
        (
            "MoleculesInfo",
            runTime.constant(),
            "runTimeInfo"/runTime.timeName(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
   );
        
   List<List<label> > actMolc_(nAllMolc, List<label>(3,0));
      
   //-- Group properties
   List<word> names(allG.size(), "");
   List<scalar> Nks(allG.size(), 0.);
   List<scalar> nuEV(allG.size(), 0.);
   List<scalar> D(allG.size(), 0.);
   List<scalar> a(allG.size(), 0.);
   List<scalar> Ls(allG.size(), 0.);
       
   {
   label ri(0);
   label gi(0);
   label jprev(0);
   List<label> tmp(3,0);
   forAll(allG, i)
    {
      int nM(0);
      allG[i].dict().lookup("nMolecules") >> nM;
      int nB = nbeadsgi[i];
      
      names[i] = allG[i].keyword();
      allG[i].dict().lookup("Nks") >> Nks[i];
      allG[i].dict().lookup("nuEV") >> nuEV[i];
      allG[i].dict().lookup("D") >> D[i];
      allG[i].dict().lookup("a") >> a[i];
      allG[i].dict().lookup("Ls") >> Ls[i];
      
      for (int j=0; j<nM; j++)
       { 
         // Molc name | Nb of beads | Group
         tmp[0] = jprev+j;
         tmp[1] = nB; 
         tmp[2] = gi; 
         
         actMolc_[ri] = tmp;
         ri++;          
       } 
       
      jprev += nM;  
      gi++;     
    } 
   } 
 
   word actM("ActiveMolecules");
   word gP("GroupProperties");
   MIDict.add(actM, dictionary());
   MIDict.add(gP, dictionary());
   MIDict.subDict(actM).add("nActive", nAllMolc);
   MIDict.subDict(actM).add("idList", actMolc_);
   MIDict.subDict(gP).add("names", names);
   MIDict.subDict(gP).add("Nks", Nks);
   MIDict.subDict(gP).add("nuEV", nuEV);
   MIDict.subDict(gP).add("D", D);
   MIDict.subDict(gP).add("a", a);
   MIDict.subDict(gP).add("Ls", Ls);
        
   // Write all fields  
   xLIO.write();
   ids.write(nAllMolc > 0);
   origProcId.write(nAllMolc > 0);
   origId.write(nAllMolc > 0);
   molcID_.write(nAllMolc > 0);
   
   // Forced removal of "runTimeInfo" to avoid accumulation of garbage
   if (exists(runTime.constant()/"runTimeInfo"))
      rmDir(runTime.constant()/"runTimeInfo");
   
   MIDict.writeObject( IOstream::ASCII, IOstream::currentVersion, runTime.writeCompression(), nAllMolc > 0); 
   spr.write();
   
   // Output statistics
   Info << "Molecules created with success: "<< nl
        << "Molecules Stats " << nl
        << tab << "Total groups" << tab << allG.size() << nl 
        << tab << "Total molecules" << tab << nAllMolc << nl 
        << tab << "Total beads" << tab << nAllbeads << nl 
        << tab << "Total springs" << tab << spr.size() << nl 
        << endl;
        
   Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << tab 
       << "ClockTime = " << runTime.elapsedClockTime() << " s"
       << nl << endl;
       
   return 0; 
}

int createMoleculesGI
(
  const int gi,
  const userInput& uI,
  const Field<vector>& p0,
  const scalar Ls,
  const int nMinc,
  List<vector>& xL,
  Field<Field<label> >& idsIN,
  List<List<label> >& sprIN
)
{
 
  Random rndSeed(std::time(NULL)); 
  
  int nBranches = uI.startBead.size();
  
  // tnb - number of beads per molecule of group i
  // (constant for all the molecules in same group)
  // nSpr - number of springs in the molecules. The main backbone
  // contributes with (nbeads - 1) springs and each branch has
  // (nbeads) springs. If the branch is closed, the number of 
  // springs increases by 1 in each case. 
  int tnb(0);
  int nSpr(0);
  forAll(uI.nbeads, Bi)
   {
      tnb += uI.nbeads[Bi];
      
      if (uI.isOpen[Bi])
       {
         if (Bi == 0)
          {
            nSpr += uI.nbeads[Bi]-1; 
          } 
         else
          {
            nSpr += uI.nbeads[Bi]; 
          }
       }
      else
       {
         if (Bi == 0)
          {
            nSpr += uI.nbeads[Bi]; 
          } 
         else
          {
            nSpr += uI.nbeads[Bi]+1; 
          }
       }
   } 
  
  // loc. ID bead i-1 | loc. ID bead i | molcID
  List<List<label> > spr(nSpr*p0.size(), List<label>(3, 0) );
 
  Field<Field<label> > ind(tnb*p0.size(), Field<label>(3, 0) );
   
  int rilast(idsIN.size()); // We start from where we finish for previous group
 
  Field<label> ids(3,0);
  
  int si(0), ri(0);
  forAll (p0, mi)
  { 
   
   List<List<vector> > txL(nBranches, List<vector>() );
  
   int bli(0);
   forAll(uI.startBead, Bi)
    {
      int pB = uI.parentBranchStart[Bi];
      int sb = uI.startBead[Bi]; 
      int nb = uI.nbeads[Bi];
      bool isOp = uI.isOpen[Bi];
      bool isRnd = uI.isRand[Bi];  
      vector gV = uI.growV[Bi];
     
      // First bead
      //- Position and Springs
      txL[Bi].setSize(nb);
      if (Bi == 0)
       {
          txL[Bi][0] = p0[mi];
       }
      else
       {
          // Position (we force no-random in the parent-child connection)
          txL[Bi][0] = txL[pB][sb] + gV*Ls;
          
          // Springs
          // We have the chain ID of the parent chain/bead for the connection
          // but we want the local ID of that bead. Since we follow the 
          // order written by the user, it is just counting the beads in
          // branches until arriving to the parent branch. 
          int locID(-1);
          for (int i=0; i<pB; i++)
           {
             locID += uI.nbeads[i]; 
           }
          locID += sb+1;
           
          spr[si][0] = locID;
          spr[si][1] = bli; 
          spr[si][2] = mi+nMinc;
          si++;
       }
       
      //- Indices
      //-- Bead name | Local ID | GroupID
      ids[0] = ri+rilast; ids[1] = bli; ids[2] = gi;  
      ind[ri] = ids;
      ri++;
      bli++;
       
      // Remaining beads
      for (int bi=1; bi<nb; bi++)
       {      
          //- Position
          if (!isRnd)
           {
             txL[Bi][bi] = txL[Bi][bi-1] + gV*Ls;
           }
          else
           {
             vector gVRnd
             (
               ( -1. + 2.*rndSeed.scalar01() )*gV.x(),
               ( -1. + 2.*rndSeed.scalar01() )*gV.y(),
               ( -1. + 2.*rndSeed.scalar01() )*gV.z() 
             );
              
             txL[Bi][bi] = txL[Bi][bi-1] + gVRnd*Ls;
           }
            
          //- Springs
          spr[si][0] = bli-1; 
          spr[si][1] = bli; 
          spr[si][2] = mi+nMinc;  
          si++;
              
          //- Indices
          //-- Bead name | Local ID | GroupID
          ids[0] = ri+rilast; ids[1] = bli; ids[2] = gi; 
          ind[ri] = ids;
          ri++;
          bli++;           
       }
         
      // Build the ending spring for closed chains. Note that this operation
      // only adds one more springs, without changing the number of beads
      if (!isOp)
       {
          int pBe = uI.parentBranchEnd[Bi];
          int sbe = uI.endBead[Bi]; 
            
          int locID(-1);
          for (int i=0; i<pBe; i++)
           {
             locID += uI.nbeads[i]; 
           }
          locID += sbe+1;
           
          //- Springs
          spr[si][0] = locID; 
          spr[si][1] = bli-1; // Subtract 1 because increment has been done inside previous for loop 
          spr[si][2] = mi+nMinc;  
          si++;      
       }
      
    }
   
    // Transfer position of this molecule to the 
    // global container
    forAll(txL, Bi)
     {
        xL.append(txL[Bi]);
     }    
  }
  
  // Transfer indexes and springs to the main container
  idsIN.append(ind);
  sprIN.append(spr);
   
  return tnb;
  
}
// ************************************************************************* //
