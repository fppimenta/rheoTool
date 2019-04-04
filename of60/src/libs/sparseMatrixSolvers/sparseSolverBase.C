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

\*---------------------------------------------------------------------------*/

#include "IPstream.H"
#include "OPstream.H"
#include "processorFvPatch.H"
#include "processorCyclicFvPatch.H"
#include "cyclicFvPatch.H"

#include "sparseSolverBase.H"

// * * * * * * * * * * * * * * * * * Static Data * * * * * * * * * * * * * * * //

bool Foam::sparseSolverBase::debug_ = false;

int Foam::sparseSolverBase::counterHypre_ = 0;
int Foam::sparseSolverBase::counterPetsc_ = 0;

// Default initialize, set latter if needed
Foam::sparseSolverBase::sharedDataS Foam::sparseSolverBase::sharedData;

// * * * * * * * * * * * * * * * * * Static functions * * * * * * * * * * * * * * * //

void Foam::sparseSolverBase::setSharedData
( 
  const fvMesh& mesh
)
{

 // If serial run, just set lower/upper row indices. Remaining  
 // elements of struct are only for parallel runs
 if (!Pstream::parRun())
 {
   sharedData.ilower = 0;
   sharedData.iupper = mesh.nCells();
   return;
 }
 
 // Reset lists if topo changes
 if (mesh.topoChanging())
 {
    sharedData.fCo.clear();
    sharedData.fCn.clear();
    sharedData.procInfo.clear();
 }
 
 // Create mapping between global cell ID and processor ordered cell ID.
 // The local cell ID mapping to the above two is the same 
 
 // 1 label p/processor
 List< label > nCells(Pstream::nProcs());

 nCells[Pstream::myProcNo()] = mesh.nCells();
 
 Pstream::gatherList(nCells);
 
 if (Pstream::myProcNo() == 0)
 {
    int sum = 0;
    for (int j = 0; j<Pstream::nProcs(); j++)
    {
      int nTemp = nCells[j]; 
      nCells[j] += sum;
      sum += nTemp;  
    }
 }

 Pstream::scatter(nCells);
 
 labelList mapProcs(mesh.nCells(), 0);
 
 forAll(mapProcs, i)
 {
  // master needs no increment
  if (Pstream::myProcNo() == 0)
     mapProcs[i] = i;   
  else
     mapProcs[i] = i + nCells[Pstream::myProcNo()-1];   
 }
 
 if (debug_ && !Pstream::parRun())
 {
   labelIOList localCellProcAddr
   (
    IOobject
    (
	"cellProcAddressing",
	mesh.facesInstance(),
	mesh.meshSubDir,
	mesh,
	IOobject::MUST_READ,
	IOobject::NO_WRITE,
	false
    )
   );
 
   Pout << "Orig ID: " << localCellProcAddr << endl;
   Pout << "Mapped ID: " << mapProcs << endl;
 }
 
 sharedData.ilower = min(mapProcs);
 sharedData.iupper = max(mapProcs);
 
  
 //-- Transfer faceCells between processors stradling interfaces 

 
 // List with size equal to number of processors
 List< List< List<label > > > gatheredData(Pstream::nProcs());

 // Gather on master processor.
 List< List<label > >  procInter;
 forAll(mesh.boundary(), patchi)
 {
    if (isType<processorFvPatch>(mesh.boundary()[patchi]))
    {
      const processorPolyPatch& pp = refCast<const processorPolyPatch>( mesh.boundaryMesh()[patchi] );
      int procN = pp.neighbProcNo();
      int procO = pp.myProcNo(); 
        
      labelList tt{procO, procN, patchi, -1}; // processors get last index equal to -1
      procInter.append(tt);   
    }
    else if (isType<processorCyclicFvPatch>(mesh.boundary()[patchi])) 
    {
      const processorCyclicPolyPatch& pp = refCast<const processorCyclicPolyPatch>( mesh.boundaryMesh()[patchi] );
      int procN = pp.neighbProcNo();
      int procO = pp.myProcNo(); 
      label cycID = pp.referPatchID();
       
      labelList tt{procO, procN, patchi, cycID}; // processorCyclic get last index equal to the cyclic patch ID
      procInter.append(tt);   
     }
 }
 
 gatheredData[Pstream::myProcNo()] = procInter;
 Pstream::gatherList(gatheredData);

 // Distribute over all procs
 Pstream::scatterList(gatheredData);

 // Assemble into single list
 List< List<label > > procInterGlobal
 (
    ListListOps::combine< List< List<label > > >
    (
      gatheredData,
      accessOp< List< List<label > > >()
    )
 );

 // Transfer faceCells between processors pairs (already
 // using mapped IDs)
 for (int k = 0; k<procInterGlobal.size(); k++)
 {
   labelList lisTri = procInterGlobal[k];
   int o = lisTri[0]; 
   int n = lisTri[1];
   int p = lisTri[2]; 
   int t = lisTri[3]; 
  
   if (Pstream::myProcNo() == o)
   {  
     // Send to neig proc the name of the cyclic referPatch that the neig contains 
     if (t != -1)
     {    
         const fvPatch& cyclicPatch = mesh.boundary()[t];
    
         const fvPatch& nbrPatch = refCast<const cyclicFvPatch>
         (
             cyclicPatch
         ).neighbFvPatch();

         OPstream toNeig(Pstream::commsTypes::blocking, n);
         toNeig << nbrPatch.name(); 
     }
              
     IPstream fromNeig(Pstream::commsTypes::blocking, n);
     labelList tt;
     fromNeig >> tt;
     
     labelList fC = mesh.boundaryMesh()[p].faceCells();
     // Convert from local to global adressing
     for (int j = 0; j<fC.size(); j++)
       fC[j] = mapProcs[fC[j]];
     
     sharedData.fCo.append(fC);
     sharedData.fCn.append(tt);
     sharedData.procInfo.append(lisTri);        
   }
   else if (Pstream::myProcNo() == n)
   {
     word procNeigName;
     label pI;
     
     // We need to send to own the faceCells on the other side of the proc* patch,
     // ie in the neighbProcNo() processor. We could use built-in interfaces. Here
     // we follow a different approach, using the name patterns of proc* patches to    
     // get their ID in the neighbProcNo() processor. Note that each proc* patch
     // is an interface and has one instance in each side (processor) of that 
     // virtual interface, with different names (pattern is the same, only the
     // the indices change).  
     
     // processor
     if (t == -1)
     {
       procNeigName = "procBoundary" + Foam::name(n) + "to" + Foam::name(o);
       pI = mesh.boundaryMesh().findPatchID(procNeigName);
     }
     // processorCyclic
     else
     {
       // Receive from own proc the name of the cyclic referPatch that the neig contains 
       IPstream fromOwn(Pstream::commsTypes::blocking, o);
       word procCycNeigName;
       fromOwn >> procCycNeigName;
     
       procNeigName = "procBoundary" + Foam::name(n) + "to" + Foam::name(o) + "through" + procCycNeigName;
       pI = mesh.boundaryMesh().findPatchID(procNeigName);
     }     
     
     labelList fC = mesh.boundaryMesh()[pI].faceCells();
     
     // Convert from local to global adressing
     for (int j = 0; j<fC.size(); j++)
       fC[j] = mapProcs[fC[j]];
    
     OPstream toOwn(Pstream::commsTypes::blocking, o);
     toOwn << fC; 
   } 
   
   if (debug_)
   {
     if (Pstream::myProcNo() == o)
     {            
        Pout <<  "Faces on my proc: " <<  sharedData.fCo << endl;
        Pout <<  "Faces on neig proc: " <<  sharedData.fCn << endl;         
     }
   }   
 }

}

// Build on demand, ie if is empty
void Foam::sparseSolverBase::buildSharedDataOnDemand
( 
  const fvMesh& mesh
)
{ 
  if (sharedData.procInfo.size() == 0) 
    setSharedData(mesh);
}

// ************************************************************************* //
