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

\*---------------------------------------------------------------------------*/

#include "spinesFvMotionSolver.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"
#include "OFstream.H"
#include "meshTools.H"
#include "mapPolyMesh.H"
#include "volPointInterpolation.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(spinesFvMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        spinesFvMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::spinesFvMotionSolver::spinesFvMotionSolver
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName),
    fvMotionSolver(mesh),
    pointLocation_(nullptr),
    frozenPointsZone_
    (
        coeffDict().found("frozenPointsZone")
      ? fvMesh_.pointZones().findZoneID(coeffDict().lookup("frozenPointsZone"))
      : -1
    ),
    pID_(fvMesh_.points().size(), 0),
    globalAdr_(nullptr),
    oneFS_(true)
{
    IOobject io
    (
        "pointLocation",
        fvMesh_.time().timeName(),
        fvMesh_,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    if (debug)
    {
        Info<< "spinesFvMotionSolver:" << nl
            << "    frozenPoints zone : " << frozenPointsZone_ << endl;
    }


    if (io.typeHeaderOk<pointVectorField>(true))
    {
        pointLocation_.reset
        (
            new pointVectorField
            (
                io,
                pointMesh::New(fvMesh_)
            )
        );

        if (debug)
        {
            Info<< "spinesFvMotionSolver :"
                << " Read pointVectorField "
                << io.name()
                << " to be used for boundary conditions on points."
                << nl
                << "Boundary conditions:"
                << pointLocation_().boundaryField().types() << endl;
        }
    }
       
    markPatches();
    
    if (Pstream::parRun())
    {
       // Global adress of points
       globalAdr_.reset
       (
         new labelIOList
         (
          IOobject
          (
	    "pointProcAddressing",
	     fvMesh_.facesInstance(),
	     fvMesh_.meshSubDir,
	     fvMesh_,
	     IOobject::MUST_READ,
	     IOobject::NO_WRITE,
	     false
          )
         )
       );
       
       // Define number of points. 
       // npz: easy
       // npx: freeSurface is in a single processor, so get its size
       // npy: needs npTot, which we can get by looking into the highest global index (+1 to get size)  
       npz_ = 2;
       
       label npTot(max(globalAdr_()));
       reduce(npTot, maxOp<label>() );
       npTot+=1;
       
       label npx(0);
       forAll(mesh.boundaryMesh(), i)
       {
         word pName(mesh.boundaryMesh()[i].name());
         if (pName == "freeSurface" || pName == "freeSurface0" || pName == "freeSurface1")
         {
           npx = pointDisplacement_.boundaryField()[i].patch().meshPoints().size();
           break;  
         }
       }
       
       reduce(npx, sumOp<label>() );
       npx_ = npx/2;
       
       npy_ = npTot/(2*npx_);
       
       // Make the local ordering for freeSurface (list is non-zero in a single proc) 
       // The free-surface patch is contained in a single processor. Get the global adressing
       // for the points on that patch and make a list with them. The list is ordered in increasing
       // order of x of the corresponding points, such element [0] is the pt shared with inlet and 
       // element[end] is the element shared with outlet. This is practical, since when we exchange
       // between processors the y-coord of these free surface points in runtime, we already know the
       // ordering.
       forAll(mesh.boundaryMesh(), i)
       {
         word pName(mesh.boundaryMesh()[i].name());
         if (pName == "freeSurface0")
         {
           if (mesh.boundaryMesh()[i].size()>0)
           {
             // Index j will provide the global adressing, i.e. the points numbering in the non-decomposed
             // mesh, which is sequential.
             label j0=npx_*npy_-npx_;
             label jend=npx_*npy_;

             for (int j = j0; j<jend; j++)
             {
               forAll(globalAdr_(), k)
               {
                 if (globalAdr_()[k] == j)
                 {
                   fsOrder0_.append(k);
                   break;
                 }
               }
             } 
           }
         }
         else if (pName == "freeSurface1")
         {
           oneFS_ = false;
           if (mesh.boundaryMesh()[i].size()>0)
           {
             // Index j will provide the global adressing, i.e. the points numbering in the non-decomposed
             // mesh, which is sequential.
             label j0=0;
             label jend=npx_;

             for (int j = j0; j<jend; j++)
             {
               forAll(globalAdr_(), k)
               {
                 if (globalAdr_()[k] == j)
                 {
                   fsOrder1_.append(k);
                   break;
                 }
               }
             } 
           }
         }
       }         
    }
}
 
    
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::spinesFvMotionSolver::
~spinesFvMotionSolver()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
 

void Foam::spinesFvMotionSolver::markPatches()
{
  // pID:
  // 0 for interior points and any patch not in the list below;
  // 1 for outlet
  // 2 for symm
  // 3 for inlet
  // >3 for freeSurface(s)
  
  //---
  // outlet
  {
    label patchID = fvMesh_.boundaryMesh().findPatchID("outlet");
    if (patchID > -1)
    {
     const labelList& patchpts = pointDisplacement_.boundaryField()[patchID].patch().meshPoints();
     forAll(patchpts, j)
     {
       pID_[patchpts[j]] = 1;  
     }
    } 
  }

  //---- Priority (untouchable) patches come last to overwrite pID
  
  // symm
  {
    label patchID = fvMesh_.boundaryMesh().findPatchID("symm");
    if (patchID > -1)
    {
     const labelList& patchpts = pointDisplacement_.boundaryField()[patchID].patch().meshPoints();
     forAll(patchpts, j)
     {
       pID_[patchpts[j]] = 2;  
     }
    } 
  }
  
  // inlet
  {
    label patchID = fvMesh_.boundaryMesh().findPatchID("inlet");
    if (patchID > -1)
    {
     const labelList& patchpts = pointDisplacement_.boundaryField()[patchID].patch().meshPoints();
     forAll(patchpts, j)
     {
       pID_[patchpts[j]] = 3;  
     } 
    
     npy_ = patchpts.size()/2.;
    }
  }
   
  // freeSurface
  {
   forAll(fvMesh_.boundaryMesh(), i)
   {
    word pName(fvMesh_.boundaryMesh()[i].name());
    if (pName == "freeSurface" || pName == "freeSurface0" || pName == "freeSurface1")
    {
     const labelList& patchpts = pointDisplacement_.boundaryField()[i].patch().meshPoints();
     forAll(patchpts, j)
     {
      pID_[patchpts[j]] = 4;  
     }
    }
   } 
  }
  
  // Remaining np
  npz_ = 2;
  npx_ = pID_.size()/(2*npy_);
}

Foam::tmp<Foam::pointField> Foam::spinesFvMotionSolver::curPointsSequential() const
{ 
  const pointField& pts = fvMesh_.points();
 
  forAll(pID_, i)
  {    
    // Only set free points
    if (pID_[i] < 2)
    {
      label zi = i/(npy_*npx_);
      label a = i-npx_*npy_*zi;
      label yi = a/npx_;
      label xi = a - yi*npx_;
      
      label p0 = xi;
      label pend = npx_*npy_ - npx_ + xi;
      
      scalar frac0 = (points0()[i].y()-points0()[p0].y())/(points0()[pend].y()-points0()[p0].y());
      scalar y = pts[p0].y()+frac0*(pts[pend].y()-pts[p0].y());
      
      pointDisplacement_[i] = vector(0,y-points0()[i].y(),0);
    }  
  }
 
  tmp<pointField> tcurrPts(new pointField(pts));
  pointField& curPoints = tcurrPts.ref();
   
  forAll(pointDisplacement_, i)
  { 
     curPoints[i] += points0()[i] + pointDisplacement_[i] - pts[i]; // Subtract pts because tcurrPts is initialized with pts
  }
   
  twoDCorrectPoints(curPoints);
  
  pointDisplacement_.primitiveFieldRef() = curPoints - points0();
  
  pointDisplacement_.correctBoundaryConditions();

  return tcurrPts;   
}


Foam::tmp<Foam::pointField> Foam::spinesFvMotionSolver::curPointsParallel() const
{ 
  const pointField& pts = fvMesh_.points();
  
  // yfsi is a vector used as container (list) where only the first 2 elements are important
  // yfsi.x() - the actual y-coordinate of the free-surface for each point on that patch
  // yfsi.y() - initial y-coordinate of the free-surface for each point on that patch
  // The loop only runs in the procs where fsOrder_ is not empty, i.e. in a single proc
  
  // yfs0 always exists, as there is always, at least, patch freeSurface0
  List<vector> yfs0;
  {
   List<vector> yfsi;
   forAll(fsOrder0_, i)
   {
      vector t(pts[fsOrder0_[i]].y(), points0()[fsOrder0_[i]].y(), 0);
      yfsi.append(t);
   }
  
   // gather-scatter-combine. Some overhead (n to n-1 should be just 1 to n-1).
   List< List<vector>  > allyfs(Pstream::nProcs());     
   allyfs[Pstream::myProcNo()] = yfsi;  
   Pstream::gatherList(allyfs);           
   Pstream::scatterList(allyfs);
     
   yfs0 = ListListOps::combine< List<vector> >
   (
     allyfs,
     accessOp< List<vector> >()
   );
  }

  // Repeat process if there is another freeSurface (no symmetry plane)
  List<vector> yfs1;
  if (!oneFS_)
  {
   List<vector> yfsi;
   forAll(fsOrder1_, i)
   {
      vector t(pts[fsOrder1_[i]].y(), points0()[fsOrder1_[i]].y(), 0);
      yfsi.append(t);
   }

   // gather-scatter-combine. Some overhead (n to n-1 should be just 1 to n-1).
   List< List<vector>  > allyfs(Pstream::nProcs());     
   allyfs[Pstream::myProcNo()] = yfsi;  
   Pstream::gatherList(allyfs);           
   Pstream::scatterList(allyfs);
      
   yfs1 = ListListOps::combine< List<vector> >
   (
     allyfs,
     accessOp< List<vector> >()
   );
  }

  // Proceed as if sequential   
  forAll(pID_, i)
  {    
    // Only set free points
    if (pID_[i] < 2)
    {
      label ig = globalAdr_()[i];
      label zi = ig/(npy_*npx_);
      label a = ig-npx_*npy_*zi;
      label yi = a/npx_;
      label xi = a - yi*npx_;
      
      // Equivalence with sequential:
      // yfs0[xi].y() = points0()[pend].y() 
      // yfs0[xi].x() = pts[pend].y()
      // yfs1[xi].y() = points0()[p0].y() 
      // yfs1[xi].x() = pts[p0].y()
      if (oneFS_)
      {
       // If there is only a single free-surface, we assume implicitly that  points0()[p0].y() = pts()[p0].y(),
       // i.e. the bottom of the spine is in the symmetry axis, so y = 0
       scalar frac0 = points0()[i].y()/yfs0[xi].y();  
       scalar y = frac0*yfs0[xi].x();                 
      
       pointDisplacement_[i] = vector(0,y-points0()[i].y(),0);
      }
      else
      {     
       // Each director of the spine is in a different freeSurface. 
       scalar frac0 = (points0()[i].y()-yfs1[xi].y())/(yfs0[xi].y()-yfs1[xi].y());
       scalar y = yfs1[xi].x()+frac0*(yfs0[xi].x()-yfs1[xi].x());
                            
       pointDisplacement_[i] = vector(0,y-points0()[i].y(),0);     
      }
    }  
  }

  tmp<pointField> tcurrPts(new pointField(pts));
  pointField& curPoints = tcurrPts.ref();
   
  forAll(pointDisplacement_, i)
  { 
     curPoints[i] += points0()[i] + pointDisplacement_[i] - pts[i]; // Subtract pts because tcurrPts is initialized with pts
  }
   
  twoDCorrectPoints(curPoints);
  
  pointDisplacement_.primitiveFieldRef() = curPoints - points0();
  
  pointDisplacement_.correctBoundaryConditions();
  
  return tcurrPts;   
}
 
Foam::tmp<Foam::pointField>
Foam::spinesFvMotionSolver::curPoints() const
{ 
    pointDisplacement_.boundaryFieldRef().updateCoeffs(); 
    
    if (Pstream::parRun())
    {   
      return curPointsParallel();     
    }
    else
    {   
      return curPointsSequential();     
    }
}


void Foam::spinesFvMotionSolver::solve()
{
  // The points have moved so before interpolation update
  // the motionSolver accordingly
  movePoints(fvMesh_.points());
}


void Foam::spinesFvMotionSolver::updateMesh
(
    const mapPolyMesh& mpm
)
{
  displacementMotionSolver::updateMesh(mpm);
}


// ************************************************************************* //
