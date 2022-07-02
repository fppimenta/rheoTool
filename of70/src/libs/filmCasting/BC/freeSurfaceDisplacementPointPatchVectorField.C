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

#include "freeSurfaceDisplacementPointPatchVectorField.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"
#include "polyMesh.H"
#include "fvsPatchField.H"
#include "surfaceFields.H"
#include "fvPatchField.H"
#include "fvcGrad.H"
  
#include "scalarMatrices.H"
#include "SortableList.H"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <functional>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

freeSurfaceDisplacementPointPatchVectorField::
freeSurfaceDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(p, iF)
{}


freeSurfaceDisplacementPointPatchVectorField::
freeSurfaceDisplacementPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict),
    anchorPoints_(2, -1),
    outletPoints_(2, -1),
    stencil_(p.size(), List<label>(2, -1)),
    method_(dict.lookup("method")), 
    limiterFunc_(method_ == "streamline" ? dict.lookup("limiterFunction") : word(" ")),  
    URF_(readScalar(dict.lookup("URF"))),
    useFlux_(method_ == "streamline" ? readBool(dict.lookup("useFlux")) : false),
    y0_(p.localPoints().component(1)),
    auxF_(this->internalField().mesh()().boundaryMesh()[p.index()].faceCentres()),
    motionDir_(method_ == "streamline" ? vector::zero : dict.lookup("motionDir"))
{   
    if (method_ != "Peric" && method_ != "streamline")
    {
      FatalErrorIn("freeSurfaceDisplacementPointPatchVectorField::freeSurfaceDisplacementPointPatchVectorField\n")
            << "\n Error: unknown 'method' in patch "<< p.name()
            << "\n Available methods are: \n"
            << "     o   Peric \n"
            << "     o   streamline \n"
            << "(Take care with lower- and upper-cases)\n"
            << abort(FatalError);
    }

    buildStencil();
 
    if (dict.found("points0"))
    {
       pts0_ = vectorField("points0", dict , p.size());
    }
    else
    {
       pts0_ = p.localPoints();
    }
     

    if (dict.found("value"))
    {
        dX_ = (*this);
    }
    else
    {
        dX_ = vectorField(p.size(), vector::zero);
        updateCoeffs();
    }
}
 


freeSurfaceDisplacementPointPatchVectorField::
freeSurfaceDisplacementPointPatchVectorField
(
    const freeSurfaceDisplacementPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper)
{}


freeSurfaceDisplacementPointPatchVectorField::
freeSurfaceDisplacementPointPatchVectorField
(
    const freeSurfaceDisplacementPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void freeSurfaceDisplacementPointPatchVectorField::buildStencil
( 
)
{  
 const pointPatch& p = this->patch();

 if (p.size()>0)
 { 
   label np = p.size();
 
   // Anchor and outlet points (boundary points)
   anchorPoints_[0] = 0;
   anchorPoints_[1] = 1;
   
   outletPoints_[0] = np - 1;
   outletPoints_[1] = np - 2;
   
   // Inner points. 
   // Stencil [previous_point next_point] or s_i = [p_i-1 p_i+1]
   // The mesh may have one or two free-surfaces, depending on whether symmetry is applied or not. Point numbering is 
   // different in each case:
   // -> Free-surface is in positive quadrant
   if (p.localPoints()[0].y() > 0)
   {
    //- Special points
    stencil_[0][0] = -1; stencil_[0][1] = 3;   
    stencil_[1][0] = -1; stencil_[1][1] = 2;
   
    stencil_[2][0] = 1; stencil_[2][1] = 4;   
    stencil_[3][0] = 0; stencil_[3][1] = 5;
   
    stencil_[np-1][0] = np-3; stencil_[np-1][1] = -1;   
    stencil_[np-2][0] = np-4; stencil_[np-2][1] = -1;
   
    //- Other points
    for (label i=4; i<np-2; i++)
    {
      stencil_[i][0] = i-2;
      stencil_[i][1] = i+2; 
    }
   }
   // -> Free-surface is in negative quadrant
   else
   {
    //- Special points
    stencil_[0][0] = -1; stencil_[0][1] = 1;   
    stencil_[3][0] = -1; stencil_[3][1] = 2;
   
    stencil_[2][0] = 3; stencil_[2][1] = 5;   
    stencil_[1][0] = 0; stencil_[1][1] = 4;
    
    stencil_[5][0] = 2; stencil_[5][1] = 7;   
    stencil_[4][0] = 1; stencil_[4][1] = 6;
   
    stencil_[np-1][0] = np-3; stencil_[np-1][1] = -1;   
    stencil_[np-2][0] = np-4; stencil_[np-2][1] = -1;
   
    //- Other points
    for (label i=6; i<np-2; i++)
    {
      stencil_[i][0] = i-2;
      stencil_[i][1] = i+2; 
    }
   } 
 }
}



void freeSurfaceDisplacementPointPatchVectorField::interpolateUonPoints
(
  vectorField& Up,
  const scalarField& absPhiP
)
{
 const pointPatch& p = this->patch();
 const polyMesh& mesh = this->internalField().mesh()();
 const polyPatch& pPatch_(this->internalField().mesh()().boundaryMesh()[p.index()]);   
 const vectorField::subField A = pPatch_.faceAreas();
 const vectorField::subField Cf = pPatch_.faceCentres();
 const vectorField& pts = p.localPoints();
 const labelListList& pfs = mesh.boundaryMesh()[p.index()].pointFaces();
 const fvPatchField<vector>& Upatch = db().lookupObject<volVectorField>("U").boundaryField()[p.index()];  

 forAll(p, pi)
 {   
   labelList faces = pfs[pi];
   scalar sumW(0.);
   vector sumU(vector::zero);
            
   forAll(faces, fi)
   {
     vector dis = Cf[faces[fi]] - pts[pi];
        
     vector Uf = Upatch[faces[fi]];
     if (useFlux_)
     {
       vector n = A[faces[fi]]/mag(A[faces[fi]]);
       Uf = Uf - (Uf&n)*n + absPhiP[faces[fi]]*n/mag(A[faces[fi]]); 
     }
      
     sumU += Uf/mag(dis);        
     sumW += 1./mag(dis);            
   }
      
   Up[pi] = sumU/sumW;
 }
 
 // Note: for the points on the outlet we know the velocity of their face on the outlet patch.
 // We can use it to make it a bit more accurate.
 if (0 == 1)
 {
   label outID = mesh.boundaryMesh().findPatchID("outlet"); 
   const polyPatch& pPatchOut_(this->internalField().mesh()().boundaryMesh()[outID]); 
   const vectorField::subField Cfout = pPatchOut_.faceCentres();   
   scalar maxY(-1);
   label pyMax(-1);
   forAll(Cfout, i)
   {
     if (mag(Cfout[i].y()) > maxY)
     {
       maxY = mag(Cfout[i].y());
       pyMax = i;
     }
   }
 
   const fvPatchField<vector>& Uoutpatch = db().lookupObject<volVectorField>("U").boundaryField()[outID]; 
   vector uPout = Uoutpatch[pyMax];
 
   Up[outletPoints_[0]] = uPout;
   Up[outletPoints_[1]] = uPout;     
 }
}

void freeSurfaceDisplacementPointPatchVectorField::integrate
(
  scalarField& y,
  const vectorField& Up,
  const vectorField& pts,
  const scalar& dt,
  scalar& maxCoRef,
  scalar& meanCoRef
)
{
 // Explicit Lax-Wendroff method with limiter (https://folk.ntnu.no/leifh/teaching/tkt4140/._main074.html#___sec158)
 
 // Define limiter
 std::function<scalar(const scalar&)> limiter;
 if (limiterFunc_ == "superBee")
 {
  limiter = [](const scalar& r){return max(max(min(2.*r, 1.), min(r, 2.)), 0);};
 }
 else if (limiterFunc_ == "vanLeer")
 {
  limiter = [](const scalar& r){return (r + mag(r))/(1.+mag(r));};
 }
 else if (limiterFunc_ == "LW")
 {
  limiter = [](const scalar& r){return 1.;};
 }
 else if (limiterFunc_ == "minMod")
 {
  limiter = [](const scalar& r){return max(min(r, 1.), 0.);};
 }
 else if (limiterFunc_ == "upwind")
 {
  limiter = [](const scalar& r){return 0;};
 }
 else
 {
   FatalErrorIn("limiters\n")
     << "\nUnknown limiter." 
     << abort(FatalError);
 }
 
 // Init Courant vars
 scalar maxCo(0);
 scalar meanCo(0);
 label  allcount(0); 
 
 // Loop over points: skip first points (cntPt0 and cntPt1) as they are anchored
 // Take into account if freeSurface 0 or 1 in points numbering
 label np(pts.size());
 label upwindPt0(pts0_[0].y() > 0 ? 2 : 2 );
 label upwindPt1(pts0_[0].y() > 0 ? 3 : 1 );
 label cntPt0(pts0_[0].y() > 0 ? 0 : 0 );
 label cntPt1(pts0_[0].y() > 0 ? 1 : 3 );
 
 for (int i = 0; i < np; i++)
 {
    if (i == cntPt0 || i == cntPt1)
     continue;
 
    scalar rpre(0.), rpos(0.);
     
    // Force upwind on second and last cols.
    if (i == upwindPt0 || i == upwindPt1 || i == np - 1 || i == np - 2)
    {
      rpos = 0; rpre = 0;
    }
    else
    {
      rpre = (y0_[stencil_[i][0]] - y0_[stencil_[stencil_[i][0]][0]])/(y0_[i] - y0_[stencil_[i][0]] + 1e-12);       
      rpos = (y0_[i] - y0_[stencil_[i][0]])/(y0_[stencil_[i][1]] - y0_[i] + 1e-12);
    }
   
    scalar phipre(limiter(rpre));  
    scalar phipos(limiter(rpos));  
  
    scalar Fpos(0.), Fpre(0.);
    scalar dx(0.);
    // Last two points need special attention as they do not have i+1 points
    if (i == np - 1 || i == np - 2)
    {
     // upwind
     dx = pts[i].x() - pts[stencil_[i][0]].x();
     scalar dydx = (pts[i].y() - pts[stencil_[i][0]].y())/dx;   
     y[i] = y0_[i] + dt*(Up[i].y() - Up[i].x()*dydx); 
   }
   else
   {
     label ipos(stencil_[i][1]);
     label iprev(stencil_[i][0]);
     
     // Pos
     dx = pts[ipos].x() - pts[i].x();
     Fpos = pts[i].y() + phipos*0.5*(1.-dt*Up[i].x()/dx)*(pts[ipos].y()-pts[i].y());
     
     // Prev
     dx = pts[i].x() - pts[iprev].x();
     Fpre = pts[iprev].y() + phipre*0.5*(1.-dt*Up[i].x()/dx)*(pts[i].y()-pts[iprev].y());
      
      
     dx = 0.5*(pts[ipos].x() + pts[i].x()) - 0.5*(pts[iprev].x() + pts[i].x());  
     y[i] = y0_[i] + (dt*Up[i].x()/dx)*(Fpre-Fpos) + dt*Up[i].y();  
   }
    
   dX_[i] = vector(0., y[i]-pts0_[i].y(), 0.); 
    
   // Courant (it excludes the first two points)
   scalar Coi = mag(Up[i])*dt/dx;
   allcount++;
   meanCo += Coi;
   if (Coi>maxCo)
     maxCo = Coi;   
 }
 
 meanCo /= allcount;
 
 maxCoRef = maxCo;
 meanCoRef = meanCo;
}

//- Solve the equation of the streamline
// Note: phi must absolute here   
void freeSurfaceDisplacementPointPatchVectorField::updatePoints
(
   const bool& isLastIter,
   const scalarField& absPhiP
)
{   
 scalar maxCo(0);
 scalar meanCo(0);
 const pointPatch& p = this->patch();
 if (p.size() > 0)
 { 
    const polyMesh& mesh = this->internalField().mesh()();
    const Time& t = mesh.time();
    const vectorField& pts = p.localPoints();
    scalar dt = t.deltaT().value();  
      
    // Interpolate velocity from faces to points
    vectorField Up(p.size(), vector::zero);
    
    interpolateUonPoints
    (
      Up,
      absPhiP
    );
      
    // Integrate the 2 streamlines 
    scalarField y(y0_);
    integrate(y, Up, pts, dt, maxCo, meanCo);
   
    // Kind of explicit under-relaxation (only acts on y0_, not on y, as dX_ is not touched)
    if (URF_ < 1-1e-6 && URF_ > 1e-6)
    {
      y = URF_*y + (1.-URF_)*y0_;
    }
    
    // Note: anchor points are not touched during integration, so no special treatment is needed 
     
    // Save current positions to be used in next time-step as previous time position. To be done 
    // on the last iteration of the inner loop. 
    if (isLastIter)
      y0_ = y;   
  }  
   
  reduce(maxCo, maxOp<scalar>() );
  reduce(meanCo, maxOp<scalar>() );
  Info << "Point Courant at " << this->patch().name() << tab << "Max: " << maxCo << tab << "Mean: " << meanCo << endl; 
}

//- Peric's method: http://dx.doi.org/10.1080/10407799708915014
// Note: phi must be the flux relative to mesh motion.   
void freeSurfaceDisplacementPointPatchVectorField::updatePoints
(
  const scalarField& phiRel
)
{   
  const polyMesh& mesh = this->internalField().mesh()();
  const Time& t = mesh.time();
  const pointPatch& p = this->patch();
  
  if (p.size() > 0)
  {    
   scalar dt = t.deltaT().value();
  
   const polyPatch& pPatch_(this->internalField().mesh()().boundaryMesh()[p.index()]);
    
   const vectorField::subField A = pPatch_.faceAreas();
   const vectorField::subField Cf = pPatch_.faceCentres();
   const labelListList& pfs = mesh.boundaryMesh()[p.index()].pointFaces();
  
   // Moving direction    
   vectorField dir;
   if (mag(motionDir_)<1e-200)
   //- Dynamic (normal to faces)
   {   
     dir = A/mag(A);
     // Need to ensure that outlet points only move in y-dir in case 'dir' is the face normal
     forAll(outletPoints_, pi)
     {    
       labelList faces = pfs[outletPoints_[pi]];
       forAll(faces, fi)
         dir[faces[fi]] = vector(0,1,0);
     }
   }
   else
   //- Fixed
   {
     dir = vectorField(A.size(), motionDir_); 
   }
   
    
   scalarField dh = URF_*phiRel*dt/( (A&dir) + SMALL); 
 
   auxF_ += dh * dir; 
     
   vectorField pts(p.localPoints());   
   forAll(p, pi)
   {    
     labelList faces = pfs[pi];
     vector wV(vector::zero);
     scalar sumW(0.);
     vector dir_av(vector::zero);
     forAll(faces, fi)
     {
       scalar w = 1./mag(Cf[faces[fi]] - pts[pi]) ;
       wV += auxF_[faces[fi]]*w;
       dir_av += dir[faces[fi]]*w;
       sumW += w;            
     }
    
     wV /= sumW;
     dir_av /= sumW;      
      
     pts[pi] -= dir_av*(dir_av & (pts[pi] - wV));     
   }
    
   dX_ = pts - pts0_;
  
   // Anchor points
   forAll(anchorPoints_, i)
     dX_[anchorPoints_[i]] = vector::zero;
  }
}

void freeSurfaceDisplacementPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
      return;
    }
 
    Field<vector>::operator=(dX_);

    fixedValuePointPatchField<vector>::updateCoeffs();
}


void freeSurfaceDisplacementPointPatchVectorField::write(Ostream& os) const
{
    pointPatchField<vector>::write(os);
    
    os.writeKeyword("method") << method_ << token::END_STATEMENT << nl;
        
    if (method_ == "streamline")
    {
      writeEntry(os, "useFlux", useFlux_);
      os.writeKeyword("limiterFunction") << limiterFunc_ << token::END_STATEMENT << nl;
    }
    else
    {
      writeEntry(os, "motionDir", motionDir_);
    }
 
    writeEntry(os, "URF", URF_);
    writeEntry(os, "points0", pts0_);
    writeEntry(os, "value", *this);
}

 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    freeSurfaceDisplacementPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
