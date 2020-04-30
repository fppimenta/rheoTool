// Solve coupled T system

if (solveCoupledT)
{
 Info << nl << "Solving T coupled for all mesh regions" << endl; 
 cpsT->solve();
}
