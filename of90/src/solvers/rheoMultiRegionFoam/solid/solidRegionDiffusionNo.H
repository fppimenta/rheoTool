scalar DiNum = -great;

forAll(solidRegions, i)
{
    const solidThermoModel& thermo = thermoSolid[i];

    scalar maxDiffNo;
    scalar meanDiffNo;
    
    thermo.DiffNo(solidRegions[i],runTime,maxDiffNo,meanDiffNo);

    DiNum = max
    (
        maxDiffNo,
        DiNum
    );
}
