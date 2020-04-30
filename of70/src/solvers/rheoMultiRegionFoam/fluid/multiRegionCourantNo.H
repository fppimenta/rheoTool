    scalar CoNum = -great;

    forAll(fluidRegions, regionI)
    {
        CoNum = max
        (
            iCourantNo
            (
                fluidRegions[regionI],
                runTime,
                phiFluid[regionI]
            ),
            CoNum
        );
    }
