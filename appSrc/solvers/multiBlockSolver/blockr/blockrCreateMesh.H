const wordList blockrRegionNames(rp["blockrRegionNames"]);

PtrList<fvMesh> blockrRegions(blockrRegionNames.size());

forAll(blockrRegionNames, zoneID)
{
    Info<< nl << blockrRegionNames[zoneID] << endl;

    blockrRegions.set
    (
        zoneID,
        new fvMesh
        (
            IOobject
            (
                blockrRegionNames[zoneID],
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        )
    );
}
