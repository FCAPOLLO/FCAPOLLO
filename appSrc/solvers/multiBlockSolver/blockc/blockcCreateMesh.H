const wordList blockcRegionNames(rp["blockcRegionNames"]);

PtrList<fvMesh> blockcRegions(blockcRegionNames.size());

forAll(blockcRegionNames, zoneID)
{
    Info<< nl << blockcRegionNames[zoneID] << endl;

    blockcRegions.set
    (
        zoneID,
        new fvMesh
        (
            IOobject
            (
                blockcRegionNames[zoneID],
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        )
    );
}
