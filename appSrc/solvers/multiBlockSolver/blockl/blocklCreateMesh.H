const wordList blocklRegionNames(rp["blocklRegionNames"]);

PtrList<fvMesh> blocklRegions(blocklRegionNames.size());

forAll(blocklRegionNames, zoneID)
{
	Info<< nl << blocklRegionNames[zoneID] << endl;

    blocklRegions.set
    (
        zoneID,
        new fvMesh
        (
            IOobject
            (
                blocklRegionNames[zoneID],
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        )
    );
}
