PtrList<fvMesh> blockcRegions(rp.blockcRegionNames().size());

forAll(rp.blockcRegionNames(), zoneID)
{
    Info<< nl << rp.blockcRegionNames()[zoneID] << endl;

    blockcRegions.set
    (
        zoneID,
        new fvMesh
        (
            IOobject
            (
                rp.blockcRegionNames()[zoneID],
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        )
    );
}
