// Read in all related dictionaries
IOdictionary transportProperties
(
    IOobject
    (
		"transportProperties",
		runTime.constant(),
		blockcRegions[zoneID],
		IOobject::MUST_READ_IF_MODIFIED,
		IOobject::NO_WRITE
	)
);
