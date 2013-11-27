// Initialise electron field pointer lists
#include "blocklFieldPointerLists.H"

// Populate blockl field pointer lists
forAll(blocklRegions, zoneID)
{
	Info<< nl << "    Reading transportProperties file in constants/" << blocklRegionNames[zoneID] << endl;
    // Read in all related dictionaries
    #include "blocklReadComponentDictionary.H"
	
    // Assign any Dimensioned Scalars for input values or field initialisation
    #include "blocklDimScalarAssign.H"
	
    // Create Field Objects for solved variables
    #include "blocklSolvedFieldObjects.H"

    // Create Field Objects for transport coefficient variables
    #include "blocklCoefficientFieldObjects.H"

    // Create Field Objects for source term variables
    #include "blocklSourceTermFieldObjects.H"

	// Create Field Objects for the calculated variables
	#include "blocklCalculatedFieldObjects.H"
}
