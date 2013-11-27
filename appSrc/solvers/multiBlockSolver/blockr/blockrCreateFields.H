// Initialise electron field pointer lists
#include "blockrFieldPointerLists.H"

// Populate blockr field pointer lists
forAll(blockrRegions, zoneID)
{
	Info<< nl << "    Reading transportProperties file in constants/" << blockrRegionNames[zoneID] << endl;
    // Read in all related dictionaries
    #include "blockrReadComponentDictionary.H"
	
    // Assign any Dimensioned Scalars for input values or field initialisation
    #include "blockrDimScalarAssign.H"
	
    // Create Field Objects for solved variables
    #include "blockrSolvedFieldObjects.H"

    // Create Field Objects for transport coefficient variables
    #include "blockrCoefficientFieldObjects.H"

    // Create Field Objects for source term variables
    #include "blockrSourceTermFieldObjects.H"

	// Create Field Objects for calculated variables
	#include "blockrCalculatedFieldObjects.H"
}
