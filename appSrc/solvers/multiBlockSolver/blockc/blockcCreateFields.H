// Initialise electron field pointer lists
#include "blockcFieldPointerLists.H"

// Populate blockc field pointer lists
forAll(blockcRegions, zoneID)
{
	Info<< nl << "    Reading transportProperties file in constants/" << blockcRegionNames[zoneID] << endl;
    // Read in all related dictionaries
    #include "blockcReadComponentDictionary.H"
	
    // Assign any Dimensioned Scalars for input values or field initialisation
    #include "blockcDimScalarAssign.H"
	
    // Create Field Objects for solved variables
    #include "blockcSolvedFieldObjects.H"

    // Create Field Objects for transport coefficient variables
    #include "blockcCoefficientFieldObjects.H"

    // Create Field Objects for source term variables
    #include "blockcSourceTermFieldObjects.H"

	// Create Field Objects for the calculated variables
	#include "blockcCalculatedFieldObjects.H"

}
