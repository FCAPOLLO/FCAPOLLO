forAll(blockcRegions, zoneID)
{
 //  	Info<< "Solving Transport in the " << blockcRegions[zoneID].name() << endl;
    #include "blockcSetFields.H"
	#include "blockcCalcFields.H"
    #include "blockcReadSIMPLEControls.H"
    #include "protonTransportEquation.H"
    #include "blockcConvergence.H"
}
