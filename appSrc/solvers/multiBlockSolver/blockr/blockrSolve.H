forAll(blockrRegions, zoneID)
{
 //  	Info<< "Solving Transport in the " << blockrRegions[zoneID].name() << endl;
    #include "blockrSetFields.H"
	#include "blockrCalcFields.H"
    #include "blockrReadSIMPLEControls.H"
    #include "electronTransportEquation.H"
	#include "protonTransportEquation.H"
    #include "blockrConvergence.H"
}
