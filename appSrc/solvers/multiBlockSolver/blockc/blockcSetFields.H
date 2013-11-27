// Create assignment to the mesh object
	fvMesh& mesh = blockcRegions[zoneID];

// Read in the transport Properties dictionary if it has been modified since last read
#include "blockcReadComponentDictionary.H"

// Assignment the value in the dictionary to the electronic conductivity
#include "blockcDimScalarAssign.H"
	// Create Assignments to Calculated Field Variables
	surfaceScalarField& protonFlux = proton_flux_blockc[zoneID];
	// Create Assignments to Solved Field Variables
	// Electron Transport Fields
	// Proton Transport Fields
	volScalarField& phi_p = phi_p_blockc[zoneID];
// Create Assignments to Coefficient Field Variables
	// Electron Transport Fields
	// Proton Transport Fields
	volTensorField& sigma_p_eff = sigma_p_eff_blockc[zoneID];
// Create Assignments to Source/Sink Field Variables
	// Electron Transport Fields
	volScalarField& source_proton_total = source_proton_total_blockc[zoneID];
	volScalarField& source_proton_SLOPE = source_proton_SLOPE_blockc[zoneID];
