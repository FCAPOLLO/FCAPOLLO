// Create assignment to the mesh object
	fvMesh& mesh = blocklRegions[zoneID];

//	labelList NCellsI(mesh.cellZones());

	labelList NCellsI(mesh.C().size());

// Read in the transport Properties dictionary if it has been modified since last read
#include "blocklReadComponentDictionary.H"

// Assignment the value in the dictionary to the electronic conductivity
#include "blocklDimScalarAssign.H"

	// Create Assignments to Calculated Field Variables
	volScalarField& eta_act = eta_act_blockl[zoneID];
	surfaceScalarField& protonFlux = proton_flux_blockl[zoneID];
	surfaceScalarField& electronFlux = electron_flux_blockl[zoneID];

	// Create Assignments to Solved Field Variables
	// Electron Transport Fields
	volScalarField& phi_e = phi_e_blockl[zoneID];
	// Proton Transport Fields
	volScalarField& phi_p = phi_p_blockl[zoneID];

// Create Assignments to Coefficient Field Variables
	// Electron Transport Fields
	volTensorField& sigma_e_eff = sigma_e_eff_blockl[zoneID];
	// Proton Transport Fields
	volTensorField& sigma_p_eff = sigma_p_eff_blockl[zoneID];
// Create Assignments to Source/Sink Field Variables
	// Electron Transport Fields
	volScalarField& reaction_rate_THV_ANODIC = reaction_rate_THV_ANODIC_blockl[zoneID];
	volScalarField& reaction_rate_THV_CATHODIC = reaction_rate_THV_CATHODIC_blockl[zoneID];
	volScalarField& reaction_rate_THV_SUM = reaction_rate_THV_SUM_blockl[zoneID];
	volScalarField& reaction_rate_BV_ANODIC = reaction_rate_BV_ANODIC_blockl[zoneID];
	volScalarField& reaction_rate_BV_CATHODIC = reaction_rate_BV_CATHODIC_blockl[zoneID];
	volScalarField& reaction_rate_BV_SUM = reaction_rate_BV_SUM_blockl[zoneID];

	volScalarField& reaction_rate_THV_SLOPE_ELECTRON = reaction_rate_THV_SLOPE_ELECTRON_blockl[zoneID];
	volScalarField& reaction_rate_THV_SLOPE_PROTON = reaction_rate_THV_SLOPE_PROTON_blockl[zoneID];

	volScalarField& reaction_rate_BV_SLOPE_ELECTRON = reaction_rate_BV_SLOPE_ELECTRON_blockl[zoneID];
	volScalarField& reaction_rate_BV_SLOPE_PROTON = reaction_rate_BV_SLOPE_PROTON_blockl[zoneID];

	volScalarField& source_electron_total = source_electron_total_blockl[zoneID];
	volScalarField& source_proton_total = source_proton_total_blockl[zoneID];

	volScalarField& source_electron_SLOPE = source_electron_SLOPE_blockl[zoneID];
	volScalarField& source_proton_SLOPE = source_proton_SLOPE_blockl[zoneID];
