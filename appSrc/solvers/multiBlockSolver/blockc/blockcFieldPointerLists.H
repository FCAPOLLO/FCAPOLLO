// Initilialise field pointer lists
//
// Calculated Field Variables
	PtrList<surfaceScalarField> proton_flux_blockc(blockcRegions.size());
// Solved Field Variables	
	// Electron Transport Fields
	// Proton Transport Fields
	PtrList<volScalarField> phi_p_blockc(blockcRegions.size());

// Coefficient Field Variables
	// Electron Transport Fields
	// Proton Transport Fields
	PtrList<volTensorField> sigma_p_eff_blockc(blockcRegions.size());
// Source/sink Field Variables
	// Electron Transport Fields
	PtrList<volScalarField> source_proton_total_blockc(blockcRegions.size());
	PtrList<volScalarField> source_proton_SLOPE_blockc(blockcRegions.size());
