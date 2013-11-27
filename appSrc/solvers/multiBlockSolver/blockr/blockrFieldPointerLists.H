// Initilialise field pointer lists

// Calculated Field Variables
	PtrList<volScalarField> eta_act_blockr(blockrRegions.size());
	PtrList<surfaceScalarField> proton_flux_blockr(blockrRegions.size());
	PtrList<surfaceScalarField> electron_flux_blockr(blockrRegions.size());

// Solved Field Variables	
	// Electron Transport Fields
	PtrList<volScalarField> phi_e_blockr(blockrRegions.size());
	// Proton Transport Fields
	PtrList<volScalarField> phi_p_blockr(blockrRegions.size());
// Coefficient Field Variables
	// Electron Transport Fields
	PtrList<volTensorField> sigma_e_eff_blockr(blockrRegions.size());
	// Proton Transport Fields
	PtrList<volTensorField> sigma_p_eff_blockr(blockrRegions.size());

// Source/sink Field Variables
	// Electron Transport Fields
	PtrList<volScalarField> reaction_rate_BV_ANODIC_blockr(blockrRegions.size());
	PtrList<volScalarField> reaction_rate_BV_CATHODIC_blockr(blockrRegions.size());
	PtrList<volScalarField> reaction_rate_BV_SUM_blockr(blockrRegions.size());

	PtrList<volScalarField> reaction_rate_BV_SLOPE_ELECTRON_blockr(blockrRegions.size());
	PtrList<volScalarField> reaction_rate_BV_SLOPE_PROTON_blockr(blockrRegions.size());

	PtrList<volScalarField> source_electron_total_blockr(blockrRegions.size());
	PtrList<volScalarField> source_proton_total_blockr(blockrRegions.size());

	PtrList<volScalarField> source_electron_SLOPE_blockr(blockrRegions.size());
	PtrList<volScalarField> source_proton_SLOPE_blockr(blockrRegions.size());
