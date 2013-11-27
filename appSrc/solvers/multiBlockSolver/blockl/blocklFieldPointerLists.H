// Initilialise field pointer lists

// Calculated Field Variables
	PtrList<volScalarField> eta_act_blockl(blocklRegions.size());
	PtrList<surfaceScalarField> proton_flux_blockl(blocklRegions.size());
	PtrList<surfaceScalarField> electron_flux_blockl(blocklRegions.size());
// Solved Field Variables	
	// Electron Transport Fields
	PtrList<volScalarField> phi_e_blockl(blocklRegions.size());
	// Proton Transport Fields
	PtrList<volScalarField> phi_p_blockl(blocklRegions.size());

// Coefficient Field Variables
	// Electron Transport Fields
	PtrList<volTensorField> sigma_e_eff_blockl(blocklRegions.size());
	// Proton Transport Fields
	PtrList<volTensorField> sigma_p_eff_blockl(blocklRegions.size());

// Source/sink Field Variables
	// Electron Transport Fields
	PtrList<volScalarField> reaction_rate_THV_ANODIC_blockl(blocklRegions.size());
	PtrList<volScalarField> reaction_rate_THV_CATHODIC_blockl(blocklRegions.size());
	PtrList<volScalarField> reaction_rate_THV_SUM_blockl(blocklRegions.size());
	PtrList<volScalarField> reaction_rate_BV_ANODIC_blockl(blocklRegions.size());
	PtrList<volScalarField> reaction_rate_BV_CATHODIC_blockl(blocklRegions.size());
	PtrList<volScalarField> reaction_rate_BV_SUM_blockl(blocklRegions.size());

	PtrList<volScalarField> reaction_rate_THV_SLOPE_ELECTRON_blockl(blocklRegions.size());
	PtrList<volScalarField> reaction_rate_THV_SLOPE_PROTON_blockl(blocklRegions.size());
	PtrList<volScalarField> reaction_rate_BV_SLOPE_ELECTRON_blockl(blocklRegions.size());
	PtrList<volScalarField> reaction_rate_BV_SLOPE_PROTON_blockl(blocklRegions.size());
	
	PtrList<volScalarField> source_electron_total_blockl(blocklRegions.size());
	PtrList<volScalarField> source_proton_total_blockl(blocklRegions.size());
	PtrList<volScalarField> source_electron_SLOPE_blockl(blocklRegions.size());
	PtrList<volScalarField> source_proton_SLOPE_blockl(blocklRegions.size());
