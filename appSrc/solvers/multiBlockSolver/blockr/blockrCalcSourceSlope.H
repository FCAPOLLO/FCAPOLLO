		reaction_rate_BV_SLOPE_ELECTRON = (reaction_rate_BV_SUM - reaction_rate_BV_SUM.prevIter())
									/(phi_e-phi_e.prevIter());
		reaction_rate_BV_SLOPE_PROTON = (reaction_rate_BV_SUM - reaction_rate_BV_SUM.prevIter())
									/(phi_p-phi_p.prevIter());

		source_electron_SLOPE = reaction_rate_BV_SLOPE_ELECTRON;
		source_proton_SLOPE = reaction_rate_BV_SLOPE_ELECTRON;

scalarField elecS_new_1 = reaction_rate_BV_SLOPE_ELECTRON*phi_e.prevIter();
scalarField elecS_new_2 = reaction_rate_BV_SLOPE_ELECTRON*phi_e;
scalarField protS_new_1 = reaction_rate_BV_SLOPE_PROTON*phi_p.prevIter();
scalarField protS_new_2 = reaction_rate_BV_SLOPE_PROTON*phi_p;

scalarField elecNewSource = source_electron_total - reaction_rate_BV_SLOPE_ELECTRON*phi_e.prevIter() + reaction_rate_BV_SLOPE_ELECTRON*phi_e;
scalarField protNewSource = source_proton_total - reaction_rate_BV_SLOPE_PROTON*phi_p.prevIter() + reaction_rate_BV_SLOPE_PROTON*phi_p;

/*
	Info<< "Block Right"
	<< nl
	<< "Electron Slope*phi_e.prevIter"
	<< nl
	<< elecS_new_1
	<< nl
	<< "Electron Slope*phi_e"
	<< nl
	<< elecS_new_2
	<< nl
	<< "Proton Slope*phi_p.prevIter"
	<< nl
	<< protS_new_1
	<< nl
	<< "Proton Slope*phi_p"
	<< nl
	<< protS_new_2
	<< nl
	<< "Source Electron"
	<< nl
	<< source_electron_total.internalField()
	<< nl
	<< "Source Electron Modified"
	<< nl
	<< elecNewSource
	<< nl
	<< "Source Proton"
	<< nl
	<< source_proton_total.internalField()
	<< nl
	<< "Source Proton Modified"
	<< nl
	<< protNewSource
	<< endl;
	*/
