/*---------------------------------------------------------------------------*\
=========                 |
\\      /  F ield         | Unsupported Contributions for OpenFOAM
 \\    /   O peration     |
  \\  /    A nd           | Copyright (C) 2013  Ballard Power Systems
   \\/     M anipulation  |
-------------------------------------------------------------------------------
    F uel           	  | FC-APOLLO: 
	C ell		          |	The Open-source Fuel Cell Application Package 
    A pplication          | 
	P ackage using        |
	O pen-source for      | Authors: David B. Harvey
	L ong                 |          Alexander Bellemare-Davis
	L ife                 |                     
	O peration            |
-------------------------------------------------------------------------------
License
    FC-APOLLO and this file are a derivative work of OpenFOAM.

	FCAPOLLO is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    FCAPOLLO is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with FCAPOLLO.  If not, see <http://www.gnu.org/licenses/>.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    - electronTransportEquationA.H

Description
	- Solves the electron transport equation with explicit/implicit terms 
	in the left block and calculates the slope using a numerical 
	estimated derivative dS_dphi. 
:
\*---------------------------------------------------------------------------*/
{
	//  Protonic Transport
	// Request previous iteration to be stored
	phi_e.storePrevIter();

	// Calculate the slope
	volScalarField dS_dphi(source_electron_total);
	dS_dphi.dimensions().reset(dimless);
	dS_dphi = 0.;
	dS_dphi.dimensions().reset(source_electron_total.dimensions()/phi_e.dimensions());

	// calculate the slope
	{
		// Setup clones of the needed volume fields
		scalar dPhi = SMALL;
		volScalarField etaActLeft_(eta_act);
		volScalarField etaActRight_(eta_act);
		volScalarField RRANODICLeft_(source_electron_total);
		volScalarField RRANODICRight_(source_electron_total);
		volScalarField RRCATHODICLeft_(source_electron_total);
		volScalarField RRCATHODICRight_(source_electron_total);
		volScalarField RRLeft_(source_electron_total);
		volScalarField RRRight_(source_electron_total);

		// Reset dims as required
		RRANODICLeft_.dimensions().reset(dimless);
		RRANODICRight_.dimensions().reset(dimless);
		RRCATHODICLeft_.dimensions().reset(dimless);
		RRCATHODICRight_.dimensions().reset(dimless);

		// Create the change in phi_e
		etaActLeft_.internalField() = 
			phi_e.internalField()-phi_e.internalField() -dPhi;
		etaActLeft_.boundaryField() = 
			phi_e.boundaryField() - phi_e.boundaryField() -dPhi;
		etaActRight_.internalField() = 
			phi_e.internalField() - phi_e.internalField() +dPhi;
		etaActRight_.boundaryField() = 
			phi_e.boundaryField() - phi_e.boundaryField() +dPhi;

		// Calculate the change in Reaction Rate for phi_e +/- SMALL
		RRANODICLeft_.internalField() = source_termBV.value()
			*(Foam::exp(2.*0.5*etaActLeft_.internalField()*96485./(8.314*353.)));
		RRANODICLeft_.boundaryField() =  source_termBV.value()
			*(Foam::exp(2.*0.5*etaActLeft_.boundaryField()*96485./(8.314*353.)));
		RRANODICRight_.internalField() = source_termBV.value()
			*(Foam::exp(2.*0.5*etaActRight_.internalField()*96485./(8.314*353.)));				
		RRANODICRight_.boundaryField() = source_termBV.value()
			*(Foam::exp(2.*0.5*etaActRight_.boundaryField()*96485./(8.314*353.)));

		RRCATHODICLeft_.internalField() = source_termBV.value()
			*(-Foam::exp(-2.*0.5*96485.*etaActLeft_.internalField()/(8.314*353.)));
		RRCATHODICLeft_.boundaryField() = source_termBV.value()
			*(-Foam::exp(-2.*0.5*96485.*etaActLeft_.boundaryField()/(8.314*353.)));
		RRCATHODICRight_.internalField() = source_termBV.value()
			*(-Foam::exp(-2.*0.5*96485.*etaActRight_.internalField()/(8.314*353.)));
		RRCATHODICRight_.boundaryField() = source_termBV.value()
			*(-Foam::exp(-2.*0.5*96485.*etaActRight_.boundaryField()/(8.314*353.)));

		// reset dims on req. variables
		RRANODICLeft_.dimensions().reset(source_electron_total.dimensions());
		RRCATHODICLeft_.dimensions().reset(source_electron_total.dimensions());
		RRANODICRight_.dimensions().reset(source_electron_total.dimensions());
		RRCATHODICRight_.dimensions().reset(source_electron_total.dimensions());

		// Set the value of the RRLeft_ and RRRight_
		RRLeft_ = RRANODICLeft_ + RRCATHODICLeft_;
		RRRight_ = RRANODICRight_ + RRCATHODICRight_;

		dS_dphi = (RRRight_-RRLeft_)/(etaActRight_-etaActLeft_);

	}
	
	// Enhance the implicit term for stability
	dS_dphi.dimensions().reset(dimless);
	dS_dphi = dS_dphi*implicitFactor;
	dS_dphi.dimensions().reset(source_electron_total.dimensions()/phi_e.dimensions());

//	Info<< "ds_dphiA_e = " << gAverage(dS_dphi) << endl;
	// End calculate the slope

	fvScalarMatrix electronic_transport
	(
	 //   fvm::ddt(phi_e)
	    -fvm::laplacian(sigma_e_eff,phi_e)
		== -source_electron_total -(fvm::SuSp(dS_dphi,phi_e) -dS_dphi*phi_e.prevIter())
	);

	electronic_transport.relax();
	electronFlux = electronic_transport.flux();
	electronic_transport.solve();
	phi_e.relax();
//	phi_e.correctBoundaryConditions();
}
