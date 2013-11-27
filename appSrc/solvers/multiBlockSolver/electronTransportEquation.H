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
    - electronTransportEquation.H

Description
	- Solves the electron transport equation with explicit source terms  
:
\*---------------------------------------------------------------------------*/
{
	//  Electronic Transport
	// Request previous iteration to be stored
	phi_e.storePrevIter();

	fvScalarMatrix electronic_transport
	(
	  	// fvm::ddt(phi_e)
   		-fvm::laplacian(sigma_e_eff,phi_e) 
		== -source_electron_total
	);

	electronic_transport.relax();
	electronFlux = electronic_transport.flux();
	electronic_transport.solve();
	phi_e.relax();
//	phi_e.correctBoundaryConditions();
}
