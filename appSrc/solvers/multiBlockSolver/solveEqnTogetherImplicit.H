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
    - solveEqnTogetherImplicit.H

Description
	- Solves the transport equations in the same sub-region loop and 
	uses implicit descriptions for the source terms
:
\*---------------------------------------------------------------------------*/
//Solve Left 
forAll(blocklRegions, zoneID)
{
	#include "blocklSetFields.H"
	source_electron_total = reaction_rate_BV_SUM;
	source_proton_total = reaction_rate_BV_SUM;
  	#include "electronTransportEquationA.H"
	#include "blocklCalcElectrochemistry.H"
	#include "protonTransportEquationA.H"
}
// Solve Centre 
forAll(blockcRegions, zoneID)
{
	#include "blockcSetFields.H"
	source_proton_total.internalField() = 0.0;
	source_proton_total.boundaryField() = 0.0;
	#include "protonTransportEquation.H"
}
// Solve Right 
forAll(blockrRegions, zoneID)
{
 	#include "blockrSetFields.H"
	source_proton_total = reaction_rate_BV_SUM;
	source_electron_total = reaction_rate_BV_SUM;
   	#include "electronTransportEquationC.H"
	#include "blockrCalcElectrochemistry.H"
	#include "protonTransportEquationC.H"
}
