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
    - solveEqnSeparateImplicit.H

Description
	- Solves the transport equations in the separate sub-region loops and 
	uses implicit descriptions for the source terms
:
\*---------------------------------------------------------------------------*/
// Solve Electrons
for(scalar electronCount=0;electronCount<int(electronSweep);electronCount++)
{
	// Solve Left 
	forAll(blocklRegions, zoneID)
	{
   		#include "blocklSetFields.H"
   		#include "blocklReadSIMPLEControls.H"
//		source_electron_total.storePrevIter();
		if(THV>0.)
		{
			source_electron_total = reaction_rate_THV_SUM;
		}
		else
		{
			source_electron_total = reaction_rate_BV_SUM;
		}				
//		source_electron_total.relax();
  		#include "electronTransportEquationA.H"
	}
	// Solve Right 
	forAll(blockrRegions, zoneID)
	{
   		#include "blockrSetFields.H"
		#include "blockrReadSIMPLEControls.H"

//		source_electron_total.storePrevIter();
		source_electron_total = reaction_rate_BV_SUM;
//		source_electron_total.relax();
  	 	#include "electronTransportEquationC.H"
	}

	// Update sources if not the first time through the loop
	if(electronCount!=0)
	{
		forAll(blocklRegions, zoneID)
		{
   			#include "blocklSetFields.H"
			#include "blocklCalcElectrochemistry.H"
		}
		// Solve Right 
		forAll(blockrRegions, zoneID)
		{
		 	#include "blockrSetFields.H"
			#include "blockrCalcElectrochemistry.H"
		}
	}
}
			
// Correct Electrochemistry for the electron step	
// Solve Left 
forAll(blocklRegions, zoneID)
{
   	#include "blocklSetFields.H"
	#include "blocklCalcElectrochemistry.H"
}
// Solve Right 
forAll(blockrRegions, zoneID)
{
 	#include "blockrSetFields.H"
	#include "blockrCalcElectrochemistry.H"
}

// Solve Protons
for(scalar protonCount=0;protonCount<int(protonSweep);protonCount++)
{
	// Solve Left 
	forAll(blocklRegions, zoneID)
	{
    	#include "blocklSetFields.H"
	 	#include "blocklReadSIMPLEControls.H"
	//	source_proton_total.storePrevIter();
		if(THV>0.)
		{
			source_proton_total = reaction_rate_THV_SUM;
		}
		else
		{
			source_proton_total = reaction_rate_BV_SUM;
		}

//		source_proton_total.relax();
		#include "protonTransportEquationA.H"
	}
	// Solve Centre 
	forAll(blockcRegions, zoneID)
	{
   		#include "blockcSetFields.H"
		#include "blockcReadSIMPLEControls.H"
		source_proton_total.internalField() = 0.0;
		source_proton_total.boundaryField() = 0.0;
		#include "protonTransportEquation.H"
	}
	// Solve Right 
	forAll(blockrRegions, zoneID)
	{
   		#include "blockrSetFields.H"
		#include "blockrReadSIMPLEControls.H"

	//	source_proton_total.storePrevIter();
		source_proton_total = reaction_rate_BV_SUM;
	//	source_proton_total.relax();

		#include "protonTransportEquationC.H"
	}
	// Update sources if not the first time through the loop
	if(protonCount!=0)
	{
		forAll(blocklRegions, zoneID)
		{
   			#include "blocklSetFields.H"
			#include "blocklCalcElectrochemistry.H"
		}
		// Solve Right 
		forAll(blockrRegions, zoneID)
		{
		 	#include "blockrSetFields.H"
			#include "blockrCalcElectrochemistry.H"
		}
	}
}
/*
forAll(blocklRegions, zoneID)
{
	#include "blocklSetFields.H"
	phi_p.correctBoundaryConditions();
	label patchIDLeft = mesh.boundaryMesh().findPatchID("wall_l");
	label patchIDRight = mesh.boundaryMesh().findPatchID("block_l_to_block_c");
	const polyPatch& cPatchLeft = mesh.boundaryMesh()[patchIDLeft];
	const polyPatch& cPatchRight = mesh.boundaryMesh()[patchIDRight];
	const surfaceScalarField& magSf = mesh.magSf();

	scalar faceAreaLeft = gSum(magSf.boundaryField()[patchIDLeft]);
	scalar faceAreaRight = gSum(magSf.boundaryField()[patchIDRight]);
	scalarField fluxLeft = protonFlux.boundaryField()[patchIDLeft];
	scalarField fluxRight = protonFlux.boundaryField()[patchIDRight];
	
	scalar interfaceCurrentLeft = gSum(fluxLeft);
	scalar interfaceCurrentRight = gSum(fluxRight);

	scalar interfaceCurrentDensityLeft = interfaceCurrentLeft/(faceAreaLeft*10.);
	scalar interfaceCurrentDensityRight = interfaceCurrentRight/(faceAreaRight*10.);

	Info<< nl 
		<< "Flux Left:"
		<< nl
		<< "Left = " << fluxLeft 
	   	<< tab
		<< "Right = " << fluxRight
		<< tab
		<< "CD Left = " << interfaceCurrentDensityLeft 
		<< tab
		<< "CD Right = " << interfaceCurrentDensityRight
		<< endl;

	scalarField totalCurrentTest = source_electron_total*mesh.V();
	scalar sumCurrentTest = gSum(totalCurrentTest);
	scalar currentDensity = sumCurrentTest/(faceAreaLeft*10.);

	Info<< nl
		<< "CURRENT TEST"
		<< nl
		<< "totalCurrentSum = " << sumCurrentTest
		<< tab
		<< "currentDensity = " << currentDensity
		<< endl;	
}

forAll(blockcRegions, zoneID)
{
	#include "blockcSetFields.H"
	phi_p.correctBoundaryConditions();
	label patchIDLeft = mesh.boundaryMesh().findPatchID("block_c_to_block_l");
	label patchIDRight = mesh.boundaryMesh().findPatchID("block_c_to_block_r");
	const polyPatch& cPatchLeft = mesh.boundaryMesh()[patchIDLeft];
	const polyPatch& cPatchRight = mesh.boundaryMesh()[patchIDRight];
	const surfaceScalarField& magSf = mesh.magSf();

	scalar faceAreaLeft = gSum(magSf.boundaryField()[patchIDLeft]);
	scalar faceAreaRight = gSum(magSf.boundaryField()[patchIDRight]);
	scalarField fluxLeft = protonFlux.boundaryField()[patchIDLeft];
	scalarField fluxRight = protonFlux.boundaryField()[patchIDRight];
	
	scalar interfaceCurrentLeft = gSum(fluxLeft);
	scalar interfaceCurrentRight = gSum(fluxRight);

	scalar interfaceCurrentDensityLeft = interfaceCurrentLeft/(faceAreaLeft*10.);
	scalar interfaceCurrentDensityRight = interfaceCurrentRight/(faceAreaRight*10.);

	Info<< nl 
		<< "Flux Centre:"
		<< nl
		<< "Left = " << fluxLeft 
	   	<< tab
		<< "Right = " << fluxRight
		<< tab
		<< "CD Left = " << interfaceCurrentDensityLeft 
		<< tab
		<< "CD Right = " << interfaceCurrentDensityRight
		<< endl;
}

forAll(blockrRegions, zoneID)
{
	#include "blockrSetFields.H"
	phi_p.correctBoundaryConditions();
	label patchIDLeft = mesh.boundaryMesh().findPatchID("block_r_to_block_c");
	label patchIDRight = mesh.boundaryMesh().findPatchID("wall_r");
	const polyPatch& cPatchLeft = mesh.boundaryMesh()[patchIDLeft];
	const polyPatch& cPatchRight = mesh.boundaryMesh()[patchIDRight];
	const surfaceScalarField& magSf = mesh.magSf();

	scalar faceAreaLeft = gSum(magSf.boundaryField()[patchIDLeft]);
	scalar faceAreaRight = gSum(magSf.boundaryField()[patchIDRight]);
	scalarField fluxLeft = protonFlux.boundaryField()[patchIDLeft];
	scalarField fluxRight = protonFlux.boundaryField()[patchIDRight];
	
	scalar interfaceCurrentLeft = gSum(fluxLeft);
	scalar interfaceCurrentRight = gSum(fluxRight);

	scalar interfaceCurrentDensityLeft = interfaceCurrentLeft/(faceAreaLeft*10.);
	scalar interfaceCurrentDensityRight = interfaceCurrentRight/(faceAreaRight*10.);

	Info<< nl 
		<< "Flux Right:"
		<< nl
		<< "Left = " << fluxLeft 
	   	<< tab
		<< "Right = " << fluxRight
		<< tab
		<< "CD Left = " << interfaceCurrentDensityLeft 
		<< tab
		<< "CD Right = " << interfaceCurrentDensityRight
		<< endl;

	scalarField totalCurrentTest = source_electron_total*mesh.V();
	scalar sumCurrentTest = gSum(totalCurrentTest);
	scalar currentDensity = sumCurrentTest/(faceAreaLeft*10.);

	Info<< nl
		<< "CURRENT TEST"
		<< nl
		<< "totalCurrentSum = " << sumCurrentTest
		<< tab
		<< "currentDensity = " << currentDensity
		<< endl;	
}
*/
