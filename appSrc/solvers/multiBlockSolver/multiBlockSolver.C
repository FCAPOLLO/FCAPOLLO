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
    multiBlockSolver

Description
	Steady State version of a CCM Electrochemical Model
:
\*---------------------------------------------------------------------------*/
#include <iostream>
#include <iomanip>
#include <unistd.h>
#include "fvCFD.H"
#include "regionProperties.H"
#include "interiorTMixedFvPatchScalarField.H"

int main(int argc, char *argv[])
{
	#include "setRootCase.H"
    #include "createTime.H"

    regionProperties rp(runTime);

    // Create Component Mesh and Field Objects
    Info<< "Creating component mesh objects at time = " << runTime.timeName() << nl << endl;
	
	// Create the Meshes
    #include "blocklCreateMesh.H"
	#include "blockcCreateMesh.H"
 	#include "blockrCreateMesh.H"
	
	// Create the fields
	#include "blocklCreateFields.H"
	#include "blockcCreateFields.H"
	#include "blockrCreateFields.H"

	// Setup the block identifiers
	#include "FCPEMDictionary.H"
	
	Info<<nl<< "Max Iteration set at " << maxSweep << endl;

	
	// Vars for looping
	scalar phi_p_left = 0;
	scalar phi_p_right = 0.;
	scalar counter = 0;
	scalar writeCounter = 0;
	scalar anodeCurrentDensity = 0.;
	scalar cathodeCurrentDensity = 0.;
	scalar currentDensityDiff = 1.e15;
	scalar stepTimeOld = 0.;
	scalar totalTime = 0.;


	// Initialize 
	#include "initialization.H"
	#include "convergenceOutput.H"

    while (runTime.run())
    {
	    runTime++;
	
		counter = 0;
		currentDensityDiff = 1.e15;

		Info<< nl 
		    << "Time = " << runTime.timeName() << endl;
		
		while(counter<int(maxSweep)&&(currentDensityDiff>currentTolerance))
		{

			// Calc the fields, this calcs the electrochemistry...
			forAll(blockrRegions, zoneID)
			{
				#include "blockrSetFields.H"
				#include "blockrCalcFields.H"
				#include "blockrCalcElectrochemistry.H"
			}
			forAll(blockcRegions, zoneID)
			{
   				#include "blockcSetFields.H"
				#include "blockcCalcFields.H"
			}

			forAll(blocklRegions, zoneID)
			{
  	 	 		#include "blocklSetFields.H"
				#include "blocklCalcFields.H"
				#include "blocklCalcElectrochemistry.H"
			}
	
			// Solve the equations
			// "implicit" adds implicit sources
			// "separate" solves the equations in individual loops
			if(apolloSolverMethod=="implicit")
			{
				if(apolloEqnSystem=="separate")
				{
					#include "solveEqnSeparateImplicit.H"
				}
				else
				{
					#include "solveEqnTogetherImplicit.H"
				}
			}
			else
			{
				if(apolloEqnSystem=="separate")
				{
					#include "solveEqnSeparate.H"
				}
				else
				{
					#include "solveEqnTogether.H"
				}
			}		
			
			#include "convergenceOutput.H"
		   	
		    counter++;
		}

		
		forAll(blockrRegions, zoneID)
		{
			#include "blockrSetFields.H"
			phi_e.correctBoundaryConditions();
			phi_p.correctBoundaryConditions();
		}
		forAll(blockcRegions, zoneID)
		{
   			#include "blockcSetFields.H"
			phi_p.correctBoundaryConditions();
		}
		forAll(blocklRegions, zoneID)
		{
  	  		#include "blocklSetFields.H"
			phi_e.correctBoundaryConditions();
			phi_p.correctBoundaryConditions();
		}
		
		Info<< endl;

		runTime.write();

		#include "convergenceOutput.H"

		Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << nl
 	       << "ClockTime = " << runTime.elapsedClockTime() << " s" << nl << endl; 
    }

 	Info<< "End\n" << endl;

 	return 0;
}


// ************************************************************************* //
