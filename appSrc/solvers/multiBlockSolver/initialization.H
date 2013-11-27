	// Determine the geometric cross-sectional Area
	scalar activeArea = 0.;
	scalar electrodeVol = 0.;
	forAll(blocklRegions, zoneID)
	{
		#include "blocklSetFields.H"
		
		label patchID = mesh.boundaryMesh().findPatchID("block_l_to_block_c");

		Info << "patchID = " << patchID << endl;

		const polyPatch& cPatch = mesh.boundaryMesh()[patchID];

		Info << "cPATCH = " << nl << cPatch << endl;

		const surfaceScalarField& magSf = mesh.magSf();

		scalarField volume = mesh.V();

		activeArea = gSum(magSf.boundaryField()[patchID]);
		electrodeVol = gSum(volume);
	}
	
	Info << "Active Area = " << nl << activeArea << "  [m^2]" << endl;
	Info << "Anode Vol = " << electrodeVol << "  [m^3]" << endl;

	// Initialise the Solution Set and Calculated Variables
	// Vars for initilization
	{
		// Setup calculated Fields, i.e. conductivities etc.
		forAll(blocklRegions, zoneID)
		{
			#include "blocklSetFields.H"
			// Set the tensor for the proton conductivity
			{
				tensor temp = tensor::zero;
				temp.component(tensor::XX) = .01;
				temp.component(tensor::YY) = .01;
				temp.component(tensor::ZZ) = .01;

				dimensionedTensor sigmaTemp("sigmaTemp", dimensionSet(sigma_p.dimensions()), temp);
				sigma_p_eff = sigmaTemp;
			}
		}
		forAll(blockcRegions, zoneID)
		{
			#include "blockcSetFields.H"
			// Set the tensor for the proton conductivity
			{
				tensor temp = tensor::zero;
				temp.component(tensor::XX) = 12.;
				temp.component(tensor::YY) = 12.;
				temp.component(tensor::ZZ) = 12.;

				dimensionedTensor sigmaTemp("sigmaTemp", dimensionSet(sigma_p.dimensions()), temp);
				sigma_p_eff = sigmaTemp;
			}
		}
		forAll(blockrRegions, zoneID)
		{
			#include "blockrSetFields.H"
			// Set the tensor for the proton conductivity
			{
				tensor temp = tensor::zero;
				temp.component(tensor::XX) = 0.5;
				temp.component(tensor::YY) = 0.5;
				temp.component(tensor::ZZ) = 0.5;

				dimensionedTensor sigmaTemp("sigmaTemp", dimensionSet(sigma_p.dimensions()), temp);
				sigma_p_eff = sigmaTemp;
			}
		}
	
		// Solution Fields
		Info<< nl
			<< "Initialisation of Solution Fields"
			<< endl;

		for(scalar i=0;i<100;i++)
		{
			forAll(blocklRegions, zoneID)
			{
				#include "blocklSetFields.H"
				phi_e.correctBoundaryConditions();
				phi_p.correctBoundaryConditions();
				scalar eta = eta_blockl_guess;
				phi_p.internalField() = phi_e.internalField() - eta;
				phi_p.boundaryField() = phi_p.boundaryField() - eta;
				phi_e.correctBoundaryConditions();
				phi_p.correctBoundaryConditions();

			}
			forAll(blockcRegions, zoneID)
			{
				#include "blockcSetFields.H"
				fvScalarMatrix protonic_transport
				(
 				//   fvm::ddt(phi_p)
    				-fvm::laplacian(sigma_p_eff,phi_p)
					== source_proton_total 
				);
				protonic_transport.solve();
				phi_p.correctBoundaryConditions();
			}
			forAll(blockrRegions, zoneID)
			{
				#include "blockrSetFields.H"
				scalar eta = eta_blockr_guess;
				phi_p.internalField() = phi_e.internalField() - eta - 1.18;
				phi_p.boundaryField() = phi_e.boundaryField() - eta - 1.18;
				phi_e.correctBoundaryConditions();
				phi_p.correctBoundaryConditions();
			}
		}

		forAll(blocklRegions, zoneID)
		{
			#include "blocklSetFields.H"
			phi_e.write();
			phi_p.write();
		}
		forAll(blockcRegions, zoneID)
		{
			#include "blockcSetFields.H"
			phi_p.write();
		}
		forAll(blockrRegions, zoneID)
		{
			#include "blockrSetFields.H"
			phi_e.write();
			phi_p.write();
		}
	
		// Calculated Fields
		Info<< nl
			<< "Initialisation of Calculated Fields"
			<< endl;

		forAll(blocklRegions, zoneID)
		{
			#include "blocklSetFields.H"
			#include "blocklCalcFields.H"
			#include "blocklCalcElectrochemistry.H"
		
		}
		forAll(blockcRegions, zoneID)
		{
			#include "blockcSetFields.H"
			#include "blockcCalcFields.H"
		}
		forAll(blockrRegions, zoneID)
		{
			#include "blockrSetFields.H"
			#include "blockrCalcFields.H"
			#include "blockrCalcElectrochemistry.H"
		}

	} // end initialization
