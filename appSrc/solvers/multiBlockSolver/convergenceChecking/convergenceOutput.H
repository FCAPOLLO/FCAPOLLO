Info<< nl << nl << nl << nl << nl << nl << nl << nl << endl;
forAll(blocklRegions, zoneID)
{
	#include "blocklSetFields.H"
	const scalar anodeCurrent = gSum(source_electron_total.internalField()*mesh.V());
	anodeCurrentDensity = (anodeCurrent)/(activeArea*10.);
	scalar reactionRateCurrentDensity = 0.;
	scalar reactionRateMax = 0.;
	scalar reactionRateMin = 0.;
	scalar reactionRateAve = 0.;

	if(THV>0.)
	{
		reactionRateCurrentDensity = gSum(reaction_rate_THV_SUM.internalField()
										*mesh.V())/(activeArea*10.);
		reactionRateMax = gMax(reaction_rate_THV_SUM);
		reactionRateMin = gMin(reaction_rate_THV_SUM);
		reactionRateAve = gAverage(reaction_rate_THV_SUM);
	}
	else
	{
		reactionRateCurrentDensity = gSum(reaction_rate_BV_SUM.internalField()
									*mesh.V())/(activeArea*10.);
		reactionRateMax = gMax(reaction_rate_BV_SUM);
		reactionRateMin = gMin(reaction_rate_BV_SUM);
		reactionRateAve = gAverage(reaction_rate_BV_SUM);
	}
	
	cout<< "\n"
		<< "******************************* Block Left *******************************"
		<< "\n"
		<< std::left
		<< std::setw(30) << "Total Curent Density" 
		<< "\n"	
		<< std::setw(30) << "Source Term" 
		<< std::setw(30) << "Reaction Rate HOR"
		<< "\n"
		<< std::left << std::setprecision(2) << std::fixed
		<< std::setw(30) << anodeCurrentDensity
		<< std::setw(30) << reactionRateCurrentDensity
		<< "\n"	
		<< resetiosflags(ios_base::floatfield)
		<< std::setprecision(6)
		<< std::setw(30) << "Max Reaction Rate HOR"
		<< std::setw(30) << "Min Reaction Rate HOR"
		<< std::setw(30) << "Ave Reaction Rate HOR"
		<< "\n"
		<< std::setw(30) << reactionRateMax
		<< std::setw(30) << reactionRateMin
		<< std::setw(30) << reactionRateAve
		<< "\n"
		<< std::setw(30) << "Average eta"
		<< std::setw(30) << "Max eta"
		<< std::setw(30) << "Min eta"
		<< "\n"
		<< std::setw(30) << gAverage(eta_act)
		<< std::setw(30) << gMax(eta_act)
		<< std::setw(30) << gMin(eta_act)
		<< "\n"
		<< std::setw(30) << "Average phi_e"
		<< std::setw(30) << "Max phi_e"
		<< std::setw(30) << "Min phi_e"
		<< "\n"
		<< std::setw(30) << gAverage(phi_e)
		<< std::setw(30) << gMax(phi_e)
		<< std::setw(30) << gMin(phi_e)
		<< "\n"
		<< std::setw(30) << "Average phi_p"
		<< std::setw(30) << "Max phi_p"
		<< std::setw(30) << "Min phi_p"
		<< "\n"
		<< std::setw(30) << gAverage(phi_p)
		<< std::setw(30) << gMax(phi_p)
		<< std::setw(30) << gMin(phi_p)
		<< "\n";

}
forAll(blockrRegions, zoneID)
{
	
	#include "blockrSetFields.H"
	const scalar cathodeCurrent = gSum(source_proton_total.internalField()*mesh.V());
	cathodeCurrentDensity = (cathodeCurrent)/(activeArea*10.);

	scalar reactionRateCurrentDensity = gSum(reaction_rate_BV_SUM.internalField()
										*mesh.V())/(activeArea*10.);
	scalar reactionRateMax = gMax(reaction_rate_BV_SUM);
	scalar reactionRateMin = gMin(reaction_rate_BV_SUM);
	scalar reactionRateAve = gAverage(reaction_rate_BV_SUM);

	cout<< "\n"
		<< "******************************* Block Right *******************************"
		<< "\n"
		<< std::left
		<< std::setw(30) << "Total Curent Density" 
		<< "\n"	
		<< std::setw(30) << "Source Term" 
		<< std::setw(30) << "Reaction Rate OOR"
		<< "\n"
		<< std::left << std::setprecision(2) << std::fixed
		<< std::setw(30) << cathodeCurrentDensity
		<< std::setw(30) << reactionRateCurrentDensity
		<< "\n"	
		<< resetiosflags(ios_base::floatfield)
		<< std::setprecision(6)
		<< std::setw(30) << "Max Reaction Rate ORR"
		<< std::setw(30) << "Min Reaction Rate ORR"
		<< std::setw(30) << "Ave Reaction Rate ORR"
		<< "\n"
		<< std::setw(30) << reactionRateMax
		<< std::setw(30) << reactionRateMin
		<< std::setw(30) << reactionRateAve
		<< "\n"
		<< std::setw(30) << "Average eta"
		<< std::setw(30) << "Max eta"
		<< std::setw(30) << "Min eta"
		<< "\n"
		<< std::setw(30) << gAverage(eta_act)
		<< std::setw(30) << gMax(eta_act)
		<< std::setw(30) << gMin(eta_act)
		<< "\n"
		<< std::setw(30) << "Average phi_e"
		<< std::setw(30) << "Max phi_e"
		<< std::setw(30) << "Min phi_e"
		<< "\n"
		<< std::setw(30) << gAverage(phi_e)
		<< std::setw(30) << gMax(phi_e)
		<< std::setw(30) << gMin(phi_e)
		<< "\n"
		<< std::setw(30) << "Average phi_p"
		<< std::setw(30) << "Max phi_p"
		<< std::setw(30) << "Min phi_p"
		<< "\n"
		<< std::setw(30) << gAverage(phi_p)
		<< std::setw(30) << gMax(phi_p)
		<< std::setw(30) << gMin(phi_p)
		<< "\n";

}	
// Check the current imbalance
{
	if(counter>1)
	{
		currentDensityDiff = 
				fabs(fabs(anodeCurrentDensity) - fabs(cathodeCurrentDensity));	
	}

	cout<< "\n" 
		<< "******************************* Iteration Summary *******************************"
		<< "\n"
		<< std::left
		<< std::setw(30) << "anodeCurrentDensity" 
		<< std::setw(30) << "cathodeCurrentDensity"
		<< "\n" 
		<< std::setw(30) << anodeCurrentDensity
		<< std::setw(30) << cathodeCurrentDensity
		<< "\n";
	const scalar stepTime = runTime.elapsedCpuTime();
	const scalar deltaTime = stepTime - stepTimeOld;
	stepTimeOld = stepTime;
	totalTime += deltaTime;
	
	cout<< "\n"
   		<< std::left
		<< std::setw(12) << " "
		<< std::setw(16) << "Convergence Checking"
		<< "\n"
		<< std::setw(12) << "Iteration" 
		<< std::setw(16) << "Current Consv."
		<< std::setw(16) << "Iteration Time" 
		<< std::setw(16) << "Total Time"
		<< "\n"
		<< std::setw(12) << counter
		<< std::setw(16) << std::setprecision(4) << currentDensityDiff
		<< std::setw(16) << currentDensityDiff
		<< std::setw(16) << deltaTime
		<< std::setw(16) << totalTime
		<< "\n";
}
cout<< resetiosflags(ios_base::floatfield);
