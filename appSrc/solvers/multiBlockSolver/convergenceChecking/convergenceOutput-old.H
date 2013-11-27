	forAll(blocklRegions, zoneID)
	{
		#include "blocklSetFields.H"
		scalar anodeCurrent = gSum(source_electron_total.internalField()*mesh.V());
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
		
		Info<< "************* Block Left *************" << nl << endl;
		Info<< "Iteration = " << counter << endl;
		Info<< "Total Current [mA/cm2]"            
			<< nl << anodeCurrentDensity << endl;
		Info<< "Reaction Rate (THV=1,BV=-1)"      
				<< nl << THV << endl;
		Info<< "Reaction Rate Current [mA/cm2]"
			<< nl << reactionRateCurrentDensity << endl;
		Info<< "MIN R.Rate Current [A/m3]" 
			<< nl << reactionRateMin << endl;
		Info<< "MAX R.Rate Current [A/m3]" 
			<< nl << reactionRateMax << endl;
		Info<< "AVE R.Rate Current [A/m3]"
			<< nl << reactionRateAve << endl;
		Info<< "AVE Eta"  
			<< nl << gAverage(eta_act) << endl;
		Info<< "MIN Eta"                           
			<< nl << gMin(eta_act) << endl;
		Info<< "MAX Eta"                           
			<< nl << gMax(eta_act) << endl;
		Info<< "AVE phi_e" 	                 	  
			<< nl << gAverage(phi_e) << endl;
		Info<< "MIN phi_e"                     
			<< nl << gMin(phi_e) << endl;
		Info<< "MAX phi_e"                    
			<< nl << gMax(phi_e) << endl;
		Info<< "AVE phi_p" 	                    
			<< nl << gAverage(phi_p) << endl;
		Info<< "MIN phi_p"                        
			<< nl << gMin(phi_p) << endl;
		Info<< "MAX phi_p"                   
			<< nl << gMax(phi_p) << endl;
	}
	forAll(blockcRegions, zoneID)
	{
		
		#include "blockcSetFields.H"
		Info<< "************* Block Centre *************" << endl;

		Info<< "AVE phi_p" 
			<< nl << gAverage(phi_p) << endl;
		Info<< "MIN phi_p" 
			<< nl << gMin(phi_p) << endl;
		Info<< "MAX phi_p"
			<< nl << gMax(phi_p) << endl;

	}
	forAll(blockrRegions, zoneID)
	{
		
		#include "blockrSetFields.H"
		scalar cathodeCurrent = gSum(source_proton_total.internalField()*mesh.V());
		cathodeCurrentDensity = (cathodeCurrent)/(activeArea*10.);
	
		scalar reactionRateCurrentDensity = gSum(reaction_rate_BV_SUM.internalField()
											*mesh.V())/(activeArea*10.);
	
		scalar reactionRateMax = gMax(reaction_rate_BV_SUM);
		scalar reactionRateMin = gMin(reaction_rate_BV_SUM);
		scalar reactionRateAve = gAverage(reaction_rate_BV_SUM);

		Info<< "************* Block Right *************" << endl;
		Info<< "Total Current [mA/cm2]" 
			<< nl << cathodeCurrentDensity << endl;
		Info<< "Reaction Rate Current [mA/cm2]"   
			<< nl << reactionRateCurrentDensity << endl;
		Info<< "MIN R.Rate Current [A/m3]"         
			<< nl << reactionRateMin << endl;
		Info<< "MAX R.Rate Current [A/m3]"        
			<< nl << reactionRateMax << endl;
		Info<< "AVE R.Rate Current [A/m3]"        
			<< nl << reactionRateAve << endl;
		Info<< "AVE Eta"             	        
			<< nl << gAverage(eta_act) << endl;
		Info<< "MIN Eta"                           
			<< nl << gMin(eta_act) << endl;
		Info<< "MAX Eta"                          
				<< nl << gMax(eta_act) << endl;
		Info<< "AVE phi_e" 	                 
				<< nl << gAverage(phi_e) << endl;
		Info<< "MIN phi_e"                       
			<< nl << gMin(phi_e) << endl;
		Info<< "MAX phi_e"                        
			<< nl << gMax(phi_e) << endl;
		Info<< "AVE phi_p" 	                    
			<< nl << gAverage(phi_p) << endl;
		Info<< "MIN phi_p"                        
			<< nl << gMin(phi_p) << endl;
		Info<< "MAX phi_p"                        
			<< nl << gMax(phi_p) << endl;
		
	}	
	
