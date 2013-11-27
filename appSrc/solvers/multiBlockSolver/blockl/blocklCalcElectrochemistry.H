// Overpotential
eta_act = (phi_e - phi_p);

// THV
reaction_rate_THV_ANODIC.internalField() = source_termTHV.value()*((4700.*(1.-Foam::exp((-2.*96485.*eta_act.internalField())/(1.2*8.314*353.))))+(100.*(Foam::exp((96485.*eta_act.internalField())/(2.*8.314*353.))-Foam::exp((-1.*96485.*eta_act.internalField())/(1.2*8.314*353.))*Foam::exp((-1.*96485.*eta_act.internalField())/(2.*8.314*353.)))));
reaction_rate_THV_ANODIC.boundaryField() = source_termTHV.value()*((4700.*(1.-Foam::exp((-2.*96485.*eta_act.boundaryField())/(1.2*8.314*353.))))+(100.*(Foam::exp((96485.*eta_act.boundaryField())/(2.*8.314*353.))-Foam::exp((-1.*96485.*eta_act.boundaryField())/(1.2*8.314*353.))*Foam::exp((-1.*96485.*eta_act.boundaryField())/(2.*8.314*353.)))));

reaction_rate_THV_CATHODIC.internalField() = source_termTHV.value()*((4700.*(Foam::exp((2.*96485.*eta_act.internalField())/(1.2*8.314*353.))-1.))+(100.*(-1*Foam::exp((-1.*96485.*eta_act.internalField())/(2.*8.314*353.))+Foam::exp((1.*96485.*eta_act.internalField())/(1.2*8.314*353.))*Foam::exp((1.*96485.*eta_act.internalField())/(2.*8.314*353.)))));
reaction_rate_THV_CATHODIC.boundaryField() = source_termTHV.value()*((4700.*(Foam::exp((2.*96485.*eta_act.boundaryField())/(1.2*8.314*353.))-1.))+(100.*(-1*Foam::exp((-1.*96485.*eta_act.boundaryField())/(2.*8.314*353.))+Foam::exp((1.*96485.*eta_act.boundaryField())/(1.2*8.314*353.))*Foam::exp((1.*96485.*eta_act.boundaryField())/(2.*8.314*353.)))));

// BV 
reaction_rate_BV_ANODIC.internalField() = source_termBV.value()*(Foam::exp(2.*0.5*eta_act.internalField()*96485./(8.314*353.)));
reaction_rate_BV_ANODIC.boundaryField() = source_termBV.value()*(Foam::exp(2.*0.5*eta_act.boundaryField()*96485./(8.314*353.)));

reaction_rate_BV_CATHODIC.internalField() = source_termBV.value()*(-Foam::exp(-2.*0.5*96485.*eta_act.internalField()/(8.314*353.)));
reaction_rate_BV_CATHODIC.boundaryField() = source_termBV.value()*(-Foam::exp(-2.*0.5*96485.*eta_act.boundaryField()/(8.314*353.)));

// Limiters
reaction_rate_THV_ANODIC = Foam::max(reaction_rate_THV_ANODIC, dimensionedScalar("reaction_rate_THV_CATHODIC_llimit", dimensionSet(0,-3, 0, 0, 0, 1, 0), 0.));
reaction_rate_THV_CATHODIC = Foam::min(reaction_rate_THV_CATHODIC, dimensionedScalar("reaction_rate_THV_ANODIC_ulimit", dimensionSet(0,-3, 0, 0, 0, 1, 0), 0.));

//reaction_rate_BV_ANODIC = Foam::max(reaction_rate_BV_ANODIC, dimensionedScalar("reaction_rate_BV_CATHODIC_llimit", dimensionSet(0,-3, 0, 0, 0, 1, 0), 0.000));
//reaction_rate_BV_CATHODIC = Foam::min(reaction_rate_BV_CATHODIC, dimensionedScalar("reaction_rate_BV_ANODIC_ulimit", dimensionSet(0,-3, 0, 0, 0, 1, 0), 0.000));

// Correct BC's
//reaction_rate_THV_ANODIC.correctBoundaryConditions();
//reaction_rate_THV_CATHODIC.correctBoundaryConditions();
//reaction_rate_BV_ANODIC.correctBoundaryConditions();
//reaction_rate_BV_CATHODIC.correctBoundaryConditions();

// Form the Sum
reaction_rate_THV_SUM = reaction_rate_THV_ANODIC + reaction_rate_THV_CATHODIC;
reaction_rate_BV_SUM = reaction_rate_BV_ANODIC + reaction_rate_BV_CATHODIC;

// Correct "Sum" BC's
//reaction_rate_THV_SUM.correctBoundaryConditions();
//reaction_rate_BV_SUM.correctBoundaryConditions();
