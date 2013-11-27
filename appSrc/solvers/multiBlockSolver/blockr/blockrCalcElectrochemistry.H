// Overpotential
dimensionedScalar OCV("OCV", dimensionSet( 1, 2, -3, 0, 0, -1, 0 ), 1.18);

eta_act = ((phi_e - phi_p) - OCV);

// Reaction Rates
reaction_rate_BV_ANODIC.internalField() = source_term.value()*(Foam::exp(2.*0.51*96485.*eta_act.internalField()/(8.314*353.)));
reaction_rate_BV_ANODIC.boundaryField() = source_term.value()*(Foam::exp(2.*0.51*96485.*eta_act.boundaryField()/(8.314*353.)));

reaction_rate_BV_CATHODIC.internalField() = source_term.value()*(-Foam::exp(-2.*0.49*96485.*eta_act.internalField()/(8.314*353.)));
reaction_rate_BV_CATHODIC.boundaryField() = source_term.value()*(-Foam::exp(-2.*0.49*96485.*eta_act.boundaryField()/(8.314*353.)));

// Form the Sum
reaction_rate_BV_SUM = reaction_rate_BV_ANODIC + reaction_rate_BV_CATHODIC;
