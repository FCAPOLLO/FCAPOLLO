L2NormElectron = convergenceCheck(phi_e, phi_e_old, L2NormElectron);
L2NormProton = convergenceCheck(phi_p, phi_p_old, L2NormProton);

contErrElectron = contErrElectron + gSum(contErr_term1);
contErrProton = contErrProton + gSum(contErr_term2);
