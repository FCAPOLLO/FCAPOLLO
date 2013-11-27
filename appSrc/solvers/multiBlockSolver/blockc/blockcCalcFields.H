//	Set up necessary Gradients (These are not stored)

volVectorField phi_p_GRAD = fvc::grad(phi_p);

// Set the tensor for the proton conductivity
{
	tensor temp = tensor::zero;
	temp.component(tensor::XX) = sigma_p.value();
	temp.component(tensor::YY) = sigma_p.value();
	temp.component(tensor::ZZ) = sigma_p.value();

	dimensionedTensor sigmaTemp("sigmaTemp", dimensionSet(sigma_p.dimensions()), temp);
	sigma_p_eff = sigmaTemp;
}
