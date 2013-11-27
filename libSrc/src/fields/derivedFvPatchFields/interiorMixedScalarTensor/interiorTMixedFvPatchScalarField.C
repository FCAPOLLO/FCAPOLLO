/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "interiorTMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace FCAPOLLO
{
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

interiorTMixedFvPatchScalarField::
interiorTMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    transportCoeffName_("undefined-neighbourFieldName")
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


interiorTMixedFvPatchScalarField::
interiorTMixedFvPatchScalarField
(
    const interiorTMixedFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(psf, p, iF, mapper),
    transportCoeffName_(psf.transportCoeffName_)
{}


interiorTMixedFvPatchScalarField::
interiorTMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    transportCoeffName_(dict.lookup("transportCoeff"))
{
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorIn
        (
            "interiorTMixedFvPatchScalarField::"
            "interiorTMixedFvPatchScalarField\n"
            "(\n"
            "    const fvPatch& p,\n"
            "    const DimensionedField<scalar, volMesh>& iF,\n"
            "    const dictionary& dict\n"
            ")\n"
        )   << "\n    patch type '" << p.type()
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << dimensionedInternalField().name()
            << " in file " << dimensionedInternalField().objectPath()
            << exit(FatalError);
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 1.0;
    }
}


interiorTMixedFvPatchScalarField::
interiorTMixedFvPatchScalarField
(
    const interiorTMixedFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(psf, iF),
    transportCoeffName_(psf.transportCoeffName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void interiorTMixedFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

	// Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

	const mappedPatchBase& mpp = 
		refCast<const mappedPatchBase>(patch().patch());
    const polyMesh& neighbourMesh = mpp.sampleMesh();
    const fvPatch& neighbourPatch = 
		refCast<const fvMesh>(neighbourMesh).boundary()[mpp.samplePolyPatch().index()];

	// Determine the unit normal vectors for the principal and neighbour patches
	// Note these should both be the same, connecting the cell centroids, 
	// but the sign should be different because of multi-block
	vectorField pUnitNormal(this->patch().Sf()/this->patch().magSf());
	vectorField nUnitNormal(neighbourPatch.Sf()/neighbourPatch.magSf());
	mpp.distribute(nUnitNormal);


	// Look up and assign internal Field values for the 
	// principal and neighbour patches
	const fvPatchScalarField& pFieldPhi = refCast
		<const fvPatchScalarField>
		(
			patch().lookupPatchField<volScalarField, scalar>
				(
				 	dimensionedInternalField().name()
				)
		);

	const fvPatchScalarField& nFieldPhi = refCast
		<const fvPatchScalarField>
		(
			neighbourPatch.lookupPatchField<volScalarField, scalar>
			(
			 	dimensionedInternalField().name()
			)
		);

	// Look up and assign transport Coefficient Field values for the 
	// principal and neighbour patches

	const fvPatchTensorField& pFieldK = refCast
		<const fvPatchTensorField>
		(
			patch().lookupPatchField<volTensorField, tensor>
			(
			 	transportCoeffName_
			)
		);

	const fvPatchTensorField& nFieldK = refCast
		<const fvPatchTensorField>
		(
			neighbourPatch.lookupPatchField<volTensorField, tensor>
			(
			 	transportCoeffName_
			)
		);

	// swap to obtain full local values
	// variable
	scalarField pIntFldPhi(pFieldPhi.patchInternalField());
	scalarField nIntFldPhi(nFieldPhi.patchInternalField());
	mpp.distribute(nIntFldPhi);
	
	// transport coefficient
	tensorField pIntTensorFldK(pFieldK.patchInternalField());
	tensorField nIntTensorFldK(nFieldK.patchInternalField());
	mpp.distribute(nIntTensorFldK);

	// Compute the transport coefficient in the normal direction
	// Note unit normals in multiBlock have a -ve sign on the exterior patches, so
	// the "interior" patch is band-aided by a double dot product
	scalarField pIntFldK(pIntTensorFldK & pUnitNormal & pUnitNormal);
	scalarField nIntFldK(nIntTensorFldK & nUnitNormal & nUnitNormal);

	// Inverse cell spacing
	scalarField pDelta(patch().deltaCoeffs());
	scalarField nDelta(neighbourPatch.deltaCoeffs());
	mpp.distribute(nDelta);

	// Apply cell spacing weighting on the conductivities
	scalarField pKDelta(pIntFldK*pDelta);
	scalarField nKDelta(nIntFldK*nDelta);

//	Info<< nl << "patchName: " << patch().name();
//	Info<< nl << "DimensionField " << dimensionedInternalField().name();
//	Info<< nl << "pDelta" << nl << pDelta << endl;
//	Info<< nl << "nDelta" << nl << nDelta << endl;
//	Info<< nl << "pKDelta" << nl << pKDelta << endl;
//	Info<< nl << "nKDelta" << nl << nKDelta << endl;
//	Info<< nl << "valueFraction" << nl << nKDelta/(nKDelta + pKDelta) << endl;

	// Determine the flux at the boundary patch based on a weighted
	// average of the principal and neighbour cells
	this->refValue() = nIntFldPhi;

	this->refGrad() = 0.;

	this->valueFraction() = (nKDelta) / (nKDelta + pKDelta);
	
    // Restore tag
    UPstream::msgType() = oldTag;
   
   	mixedFvPatchScalarField::updateCoeffs();
}


void interiorTMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
	os.writeKeyword("transportCoeff") << transportCoeffName_ 
		<< token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    interiorTMixedFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace FCAPOLLO
} // End namespace Foam


// ************************************************************************* //
