/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


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
	// Note these should both be the same, connecting the cell centroids, but the sign should be different
	vectorField pUnitNormal(this->patch().Sf()/this->patch().magSf());
	vectorField nUnitNormal(neighbourPatch.Sf()/neighbourPatch.magSf());
	mpp.distribute(nUnitNormal);

	Info << "pUnitNormal" << nl << pUnitNormal << endl;
	Info << "nUnitNormal" << nl << nUnitNormal << endl;

	// note that the vector field for the neighbour unit normal errors when
	// mpp.distribute(nUnitNormal) is called

	// Look up and assign internal Field values and transport coefficients for the principal and neighbour patches
	const fvPatchScalarField& pFieldPhi = patch().lookupPatchField<volScalarField, scalar>(dimensionedInternalField().name());
	const fvPatchTensorField& pFieldK = patch().lookupPatchField<volTensorField, tensor>(transportCoeffName_);
    const fvPatchScalarField& nFieldPhi = neighbourPatch.lookupPatchField<volScalarField, scalar>(dimensionedInternalField().name());
	const fvPatchTensorField& nFieldK = neighbourPatch.lookupPatchField<volTensorField, tensor>(transportCoeffName_);

    // Swap to obtain full local values of neighbour internal field and compute the transport coefficient in the unit normal direction
	// Note: There is an issue with the unit Normals in multi-block, because it is an exterior patch the unit normal may be
	// negative i.e. pointing the wrong way for a interior style coupling.  The double dot product band-aids this problem.

    scalarField nIntFldPhi(nFieldPhi.patchInternalField());
	mpp.distribute(nIntFldPhi);
	scalarField nIntFldK((nFieldK.patchInternalField() & nUnitNormal) & nUnitNormal);
	mpp.distribute(nIntFldK);

	// Swap to obtain full local values of principal internal field and compute the transport coefficient in the unit normal direction
	// Note: There is an issue with the unit Normals in multi-block, because it is an exterior patch the unit normal may be
	// negative i.e. pointing the wrong way for a interior style coupling.  The double dot product band-aids this problem.
	const scalarField pIntFldPhi(pFieldPhi.patchInternalField());
	const scalarField pIntFldK((pFieldK.patchInternalField() & pUnitNormal) & pUnitNormal);

	// Determine the product of the transport coefficient and the inverse cell spacing
	scalarField pDelta(patch().deltaCoeffs());
	scalarField nDelta(neighbourPatch.deltaCoeffs());
	mpp.distribute(nDelta);

	scalarField pKDelta(pIntFldK*patch().deltaCoeffs());
	scalarField nKDelta(nIntFldK*neighbourPatch.deltaCoeffs());
	mpp.distribute(nKDelta);

	// Determine the value of phi at the boundary patch based on a weighted average of the principal and neighbour cells
	this->refValue() = nIntFldPhi;

	this->refGrad() = 0.0;

	this->valueFraction() = nKDelta / (nKDelta + pKDelta);

//		Info<< "Working Block" << nl << patch().boundaryMesh().mesh().name() << nl
//			<< "Name of Patch" << nl << patch().name() << nl
//			<< "Field Variable" << nl << this->dimensionedInternalField().name() << nl
//			<< neighbourMesh.name() << ':'
//			<< neighbourPatch.name() << ':'
//			<< this->dimensionedInternalField().name() << " :" << nl
//			<< "Patch Index" << nl << patch().index() << nl
//			<< "Neighbour Index" << nl << neighbourPatch.index()
		//	<< "Coeff Name" << nl << transportCoeffName_ << nl
//			<< "pFieldK" << nl << pFieldK.patchInternalField() << nl
//			<< "pIntFldk" << nl << pIntFldK << nl
//			<< "nFieldK" << nl << nFieldK.patchInternalField() << nl
//			<< "nIntFldK" << nl << nIntFldK << nl
//			<< "pUnitNormal" << nl << pUnitNormal << nl
//			<< "nUnitNormal" << nl << nUnitNormal << nl
		//	<< "INT location" << nl << patch().boundaryMesh().mesh().C() << nl
		//	<< "NBR location" << nl << neighbourPatch.boundaryMesh().mesh().C() << nl
		//	<< "INT DeltaCoeff" << nl << pDelta << nl
		//	<< "NBR DeltaCeoff" << nl << nDelta << nl
//			<< "refValue" << nl << refValue() << nl
//			<< "refGrad" << nl << refGrad() << nl
//			<< "valueFraction" << nl << valueFraction() << nl
//			<< endl;

/*	if (debug)
	{
		Info<< patch().boundaryMesh().mesh().name() << ':'
			<< patch().name() << ':'
			<< this->dimensionedInternalField().name() << " <- "
			<< neighbourMesh.name() << ':'
			<< neighbourPatch.name() << ':'
			<< this->dimensionedInternalField().name() << " :"
			<< "Coeff Name" << nl << transportCoeffName_
			<< "pFieldK" << nl << pFieldK.patchInternalField()
			<< "pIntFldk" << nl << pIntFldK
			<< "nFieldK" << nl << nFieldK.patchInternalField()
			<< "nIntFldK" << nl << nIntFldK
			<< "pUnitNormal" << nl << unitNormal
			<< "nUnitNormal" << nl << unitNormal
			<< "INT location" << nl << patch().boundaryMesh().mesh().C()
			<< "NBR location" << nl << neighbourPatch.boundaryMesh().mesh().C()
			<< "INT DeltaCoeff" << nl << patch().deltaCoeffs()
			<< "NBR DeltaCeoff" << nl << neighbourPatch.deltaCoeffs()
			<< endl;
	}
*/

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
} // End namespace Foam


// ************************************************************************* //
