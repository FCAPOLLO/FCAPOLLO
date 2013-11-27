/*---------------------------------------------------------------------------*\
=========                 |
\\      /  F ield         | Unsupported Contributions for OpenFOAM
 \\    /   O peration     |
  \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation	
   \\/     M anipulation  |
-------------------------------------------------------------------------------
2013-11-13 David B. Harvey:

Notes:
	- Re-wrote turbulentTemperatureCoupleBaffleMixed to apply an interior
	patch coupling for a given scalar field and tensor transport property
	- Re-wrote modified boundary condition to work in parallel
	- Moved modified boundary condition to operate in FCAPOLLO namespace
-------------------------------------------------------------------------------
License
    FC-APOLLO and this file are a derivative work of OpenFOAM.

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

#include "interiorTFvPatchScalarField.H"
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

interiorTFvPatchScalarField::
interiorTFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    transportCoeffName_("undefined-transportCoeff")

{}

interiorTFvPatchScalarField::
interiorTFvPatchScalarField
(
    const interiorTFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    transportCoeffName_(ptf.transportCoeffName_)
{}

interiorTFvPatchScalarField::
interiorTFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    transportCoeffName_(dict.lookup("transportCoeff"))

{}

interiorTFvPatchScalarField::
interiorTFvPatchScalarField
(
    const interiorTFvPatchScalarField& ptf 
)
:
    fixedValueFvPatchScalarField(ptf)
{}

interiorTFvPatchScalarField::
interiorTFvPatchScalarField
(
    const interiorTFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    transportCoeffName_(ptf.transportCoeffName_)

{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void interiorTFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
}


void interiorTFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const interiorTFvPatchScalarField& hfptf =
        refCast<const interiorTFvPatchScalarField>(ptf);

}


void interiorTFvPatchScalarField::updateCoeffs()
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

	// Determine the value of phi at the boundary patch based on a weighted average of the principal and neighbour cells
    operator==( (pKDelta*pIntFldPhi + nKDelta*nIntFldPhi)/((pKDelta + nKDelta)) );


	// Restore tag
    UPstream::msgType() = oldTag;

    fixedValueFvPatchScalarField::updateCoeffs();

}

void interiorTFvPatchScalarField::write
(
    Ostream& os
) const
{
    fixedValueFvPatchScalarField::write(os);
	os.writeKeyword("transportCoeff") << transportCoeffName_ 
		<< token::END_STATEMENT << nl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
	fvPatchScalarField,
	interiorTFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace FCAPOLLO
} // End namespace Foam


// ************************************************************************* //
