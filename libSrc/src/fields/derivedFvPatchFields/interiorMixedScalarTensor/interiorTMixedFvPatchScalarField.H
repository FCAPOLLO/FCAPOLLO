/*---------------------------------------------------------------------------*\
=========                 |
\\      /  F ield         | Unsupported Contributions for OpenFOAM
 \\    /   O peration     |
  \\  /    A nd           | Copyright (C) 2013	OpenFOAM Foundation
   \\/     M anipulation  |
-------------------------------------------------------------------------------
2013-01-13 David B. Harvey:

Notes:
	- Re-wrote turbulentTemperatureCoupleBaffleMixed to apply an interior
	coupling for a given scalar field and tensor transport property
	- Re-wrote modified boundary condition to work in parallel
	- Moved modified boundary condition to operate in FCAPOLLO namespace
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

Class
	FCAPOLLO::interiorTMixedFvPatchScalarField

Description
    This is a mixed boundary condition intended for a scalar field with
	a tensor transport ceofficient.

	This follows the convention of the mixed boundary conditions as
	written for OpenFOAM's chtMultiRegionHeater tutorial

    Specifies gradient and potential to equate the flux on both sides
	- refGradient = zero gradient
    - refValue = neighbour value
    - mixFraction = nKDelta / (nKDelta + pKDelta)

    where pKDelta and nKDelta are the transport coefficient * inverse cell spacing

    Example usage:
        pInterfaceName_to_nInterfaceName
        {
            type        		FCAPOLLO::interiorMixedScalarTensor;
            transportCoeff		transportCoefficientName;
			value				initialization value for the face;
        }

    Needs to be on underlying mapped(Wall)FvPatch.

SourceFiles
    interiorTMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef interiorTMixedFvPatchScalarField_H
#define interiorTMixedFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace FCAPOLLO
{
/*---------------------------------------------------------------------------*\
   Class interiorTMixedFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class interiorTMixedFvPatchScalarField
:
public mixedFvPatchScalarField
{
    // Private data

        //- Name of field on the neighbour region
        const word transportCoeffName_;

public:

    //- Runtime type information
    TypeName("FCAPOLLO::interiorMixedScalarTensor");

    // Constructors

        //- Construct from patch and internal field
        interiorTMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        interiorTMixedFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  interiorTMixedFvPatchScalarField onto a
        //  new patch
        interiorTMixedFvPatchScalarField
        (
            const interiorTMixedFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new interiorTMixedFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        interiorTMixedFvPatchScalarField
        (
            const interiorTMixedFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new interiorTMixedFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }

    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace FCAPOLLO
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
