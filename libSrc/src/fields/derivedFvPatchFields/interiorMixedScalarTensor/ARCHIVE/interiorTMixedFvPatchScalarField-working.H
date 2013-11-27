/*---------------------------------------------------------------------------*\
    F uel            | FC-APOLLO: 
	C ell		     | 	The Open-source Fuel Cell Application Package 
    A pplication     | 
	P ackage using   |
	O pen-source for |  Copyright (C) 2013  Ballard Power Systems
	L ong            |
	L ife            |
	O peration       |
-------------------------------------------------------------------------------
License
    This file is part of FC-APOLLO.

    FC-APOLLO is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    FC-APOLLO is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with FC-APOLLO.  If not, see <http://www.gnu.org/licenses/>.

Class
	interiorTMixedFvPatchScalarField

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
            type        		interiorMixedScalarTensor;
            transportCoeff		transportCoefficientName;
			value				initialization value for the face;
        }

    Needs to be on underlying mapped(Wall)FvPatch.

    Note: This does not currently run in parallel, code fails when interior patches
	are returned with zero faces (i.e. nFaces = 0).

SourceFiles
    interiorTMixedFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef interiorTMixedFvPatchScalarField_H
#define interiorTMixedFvPatchScalarField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
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
    TypeName("interiorMixedScalarTensor");

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

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
