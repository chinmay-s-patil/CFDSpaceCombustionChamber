/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
    Copyright (C) YEAR AUTHOR, AFFILIATION
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

#include "fixedValueFvPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
#include "PatchFunction1.H"

//{{{ begin codeInclude

//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// dynamicCode:
// SHA1 = 20ec42ca85ea0edcbbe96f36ed516a1ad0b6ed99
//
// unique function name that can be checked if the correct library version
// has been loaded
extern "C" void chamberUpT_20ec42ca85ea0edcbbe96f36ed516a1ad0b6ed99(bool load)
{
    if (load)
    {
        // Code that can be explicitly executed after loading
    }
    else
    {
        // Code that can be explicitly executed before unloading
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatchScalarField,
    chamberUpTFixedValueFvPatchScalarField
);

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
chamberUpTFixedValueFvPatchScalarField::
chamberUpTFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    parent_bctype(p, iF)
{
    if (false)
    {
        printMessage("Construct chamberUpT : patch/DimensionedField");
    }
}


Foam::
chamberUpTFixedValueFvPatchScalarField::
chamberUpTFixedValueFvPatchScalarField
(
    const chamberUpTFixedValueFvPatchScalarField& rhs,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    parent_bctype(rhs, p, iF, mapper)
{
    if (false)
    {
        printMessage("Construct chamberUpT : patch/DimensionedField/mapper");
    }
}


Foam::
chamberUpTFixedValueFvPatchScalarField::
chamberUpTFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    parent_bctype(p, iF, dict)
{
    if (false)
    {
        printMessage("Construct chamberUpT : patch/dictionary");
    }
}


Foam::
chamberUpTFixedValueFvPatchScalarField::
chamberUpTFixedValueFvPatchScalarField
(
    const chamberUpTFixedValueFvPatchScalarField& rhs
)
:
    parent_bctype(rhs),
    dictionaryContent(rhs)
{
    if (false)
    {
        printMessage("Copy construct chamberUpT");
    }
}


Foam::
chamberUpTFixedValueFvPatchScalarField::
chamberUpTFixedValueFvPatchScalarField
(
    const chamberUpTFixedValueFvPatchScalarField& rhs,
    const DimensionedField<scalar, volMesh>& iF
)
:
    parent_bctype(rhs, iF)
{
    if (false)
    {
        printMessage("Construct chamberUpT : copy/DimensionedField");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::
chamberUpTFixedValueFvPatchScalarField::
~chamberUpTFixedValueFvPatchScalarField()
{
    if (false)
    {
        printMessage("Destroy chamberUpT");
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::
chamberUpTFixedValueFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        printMessage("updateCoeffs chamberUpT");
    }

//{{{ begin code
    #line 45 "/home/lavender/OpenFoamUbu/CFDSpaceCombustionChamber/openFOAMTry/simFile/try2/0/T/boundaryField/chamber_wall_up"
const vectorField& Cf = patch().Cf();
            scalarField& field = *this;
            const scalar zOffset = 0.029;
            static const scalar zTab[34] = {
                0.001, 0.011, 0.021, 0.031, 0.041, 0.051, 0.061, 0.071, 0.081, 0.091,
                0.101, 0.111, 0.121, 0.131, 0.141, 0.151, 0.161, 0.171, 0.181, 0.191,
                0.201, 0.211, 0.221, 0.231, 0.241, 0.251, 0.261, 0.271, 0.281, 0.291,
                0.301, 0.311, 0.321, 0.331
            };
            static const scalar TupTab[34] = {
                294.9, 296.1, 300.9, 323.4, 337.8, 347.3, 355.3, 362.8, 369.9, 376.9,
                383.6, 390.1, 396.5, 402.7, 408.6, 414.4, 420.0, 425.4, 430.6, 435.6,
                440.5, 445.2, 449.7, 454.0, 458.1, 462.1, 466.0, 469.7, 473.4, 477.6,
                483.3, 494.2, 527.3, 617.2
            };
            forAll(Cf, i)
            {
                scalar z = Cf[i].z() + zOffset;
                if (z <= zTab[0]) field[i] = TupTab[0];
                else if (z >= zTab[33]) field[i] = TupTab[33];
                else
                {
                    for (int k=0; k<33; ++k)
                    {
                        if (z >= zTab[k] && z <= zTab[k+1])
                        {
                            scalar frac = (z - zTab[k]) / (zTab[k+1] - zTab[k]);
                            field[i] = TupTab[k] + frac * (TupTab[k+1] - TupTab[k]);
                            break;
                        }
                    }
                }
            }
//}}} end code

    this->parent_bctype::updateCoeffs();
}


// ************************************************************************* //

