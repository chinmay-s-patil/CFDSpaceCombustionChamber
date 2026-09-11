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
// SHA1 = bec84168b91024544db16f1c0f1bfd3ebf37d325
//
// unique function name that can be checked if the correct library version
// has been loaded
extern "C" void chamberSideT_bec84168b91024544db16f1c0f1bfd3ebf37d325(bool load)
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
    chamberSideTFixedValueFvPatchScalarField
);

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
chamberSideTFixedValueFvPatchScalarField::
chamberSideTFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    parent_bctype(p, iF)
{
    if (false)
    {
        printMessage("Construct chamberSideT : patch/DimensionedField");
    }
}


Foam::
chamberSideTFixedValueFvPatchScalarField::
chamberSideTFixedValueFvPatchScalarField
(
    const chamberSideTFixedValueFvPatchScalarField& rhs,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    parent_bctype(rhs, p, iF, mapper)
{
    if (false)
    {
        printMessage("Construct chamberSideT : patch/DimensionedField/mapper");
    }
}


Foam::
chamberSideTFixedValueFvPatchScalarField::
chamberSideTFixedValueFvPatchScalarField
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
        printMessage("Construct chamberSideT : patch/dictionary");
    }
}


Foam::
chamberSideTFixedValueFvPatchScalarField::
chamberSideTFixedValueFvPatchScalarField
(
    const chamberSideTFixedValueFvPatchScalarField& rhs
)
:
    parent_bctype(rhs),
    dictionaryContent(rhs)
{
    if (false)
    {
        printMessage("Copy construct chamberSideT");
    }
}


Foam::
chamberSideTFixedValueFvPatchScalarField::
chamberSideTFixedValueFvPatchScalarField
(
    const chamberSideTFixedValueFvPatchScalarField& rhs,
    const DimensionedField<scalar, volMesh>& iF
)
:
    parent_bctype(rhs, iF)
{
    if (false)
    {
        printMessage("Construct chamberSideT : copy/DimensionedField");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::
chamberSideTFixedValueFvPatchScalarField::
~chamberSideTFixedValueFvPatchScalarField()
{
    if (false)
    {
        printMessage("Destroy chamberSideT");
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::
chamberSideTFixedValueFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        printMessage("updateCoeffs chamberSideT");
    }

//{{{ begin code
    #line 134 "/home/lavender/OpenFoamUbu/CFDSpaceCombustionChamber/openFOAMTry/simFile/0/T/boundaryField/chamber_wall_side"
const vectorField& Cf = patch().Cf();
            scalarField& field = *this;
            const scalar zOffset = 0.029;
            static const scalar zTab[34] = {
                0.001, 0.011, 0.021, 0.031, 0.041, 0.051, 0.061, 0.071, 0.081, 0.091,
                0.101, 0.111, 0.121, 0.131, 0.141, 0.151, 0.161, 0.171, 0.181, 0.191,
                0.201, 0.211, 0.221, 0.231, 0.241, 0.251, 0.261, 0.271, 0.281, 0.291,
                0.301, 0.311, 0.321, 0.331
            };
            static const scalar TsideTab[34] = {
                294.9, 296.1, 300.8, 323.3, 337.8, 347.2, 355.2, 362.7, 369.8, 376.8,
                383.5, 390.0, 396.3, 402.5, 408.5, 414.2, 419.8, 425.2, 430.4, 435.4,
                440.3, 444.9, 449.4, 453.8, 457.9, 461.9, 465.7, 469.4, 473.2, 477.2,
                482.8, 492.5, 517.5, 582.7
            };
            forAll(Cf, i)
            {
                scalar z = Cf[i].z() + zOffset;
                if (z <= zTab[0]) field[i] = TsideTab[0];
                else if (z >= zTab[33]) field[i] = TsideTab[33];
                else
                {
                    for (int k=0; k<33; ++k)
                    {
                        if (z >= zTab[k] && z <= zTab[k+1])
                        {
                            scalar frac = (z - zTab[k]) / (zTab[k+1] - zTab[k]);
                            field[i] = TsideTab[k] + frac * (TsideTab[k+1] - TsideTab[k]);
                            break;
                        }
                    }
                }
            }
//}}} end code

    this->parent_bctype::updateCoeffs();
}


// ************************************************************************* //

