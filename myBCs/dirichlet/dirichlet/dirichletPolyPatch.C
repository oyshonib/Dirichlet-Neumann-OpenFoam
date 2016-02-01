/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "dirichletPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dirichletPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, dirichletPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, dirichletPolyPatch, dictionary);
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::dirichletPolyPatch::dirichletPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm
)
:
    polyPatch(name, size, start, index, bm)
{}


Foam::dirichletPolyPatch::dirichletPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
)
:
    polyPatch(name, dict, index, bm),
    actualTypeName_(dict.lookup("type")),
    dict_(dict)
{}


Foam::dirichletPolyPatch::dirichletPolyPatch
(
    const dirichletPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    polyPatch(pp, bm),
    actualTypeName_(pp.actualTypeName_),
    dict_(pp.dict_)
{}


Foam::dirichletPolyPatch::dirichletPolyPatch
(
    const dirichletPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    polyPatch(pp, bm, index, newSize, newStart),
    actualTypeName_(pp.actualTypeName_),
    dict_(pp.dict_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dirichletPolyPatch::~dirichletPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dirichletPolyPatch::write(Ostream& os) const
{
    os.writeKeyword("type") << actualTypeName_ << token::END_STATEMENT << nl;
    patchIdentifier::write(os);
    os.writeKeyword("nFaces") << size() << token::END_STATEMENT << nl;
    os.writeKeyword("startFace") << start() << token::END_STATEMENT << nl;

    for
    (
        dictionary::const_iterator iter = dict_.begin();
        iter != dict_.end();
        ++iter
    )
    {
        if
        (
            iter().keyword() != "type"
         && iter().keyword() != "nFaces"
         && iter().keyword() != "startFace"
        )
        {
            iter().write(os);
        }
    }
}


// ************************************************************************* //
