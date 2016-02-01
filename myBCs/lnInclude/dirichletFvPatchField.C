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

#include "dirichletFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
dirichletFvPatchField<Type>::dirichletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
  fvPatchField<Type>(p, iF),
     cyclicPatch_(refCast<const cyclicFvPatch>(p))
  //cyclicPatch_(const cyclicFvPatch(p))
{}


template<class Type>
dirichletFvPatchField<Type>::dirichletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
  fvPatchField<Type>(p, iF, dict, true),
    cyclicPatch_(refCast<const cyclicFvPatch>(p))
{}


template<class Type>
dirichletFvPatchField<Type>::dirichletFvPatchField
(
    const dirichletFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
  fvPatchField<Type>(ptf, p, iF, mapper),
  cyclicPatch_(refCast<const cyclicFvPatch>(p))
{}


template<class Type>
dirichletFvPatchField<Type>::dirichletFvPatchField
(
    const dirichletFvPatchField<Type>& ptf
)
:
  fvPatchField<Type>(ptf),
  cyclicPatch_(ptf.cyclicPatch_)
{}


template<class Type>
dirichletFvPatchField<Type>::dirichletFvPatchField
(
    const dirichletFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
  fvPatchField<Type>(ptf, iF),
  cyclicPatch_(ptf.cyclicPatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Type>
tmp<Field<Type> > dirichletFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<Type> >
    (
        new Field<Type>(this->size(), pTraits<Type>::zero)
    );
}


template<class Type>
tmp<Field<Type> > dirichletFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return *this;
}


template<class Type>
tmp<Field<Type> > dirichletFvPatchField<Type>::gradientInternalCoeffs() const
{
    return -pTraits<Type>::one*this->patch().deltaCoeffs();
}


template<class Type>
tmp<Field<Type> > dirichletFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    return this->patch().deltaCoeffs()*(*this);
}


template<class Type>
void dirichletFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
