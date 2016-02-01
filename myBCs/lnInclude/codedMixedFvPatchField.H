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

Class
    Foam::mixedFvPatchField

Description
    Foam::mixedFvPatchField

SourceFiles
    mixedFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef codedMixedFvPatchField_H
#define codedMixedFvPatchField_H

//#include "fvPatchField.H"
#include "mixedFvPatchFields.H"
#include "codedBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


  // Forward declaration of classes
  class dynamicCode;
  class dynamicCodeContext;
  class IOdictionary;

/*---------------------------------------------------------------------------*\
                           Class mixedFvPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class mixedFvPatchField
:
    public mixedFvPatchField<Type>,
    public codedBase
{
    // Private data
  //- Dictionary contents for the boundary condition
  mutable dictionary dict_;

  const word redirectType_;

  mutable autoPtr<mixedFvPatchField<Type> > redirectPatchFieldPtr_;

  // Private Member Functions

  const IOdictionary& dict() const;

  //- Set the rewrite vars controlling the Type
  static void setFieldTemplates(dynamicCode& dynCode);

  //- get the loaded dynamic libraries
  virtual dlLibraryTable& libs() const;

  //- adapt the context for the current object
  virtual void prepare(dynamicCode&, const dynamicCodeContext&) const;

  // Return a description (type + name) for the output
  virtual string description() const;

  // Clear the ptr to the redirected object
  virtual void clearRedirect() const;

  // Get the dictionary to initialize the codeContext
  virtual const dictionary& codeDict() const;




public:

  // Static data members

  //- Name of the C code template to be used
  static const word codeTemplateC;

  //- Name of the H code template to be used
  static const word codeTemplateH;




    //- Runtime type information
    TypeName("codedMixed");


    // Constructors

        //- Construct from patch and internal field
        codedMixedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        codedMixedFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given mixedFvPatchField onto a new patch
        codedMixedFvPatchField
        (
            const mixedFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        codedMixedFvPatchField
        (
            const mixedFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new mixedFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        codedMixedFvPatchField
        (
            const mixedFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new mixedFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions
  //- Get reference to the underlying patchField
  const mixedFvPatchField<Type>& redirectPatchField() const;

  //- Update the coefficients associated with the patch field
  virtual void updateCoeffs();

  //- Evaluate the patch field
  //  This is only needed to set the updated() flag of the redirectType
  //  to false.
        virtual void evaluate
        (
	 const Pstream::commsTypes commsType=Pstream::blocking
	 );

  //- Write
  virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "codedMixedFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
