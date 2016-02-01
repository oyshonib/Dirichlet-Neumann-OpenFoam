    Info<< "Reading field U\n" << endl;

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );


//volScalarField x = exp(U.mesh().C() & vector(1,0,0));
dimensionedVector one
(
 "one",
 dimensionSet(0, -1, 0, 0, 0, 0 ,0),
 vector(1,0,0)
 );

dimensionedVector adjust
(
 "adjust",
 dimensionSet(0, -1, -1, 0, 0, 0 ,0),
 vector(1,0,0)
 );
volScalarField x = exp(U.mesh().C() & one);
//volVectorField f(x * ); 
volVectorField f(x * adjust);

//Info << x << endl;

//  Info<< "Reading diffusivity D\n" << endl;
/*    Info<< "Reading speed c\n" << endl;

    dimensionedScalar c
    (
        transportProperties.lookup("c")
    );
*/



//#   include "createPhi.H"
