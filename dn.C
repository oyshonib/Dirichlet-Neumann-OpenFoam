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

Application
    scalarTransportFoam

Description
    Solves a transport equation for a passive scalar

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include <mpi.h>

#include "OPstream.H"
#include "IPstream.H"

#include "mixedFvPatchField.H"
//#include "simpleControl.H"


/*#if defined(WM_SP)
#   define MPI_SCALAR MPI_FLOAT
#elif defined(WM_DP)
#   define MPI_SCALAR MPI_DOUBLE
#endif*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

  //  simpleControl simple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    //    #include "CourantNo.H"


    //vector lambda(0,0,0); //createfields
    vectorField lambda(1, vector(0,0,0));

    //Info<< "lambda " << lambda << endl;

    scalar theta = 0.6; //0.6;










    //label patchID_to = mesh.boundaryMesh().findPatchID("to");
    //const fvPatchVectorField fvp = U.boundaryField()[patchID_to];



//fileName fvp(runTime.path()/runTime.timeName());
//OFstream os(fvp);


    //while (simple.loop())
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

	#       include "readSIMPLEControls.H"



         

    //modify processorBC
    //identify patchID
    label patchID_from = mesh.boundaryMesh().findPatchID("from");
    label patchID_to = mesh.boundaryMesh().findPatchID("to");

    


//declaring Field<vector> given size
vectorField u1(1);

//surfaceVectorField second_exchange(vector(0,0,0));
//fvPatchField<vector> second_exchange = U.boundaryField()[patchID_to];
vectorField u2(1);




if (Pstream::myProcNo() == 0) {

    
   

    //receive first_exchange from processor 1    
    IPstream::read(Pstream::blocking, 1,
         reinterpret_cast<char*>(u2.begin()), 
         u2.byteSize());


    //updating velocity at processor patch
    fixedValueFvPatchVectorField& procVel = 
      refCast<fixedValueFvPatchVectorField>(U.boundaryField()[patchID_from]);


     u1 = -(U.boundaryField()[patchID_from]).snGrad();

     OPstream::write(Pstream::nonBlocking, 1,
         reinterpret_cast<char*>(u1.begin()), 
         u1.byteSize());

    
      
      lambda = (theta * u2 ) +
                  (1.0-theta)*lambda;

          
      
        Pout<< "grad(u2): " << u1 << endl;


      forAll(procVel, i)
          {
            procVel[i] = lambda[i];
          }

          //
         

    
}
else if (Pstream::myProcNo() == 1) {
    
    //u2
    u2 = U.boundaryField()[patchID_to];

    //send first_exchange to processor 0
    OPstream::write(Pstream::nonBlocking, 0,
         reinterpret_cast<char*>(u2.begin()), 
         u2.byteSize());








     //Pout<< second_exchange << endl;

    //receive second_exchange from processor 0
    IPstream::read(Pstream::blocking, 0,
         reinterpret_cast<char*>(u1.begin()), 
         u1.byteSize());

     //Pout<< second_exchange << endl;

     
    mixedFvPatchField<vector>& procVel = 
      refCast<mixedFvPatchField<vector> >(U.boundaryField()[patchID_to]);
     
      

     procVel.refGrad() = u1;

} 







    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
   
        solve( -fvm::laplacian(U) == f);
      //U.correctBoundaryConditions();
        }

        runTime.write();







    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
