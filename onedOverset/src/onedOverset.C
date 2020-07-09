/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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
    onedOverset


Description
    -Poisson equation solver in 1D using an overset mesh
    -Solves the Poisson Equation Del^T = Del.q analogous to the pressure-Poisson equation
    -Conservative correction by D.Chandar(New approach -**unpublished**, work in progress)


    -Overset search is based on the following

    @conference{behley2015icra,
     author = {Jens Behley and Volker Steinhage and Armin B. Cremers},
      title = {{Efficient Radius Neighbor Seach in Three-dimensional Point Clouds}},
  booktitle = {Proc. of the IEEE International Conference on Robotics and Automation (ICRA)},
       year = {2015}
       https://github.com/jbehley/octree
       
       
}


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "regionProperties.H"
#include "simpleControl.H"
#include "simpleMatrix.H"
#include "octree.H"
#define  CELLMANAGERINTERP      -1E+30
#define  CELLMANAGERINTERPCOMP  -1E+29

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
class Point3f
{
    public:
        
    Point3f(double x, double y, double z) : x(x), y(y), z(z)
    {
  
        
    }

    double x, y, z;
};

class overlap
{
    public:
    
    overlap(){}    
        
    overlap(const fvMesh & mesh) : mesh_(&mesh)
    {
        donorCells.resize(mesh.nCells());
        donorCells=-1;
        
        donorCellCells.resize(mesh.nCells());
        iweights.resize(mesh.nCells());
    }
    
    ~overlap()
    {
        octree.clear();
    }
        
    // Initialize the octree search based on Jens Behley
    unibn::Octree<Point3f>  octree;
    unibn::OctreeParams     params; 
    
    // Cell-Centres to search from (donor mesh)
    std::vector<Point3f> points;
    
    // Main donor cells for the current mesh (-1) if no donor
    labelList donorCells;
    
    // Main donor cells + neighbouring cells
    List<labelList> donorCellCells;
    
    // Interpolation weights
    List<scalarList> iweights;
    
    // Current mesh
    const fvMesh* mesh_;
    
    // Overset interfaces
    labelList oversetInterface;
    
    // Access
    const fvMesh & mesh()
    {
            return *mesh_;
    }
    
    inline void initialize()
    {
        if ( points.size() > 0 )
            octree.initialize(points);
        else
            Info << "Search not initalized yet" << endl;
    }
    
    inline void searchFromMeshAndAssign(const fvMesh & mesh)
    {
        //Set the donor search points
        forAll(mesh.C(),i)
        {
            points.push_back(Point3f(mesh.C()[i].x(),mesh.C()[i].y(),mesh.C()[i].z()));
        }
    }
    
    inline void donorSearch(const fvMesh & donorMesh)
    {
        std::vector<uint32_t> results;
        forAll(mesh().C(),i)
        {
                        
            // Search around dist
            double dist = mesh().V()[i]/(0.1*0.1);
            octree.radiusNeighbors<unibn::L2Distance<Point3f> >(Point3f(mesh().C()[i].x(),mesh().C()[i].y(),mesh().C()[i].z()),dist,results);
       
            // Now find the closest among the results;
            double minDist=1E+10;
        
            //Index of the closest cell.
            int minDistI=-1;
        
            for ( uint32_t j = 0 ; j < results.size() ; j++ )
            {
                double localDist = std::sqrt(unibn::L2Distance<Point3f>::compute(points[results[j]],Point3f(mesh().C()[i].x(),mesh().C()[i].y(),mesh().C()[i].z()))) ;
                if (localDist <= minDist )
                {
                    minDist = localDist;
                    minDistI = results[j];
                }
    
            }
            donorCells[i] = minDistI;
        
            if ( donorCells[i] >= 0 )
            {
                donorCellCells[i].resize(1+donorMesh.cellCells()[donorCells[i]].size());
                label dSize = donorCellCells[i].size();
                
                // Main donor
                donorCellCells[i][0] = donorCells[i];
                
                // Neighbor donors
                for(label j = 1; j < dSize ; j++)
                {
                    donorCellCells[i][j] = donorMesh.cellCells()[donorCells[i]][j-1];
                }
                
            }
            
            
        } 
    }
    
    inline void getInterpolationWeights(const fvMesh & donorMesh, volScalarField & donorBlanks)
    {
        forAll(mesh().C(),i)
        {
            if ( donorCells[i] >=0 )
            {
                iweights[i].resize(donorCellCells[i].size());
                iweights[i]=0.0;
                scalar sumWeights=0.0;
                
                // Distance between main donor and the recipient cell that is 'i'
                scalar d2 =  (mesh().C()[i].x()- donorMesh.C()[donorCells[i]].x())*(mesh().C()[i].x()- donorMesh.C()[donorCells[i]].x()) +
                             (mesh().C()[i].y()- donorMesh.C()[donorCells[i]].y())*(mesh().C()[i].y()- donorMesh.C()[donorCells[i]].y()) +
                             (mesh().C()[i].z()- donorMesh.C()[donorCells[i]].z())*(mesh().C()[i].z()- donorMesh.C()[donorCells[i]].z());
                if ( d2 <= 1E-10 )
                {
                    iweights[i][0]=1.0;
                }
                else
                {
                    iweights[i][0] = 1.0/(d2);
                    sumWeights+=iweights[i][0];
                    
                    for(label j = 1; j < donorCellCells[i].size(); j++)
                    {
                        d2 =  (mesh().C()[i].x()- donorMesh.C()[donorCellCells[i][j]].x())*(mesh().C()[i].x()- donorMesh.C()[donorCellCells[i][j]].x()) +
                              (mesh().C()[i].y()- donorMesh.C()[donorCellCells[i][j]].y())*(mesh().C()[i].y()- donorMesh.C()[donorCellCells[i][j]].y()) +
                              (mesh().C()[i].z()- donorMesh.C()[donorCellCells[i][j]].z())*(mesh().C()[i].z()- donorMesh.C()[donorCellCells[i][j]].z());
             
                        // Skip if nei donors are interp cells
                        if ( donorBlanks[donorCellCells[i][j]] > 0 )      
                        {      
                            iweights[i][j] = 1.0/(d2);
                        }
                        else
                            iweights[i][j]=0.0;
                        
                        sumWeights+=iweights[i][j];
                        
                    }
                    
                    
                    for(label j = 0; j < donorCellCells[i].size(); j++)
                    {
                        iweights[i][j]/=sumWeights;
                    }
                    
                }
                             
                             
            }
        }
    }
    
    inline void finalizeClassification(const fvMesh & donorMesh, const overlap* donorConnect, const volScalarField & cellmanager)
    {
        // Check if interpolation cells's donors are themselves interpolation cells
        forAll(mesh().C(),i)
        {   
            bool flip=false;
            // Is a possible interpolation cell (exclude mandatory interpolation cells)
            if ( donorCells[i] >=0 && cellmanager.primitiveField()[i] > CELLMANAGERINTERPCOMP)
            {
                //Check the primary donor;
                //label primaryDonor = donorCells[i];
                // If donor cell is an interpolation cell
                 for (label j = 0 ; j < donorCellCells[i].size(); j++ )
                 {
                    if ( donorConnect->donorCells[donorCellCells[i][j]] >=0)
                    {
                        flip=true;
                        
                    }
                 }
                 
                 if ( flip)
                 {
                    donorCells[i]=-1;
                    donorCellCells[i]=-1;   
                 }
                 
            }
        }
    }
    
    inline void setCellBlanks(volScalarField & cellBlank)
    {
        forAll(mesh().C(),i)
        {
            cellBlank[i]*= !(donorCells[i]>=0);
        }
    }
    
    inline void incrementOverlap(volScalarField & cellBlank,volScalarField & cellmanager)
    {
        labelList cellblankcopy(mesh().nCells());
        forAll(cellblankcopy, i)
        {
            cellblankcopy[i]=cellBlank[i];
        }
        
        forAll(mesh().faces(),i)
        {
            if (mesh().isInternalFace(i))
            {
                label own = mesh().owner()[i];
                label nei = mesh().neighbour()[i];
                if (cellBlank[own]==0 && cellBlank[nei]==1 && cellmanager.primitiveFieldRef()[own] > CELLMANAGERINTERPCOMP)
                {
                    cellblankcopy[own]=1;
                    donorCells[own]=-1;
                    donorCellCells[own]=-1;
                }
                else if (cellBlank[own]==1 && cellBlank[nei]==0 && cellmanager.primitiveFieldRef()[nei] > CELLMANAGERINTERPCOMP)
                {
                    cellblankcopy[nei]=1;
                    donorCells[nei]=-1;
                    donorCellCells[nei]=-1;
                }
                    
            }
        }
        
        // Copy back
        forAll(cellblankcopy, i)
        {
            cellBlank.primitiveFieldRef()[i]=cellblankcopy[i];
        }
        
        
    }
    
    inline  void setOversetInterfaces(volScalarField & cellBlank)
    {
        const labelList & own= mesh().owner();
        const labelList & nei= mesh().neighbour();
    
        labelHashSet oint;
        
        forAll(mesh().faces(),i)
        {
            if (mesh().isInternalFace(i))
            {
                bool ownInterp = (cellBlank[own[i]]==0 && cellBlank[nei[i]]==1) ;
                bool neiInterp = (cellBlank[own[i]]==1 && cellBlank[nei[i]]==0) ;
                
                if ( ownInterp || neiInterp )
                {
                    oint.insert(i);
                }
                
            }
            
        }
        
        oversetInterface=oint.toc();
        
    }
    
    inline label locateSingleOversetInterfaceCell(volScalarField & cellBlank)
    {
        if ( cellBlank[mesh().owner()[oversetInterface[0]]] == 1 )
        {
            //Info << "From locate : " << "interface = " << oversetInterface << " own = " << mesh().owner()[ oversetInterface[0] ] << endl;
            return mesh().neighbour()[oversetInterface[0]];
        }
        else
            return mesh().owner()[oversetInterface[0]];
    }
    
    
    
};

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Laplace equation solver for a scalar quantity."
    );

    argList::addBoolOption
    (
		"cons",
		"Consistent Flux Correction"
    );


    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMeshes.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    // Hard-coded for now:
    const label nMeshes = 2;
    label nCellsTotal = 0;
    label alphaMesh =0;
    
    // First need to flag the cells adjacent to overset boundaries as mandatory
    for( label mi = 0 ; mi < nMeshes ; mi++)
    {
        const fvMesh & mesh = componentMesh[mi];
        forAll(mesh.boundary(),patchI)
        {
            if ( mesh.boundaryMesh()[patchI].physicalType()=="Overset")
            {
            
                forAll(mesh.boundary()[patchI],faceI)
                {
                    cellManager[mi].primitiveFieldRef()[mesh.boundary()[patchI].faceCells()[faceI]]=CELLMANAGERINTERP;
                }
            
            }
            
        }
        
    }
    
    
    // Overset object
    overlap** connect = new overlap* [nMeshes];
    
    for ( label i = 0 ; i < nMeshes; i++ )
    {
        connect[i]= new overlap(componentMesh[i]);
        
            for ( label j = 0 ; j < nMeshes ; j++ )
            {
                if ( i!=j)
                {
                    connect[i]->searchFromMeshAndAssign(componentMesh[j]);
                    connect[i]->initialize();
                    connect[i]->donorSearch(componentMesh[j]);
                }
            }
            
            nCellsTotal+=componentMesh[i].nCells();
    }
    
    // Donor-interp finalize Classification
    for ( label i = 0 ; i < nMeshes; i++ )
    {        
            for ( label j = 0 ; j < nMeshes ; j++ )
            {
                if ( i!=j)
                {
                    connect[i]->finalizeClassification(componentMesh[j],connect[j],cellManager[i]);
                }
            }
            
    }
    
   
    
    
    // Set setCellBlanks
    for ( label i = 0 ; i < nMeshes; i++ )
    {  
            connect[i]->setCellBlanks(cellBlank[i]);
    }
    
    // Increment overlap if required
    for ( label i = 0 ; i < nMeshes; i++ )
    {  
      for (label nTimes = 0 ; nTimes < 1 ; nTimes++)
      {    
        connect[i]->incrementOverlap(cellBlank[i],cellManager[i]);
      }
          
    }
    
    // Set the interfaces
    for ( label i = 0 ; i < nMeshes; i++ )
    {  
            connect[i]->setOversetInterfaces(cellBlank[i]);
    }
    
    
    // Interpolation weights
    for ( label i = 0 ; i < nMeshes; i++ )
    {  
            for ( label j = 0 ; j < nMeshes ; j++ )
            {
                if ( i!=j)
                {
                    connect[i]->getInterpolationWeights(componentMesh[j],cellBlank[j]);
                }
            }
            
    }
    
    
    
   
    forAll(connect[0]->donorCells, i)
    {
        Info << "cell[" <<i<<"] :: donors = " << connect[0]->donorCellCells[i] << "   Weights = " << connect[0]->iweights[i] << endl;
    }
    
    Info <<" " << endl;
    
    forAll(connect[1]->donorCells, i)
    {
        Info << "cell[" <<i<<"] :: donors = " << connect[1]->donorCellCells[i] << "   Weights = " << connect[1]->iweights[i] << endl;
    }
    
    

    PtrList<volScalarField> RHSComposite(nMeshes);
    PtrList<fvScalarMatrix> TEqnComposite(nMeshes);
    PtrList<surfaceScalarField> phiComposite(nMeshes);
    
    for( label mi = 0 ; mi < nMeshes ; mi++)
    {
        const fvMesh & mesh = componentMesh[mi];
        // Define the RHS 
        surfaceVectorField rhsFace("cf",mesh.Cf()*scalar(1));
    
        // Internal faces
        forAll( mesh.faces(), i)
        {
            if  (mesh.isInternalFace(i))
            {
                rhsFace[i].x()=2.0*M_PI*Foam::cos(2*M_PI*rhsFace[i].x());
                rhsFace[i].y()=0;
                rhsFace[i].z()=0;
            }
        }
    
        surfaceVectorField::Boundary & brhs = rhsFace.boundaryFieldRef(); 
    
        // Boundary faces
        forAll(mesh.boundary(),patchI)
        {
        
            fvsPatchVectorField  & bTestVec = brhs[patchI];
            vectorField assign(bTestVec.size());
        
            const fvPatch & pp = bTestVec.patch();
            vectorField fc(pp.Cf());
        
            forAll(mesh.boundary()[patchI],faceI)
            {
                assign[faceI].x()=2.0*M_PI*Foam::cos(2*M_PI*fc[faceI].x());
                assign[faceI].y()=0;
                assign[faceI].z()=0;
            }
        
            bTestVec==assign;
        }
    

        surfaceScalarField phi(rhsFace & mesh.Sf() );
    
        RHSComposite.set(mi, new volScalarField("RHS",fvc::div(phi)));
        
        TEqnComposite.set(mi, new fvScalarMatrix(fvm::laplacian(DT,TComposite[mi])));
        
        phiComposite.set(mi, new surfaceScalarField("phi",phi));
        
    }
    
    // Frame the Global Matrix of all meshes
    // First need to create local mesh to global matrix mapping
    labelList cumulativeSum(nMeshes);
    cumulativeSum[0] = 0;
    
    forAll(cumulativeSum,j)
    {
        if (j > 0)
        {
            cumulativeSum[j]=cumulativeSum[j-1]+componentMesh[j-1].nCells();
        }
    }
    
    
    simpleMatrix<scalar> Amat(nCellsTotal); 
    
    // Set the diagonal and RHS
    for(label i = 0 ; i < nMeshes ; i++)
    {
        const fvMesh & mesh = componentMesh[i];
        const fvScalarMatrix & TEqn = TEqnComposite[i];
        const surfaceScalarField & deltaCoeffs = mesh.surfaceInterpolation::deltaCoeffs();
        
        forAll(TComposite[i],rowI)
        {
            label gRow = rowI + cumulativeSum[i];
            Amat.source()[gRow] = mesh.V()[rowI]*RHSComposite[i][rowI];
            Amat[gRow][gRow]=TEqn.diag()[rowI];
                    
        }
        
        
        // Upper and Lower Coefficients
        for(label faceI= 0 ; faceI < mesh.nInternalFaces() ; faceI++)
        {
            label nei = TEqn.lduAddr().lowerAddr()[faceI] + cumulativeSum[i];
            label own = TEqn.lduAddr().upperAddr()[faceI] + cumulativeSum[i];
        
            Amat[own][nei] = TEqn.upper()[faceI];
            Amat[nei][own] = TEqn.upper()[faceI];
        }
        
        
        // Boundary Conditions
        forAll(TComposite[i].boundaryFieldRef(),patchI)
        {
            const fvPatch & pp = TComposite[i].boundaryFieldRef()[patchI].patch();
            forAll(pp,faceI)
            {
                label cellI = pp.faceCells()[faceI] + cumulativeSum[i];
                Amat[cellI][cellI] += TEqn.internalCoeffs()[patchI][faceI];
                Amat.source()[cellI]+= TEqn.boundaryCoeffs()[patchI][faceI];
            }
        }
        
        
        // Interpolation Coefficients (Overset)
        forAll(TComposite[i],rowI)
        {   
            // On interpolation cells only
            if ( connect[i]->donorCells[rowI]>=0 )
            {
                label gRow = rowI + cumulativeSum[i];
                
                // Make zero all columns
                for (label jc = 0 ; jc < nCellsTotal ; jc++)
                {
                    Amat[gRow][jc]=0;
                }
                
                // Diagonal is one
                Amat[gRow][gRow]=1;//TEqn.diag()[rowI];
                
                // No RHS
                Amat.source()[gRow] = 0;
                
                // Insert weights as coefficients
                forAll(connect[i]->donorCellCells[rowI], in_)
                {
                    label gCol = connect[i]->donorCellCells[rowI][in_]+cumulativeSum[i==0?1:0];
                    Amat[gRow][gCol] = -connect[i]->iweights[rowI][in_];
                }
            }
            
            
        }
     
    bool conservativeFix = false;

    if ( args.found("cons") )
    {
	conservativeFix = true;    
    }    

    if ( conservativeFix )
    {
        // overwrite coefficients for conservation only on the alpha mesh
        if (i == alphaMesh)
        {
            label gRow = connect[i]->locateSingleOversetInterfaceCell(cellBlank[i]);
            
           
            for ( label colIndex = 0; colIndex < nCellsTotal ; colIndex++ )
            {
                    Amat[gRow][colIndex]=0;
            }
            
            Amat.source()[gRow] = 0;
            //Info << " alphaCell = " << gRow << endl;
            //Info << " interface = " << connect[i]->oversetInterface << endl;
            
            
            // All interfaces on the current mesh
            forAll(connect[i]->oversetInterface,of_)
            {
                label faceID = connect[i]->oversetInterface[of_];
                label own    = componentMesh[i].owner()[faceID];
                label nei    = componentMesh[i].neighbour()[faceID];
            
                
                // Owner is c0 (flip)-refer notes
                if ( own == gRow )
                    Amat.source()[gRow] -= phiComposite[i][faceID];
                else
                    Amat.source()[gRow] += phiComposite[i][faceID];
        
                if ( of_ == 0 )
                {
                    // Owner is c0. So flip sign of nei and gRow (usually nei has +ve coeff on the snGrad) - refer notes
                    // Later change mesh to componentMesh[alphaMesh]
                    if (own == gRow)
                    {
                        Amat[gRow][nei]  +=  -deltaCoeffs[faceID]*mesh.magSf()[faceID];
                        Amat[gRow][gRow] +=   deltaCoeffs[faceID]*mesh.magSf()[faceID];
                    }
                    else
                    {
                        Amat[gRow][own]   +=  -deltaCoeffs[faceID]*mesh.magSf()[faceID];
                        Amat[gRow][gRow]  +=  deltaCoeffs[faceID]*mesh.magSf()[faceID];//this is the nei-so no flipping
                    }
                    
                }
                else
                {
                    if ( cellBlank[i][own] == 0 )
                    {
                        Amat[gRow][nei] += -deltaCoeffs[faceID]*mesh.magSf()[faceID];
                        Amat[gRow][own] +=  deltaCoeffs[faceID]*mesh.magSf()[faceID];
                    }
                    else
                    {
                        Amat[gRow][nei] +=  deltaCoeffs[faceID]*mesh.magSf()[faceID];
                        Amat[gRow][own] += -deltaCoeffs[faceID]*mesh.magSf()[faceID];
                    }
                     
                }
                
            
            }
            
            // All interfaces on other meshes
            for (int jMesh = 0 ; jMesh < nMeshes ; jMesh++ )
            {
                
                if (jMesh!=alphaMesh)
                {
                    const fvMesh & meshj = componentMesh[jMesh];
                    const surfaceScalarField & deltaCoeffsj = componentMesh[jMesh].surfaceInterpolation::deltaCoeffs();
                    
                    for (int oi=0; oi < connect[jMesh]->oversetInterface.size() ; oi++ )
                    {
                        label faceID = connect[jMesh]->oversetInterface[oi];
                        label own = componentMesh[jMesh].owner()[faceID];
                        label nei = componentMesh[jMesh].neighbour()[faceID];
                        
                        label ownCol = own + cumulativeSum[i==0?1:0];
                        label neiCol = nei + cumulativeSum[i==0?1:0];
                        
                        if (cellBlank[jMesh][own] == 1)
                        {
                            Amat[gRow][neiCol] +=   deltaCoeffsj[faceID]*meshj.magSf()[faceID];
                            Amat[gRow][ownCol] +=   -deltaCoeffsj[faceID]*meshj.magSf()[faceID];
                        
                            Amat.source()[gRow] += phiComposite[jMesh][faceID];
                        }
                        else
                        {
                            Amat[gRow][neiCol] +=  -deltaCoeffsj[faceID]*meshj.magSf()[faceID];
                            Amat[gRow][ownCol] +=   deltaCoeffsj[faceID]*meshj.magSf()[faceID];
                        
                            Amat.source()[gRow] -= phiComposite[jMesh][faceID];
                        }
                    }
                }
                    
            }
            
            
        }// end i==alphamesh
    }  
    
}
         
    OFstream matStream("matrix.txt");     
    for ( label i = 0 ;i < nCellsTotal ; i++ )
    {
        for (label j = 0 ; j < nCellsTotal ; j++)
        {
            if ( mag(Amat(i,j)) > 0 )
            {
                matStream << "(" << i <<","<<j<<")      " << Amat(i,j) <<"      "; 
            }
        }
        matStream << "      Source = " << Amat.source()[i];
        
        matStream << endl;
    }
        
    scalarField solution(Amat.LUsolve());
    
    
    label m=0;
    scalar l2Error=0;
    
    for(label i = 0 ; i < nMeshes ; i++)
    {
        OFstream outStream("solution."+Foam::name(i));
        
        forAll(TComposite[i],j)
        {
            TComposite[i].primitiveFieldRef()[j]=solution[m];
            TExact[i].primitiveFieldRef()[j]=Foam::sin(2*M_PI*componentMesh[i].C()[j].x());
            scalar localError=cellBlank[i][j]*(TExact[i].primitiveFieldRef()[j]-TComposite[i].primitiveFieldRef()[j]);
            
            l2Error += (localError*localError);
            
            outStream << componentMesh[i].C()[j].x() <<"  "<<TComposite[i].primitiveFieldRef()[j]<<"  "<<TExact[i].primitiveFieldRef()[j]<< "  " <<cellBlank[i][j]<<endl;
            
            m++;
        }
        
        TComposite[i].write();
        TExact[i].write();
        cellBlank[i].write();
        
    }
    
    // Check for Conservation Error: (sum(gradp)-sum(rhs)) on all boundaries except mandatory overset boundaries
    scalar unbalancedSum=0;
    for (label i = 0 ; i < nMeshes ; i++)
    {
        const fvMesh & mesh = componentMesh[i];
        
        surfaceScalarField gradTDotS(fvc::snGrad(TComposite[i])*mesh.magSf());
        
        forAll(mesh.boundary(),patchI)
        {
            if ( mesh.boundaryMesh()[patchI].physicalType()!="Overset")
            {
                unbalancedSum+= sum(gradTDotS.boundaryFieldRef()[patchI]-phiComposite[i].boundaryFieldRef()[patchI]);
            }
        }
        
    }
    
    
    Info << "Mesh.size = " << nCellsTotal << nl
         <<" L2 error  = " << Foam::sqrt(l2Error/nCellsTotal) << nl
         <<" C error   = " << unbalancedSum<<endl;
    
    
    delete [] connect;
    
    

    return 0;
}


// ************************************************************************* //
