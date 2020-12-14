/*---------------------------------------------------------------------------*\
  A CONSISTENT FLUX CORRECTION APPROACH FOR THE POISSON EQUATION
  ON OVERSET MESHES
                    
  The base code has the following dependencies
  
  1. OpenFOAM-v1812+ (ESI-OpenCFD) 
        (www.openfoam.com)
  
  2. Octree search for overset 
        (https://github.com/jbehley/octree):
            @conference{behley2015icra,
            author = {Jens Behley and Volker Steinhage and Armin B. Cremers},
             title = {{Efficient Radius Neighbor Seach in Three-dimensional Point Clouds}},
         booktitle = {Proc. of the IEEE International Conference on Robotics and Automation (ICRA)},
            year = {2015}
            
  3. The sparse matrix representation classes 
        (https://github.com/uestla/Sparse-Matrix)
  
  4. Saad's Sparse Matrix solvers (BiCGStab) 
        (https://www-users.cs.umn.edu/~saad/software/ITSOL)

        
  Dominic D. J Chandar
  Queen's University of Belfast
  School of Mechanical and Aerospace Engineering
  
  d.chandar@qub.ac.uk, dominic.chandar@gmail.com
  https://github.com/domi0002/geoOverlap
  
  For issues refer to the above GitHub page.

  Updated: July 2020

\*---------------------------------------------------------------------------*/

#include "sparseMatrix.H"
#include "itsol.h"
#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "regionProperties.H"
#include "octree.H"
#include "simpleMatrix.H"
#define  CELLMANAGERINTERP      -1E+10
#define  CELLMANAGERINTERPCOMP  -1E+9

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
class Point3f
/*
    Class to represent the point structure for the octree search
*/
{
    public:
        
    Point3f(double x, double y, double z) : x(x), y(y), z(z)
    {}

    double x, y, z;
};

class overlap
{
    public:
    
    overlap(){}    
        
    overlap(dynamicFvMesh & mesh, const label gridID) : mesh_(&mesh)
    {
        gridID_=gridID;
        
        donorCells.resize(mesh.nCells());
        donorCells=-1;
        
        donorCellCells.resize(mesh.nCells());
        iweights.resize(mesh.nCells());
    }
    
    ~overlap()
    {
        octree.clear();
    }
        
    // Initialize the octree search
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
    dynamicFvMesh* mesh_;
    
    // Current grid
    label gridID_;
    
    // Overset interfaces
    labelList oversetInterface;
    
    //How many interpolation cells
    label nInterpolationCells;
    
    // Access
    dynamicFvMesh & mesh()
    {
            return *mesh_;
    }
    
    const label & gridID()
    {
            return gridID_;
    }
    
    inline void initialize()
    {
        if ( points.size() > 0 )
            octree.initialize(points);
        else
            Info << "Search not initalized yet" << endl;
    }
    
    inline void searchFromMeshAndAssign(const dynamicFvMesh & mesh)
    {
        points.clear();
        
        //Set the donor search points
        forAll(mesh.C(),i)
        {
            points.push_back(Point3f(mesh.C()[i].x(),mesh.C()[i].y(),mesh.C()[i].z()));
        }
    }
     
    inline void donorSearch(const dynamicFvMesh & donorMesh)
    {
        std::vector<uint32_t> results;
        forAll(mesh().C(),i)
        {
            
            
            // Search around dist
            double dist = 10*Foam::sqrt(mesh().V()[i]/(0.1));
            //Info << "search radius = " << dist << endl;
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
            
			
            // Discard the main donor if it is not in cell
            if (donorCells[i]>=0 && !donorMesh.pointInCell(mesh().C()[i],donorCells[i]))
            donorCells[i]=-1;
            
			double minDistance = 0;
			double maxDistance = Foam::sqrt(max(donorMesh.V())/0.1).value();
			int32_t dCell = octree.findNeighbor<unibn::L2Distance<Point3f> >(Point3f(mesh().C()[i].x(),mesh().C()[i].y(),mesh().C()[i].z()),minDistance,maxDistance);
			//Info << " Correct donor = " << donorCells[i] << "  other donor = " << dCell << endl; 
		
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
	
	
	inline void donorSearch(const dynamicFvMesh & donorMesh, label nDonors)
    {
        std::vector<uint32_t> results;
        forAll(mesh().C(),i)
        {
            
            
            // Search around dist
            double dist = 10*Foam::sqrt(mesh().V()[i]/(0.1));
            //Info << "search radius = " << dist << endl;
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
            
			
            // Discard the main donor if it is not in cell
            if (donorCells[i]>=0 && !donorMesh.pointInCell(mesh().C()[i],donorCells[i]))
            donorCells[i]=-1;
            
			double minDistance = 0;
			double maxDistance = Foam::sqrt(max(donorMesh.V())/0.1).value();
			int32_t dCell = octree.findNeighbor<unibn::L2Distance<Point3f> >(Point3f(mesh().C()[i].x(),mesh().C()[i].y(),mesh().C()[i].z()),minDistance,maxDistance);
 
						
			// Proceed only if donorCells[i] > 0 to find other donors			
			if ( donorCells[i] >= 0 )
			{
				if ( results.size() < nDonors )
				{
					Info << " Not enough RBF support cells, Try increasing search radius " << endl;
					Info << " nDonorCells located = " << results.size() << endl;
					std::exit(0);
				}
					
							
				HashTable<label,scalar,Hash<scalar> > resultTable;
				for(label j = 0 ; j < results.size() ; j++ )
				{
					scalar localDist = std::sqrt(unibn::L2Distance<Point3f>::compute(points[results[j]],Point3f(mesh().C()[i].x(),mesh().C()[i].y(),mesh().C()[i].z()))) ;
					resultTable.insert(localDist,results[j]);
				}
				
				scalarList sortedDistance( resultTable.sortedToc() );
				
				// donorCellCells size
				donorCellCells[i].resize(nDonors);
				forAll(donorCellCells[i], j)
				{
					donorCellCells[i][j] = resultTable[ sortedDistance[j] ];
				}

			}
			
			 
            
        } 
    }
	
    
    inline void getInterpolationWeights(const dynamicFvMesh & donorMesh, volScalarField & donorBlanks)
    {
        forAll(mesh().C(),i)
        {
            if ( donorCells[i] >=0 )
            {
                iweights[i].resize(donorCellCells[i].size());
                iweights[i]=0.0;
                scalar sumWeights=0.0;
                
                if (donorBlanks[donorCells[i]] < 1 )
                {
                    Info << " Donor Cell for " << i << " is " << donorCells[i] <<" and is an interpolation cell " << endl;
                }
                
                
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
                    iweights[i][0] = 1.0/d2;
                    sumWeights+=iweights[i][0];
                    
                    for(label j = 1; j < donorCellCells[i].size(); j++)
                    {
                        d2 =  (mesh().C()[i].x()- donorMesh.C()[donorCellCells[i][j]].x())*(mesh().C()[i].x()- donorMesh.C()[donorCellCells[i][j]].x()) +
                              (mesh().C()[i].y()- donorMesh.C()[donorCellCells[i][j]].y())*(mesh().C()[i].y()- donorMesh.C()[donorCellCells[i][j]].y()) +
                              (mesh().C()[i].z()- donorMesh.C()[donorCellCells[i][j]].z())*(mesh().C()[i].z()- donorMesh.C()[donorCellCells[i][j]].z());
             
                        // Skip if nei donors are interp cells
                        if ( donorBlanks[donorCellCells[i][j]] > 0 )      
                        {      
                            iweights[i][j] = 1.0/d2;
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
    
    inline void getInterpolationWeightsLLS(const dynamicFvMesh & donorMesh, volScalarField & donorBlanks)
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
                    // Standard 5x5 matrix
                    simpleMatrix<scalar> iMat(5);
                    
                    // Row 0 (1)
                    for (label j =0 ; j < 5 ; j++)
                    {
                        iMat(0,j) = 1.0;
                        if (donorBlanks[donorCellCells[i][j]] < 1 )
                        {
                            Info << " donor cells is an interpolation cell " << endl;
                        }
                    }
                    
                    // Row 1 (x)
                    for(label j = 0; j < donorCellCells[i].size(); j++)
                    {
                        iMat(1,j) = donorMesh.C()[donorCellCells[i][j]].x();
                    }
                    
                    // Row 2 (y)
                    for(label j = 0; j < donorCellCells[i].size(); j++)
                    {
                        iMat(2,j) = donorMesh.C()[donorCellCells[i][j]].y();
                    }
                    
                    // Row 3 (x^2)
                    for(label j = 0; j < donorCellCells[i].size(); j++)
                    {
                        iMat(3,j) = iMat(1,j)*iMat(1,j);
                    }
                    
                    // Row 4 (y^2)
                    for(label j = 0; j < donorCellCells[i].size(); j++)
                    {
                        iMat(4,j) = iMat(2,j)*iMat(2,j);
                    }
                    
                    iMat.source()[0] = 1;
                    iMat.source()[1] = mesh().C()[i].x();
                    iMat.source()[2] = mesh().C()[i].y();
                    iMat.source()[3] = iMat.source()[1]*iMat.source()[1];                    
                    iMat.source()[4] = iMat.source()[2]*iMat.source()[2];
                    
                    iweights[i] = iMat.LUsolve();
                    
                    
                }
                             
                             
            }
        }
    }
    
    /*
 standard O(n^3) gaussian elimination matrix solver
*/
inline void d_solvec(double **a,double *b,int *iflag,int n)
{
  int i,j,k,l,flag;
  double fact;
  double temp;
  double sum;
  double eps=1e-8;

  
  for(i=0;i<n;i++)
    {
      if (fabs(a[i][i]) < eps)
	{
	  flag=1;
	  for(k=i+1;k<n && flag;k++)
	    {
	      if (a[k][i]!=0)
                {
		  flag=0;
		  for(l=0;l<n;l++)
		    {
		      temp=a[k][l];
		      a[k][l]=a[i][l];
		      a[i][l]=temp;
		    }
		  temp=b[k];
		  b[k]=b[i];
		  b[i]=temp;
                }
	    }
	  if (flag) {*iflag=0;return;}
	}
      for(k=i+1;k<n;k++)
	{
	  if (i!=k)
	    {
	      fact=-a[k][i]/a[i][i];
	      for(j=0;j<n;j++)
		{
		  a[k][j]+=fact*a[i][j];
		}
	      b[k]+=fact*b[i];
	    }
	}
    }

  for(i=n-1;i>=0;i--)
    {
      sum=0;
      for(j=i+1;j<n;j++)
	sum+=a[i][j]*b[j];
      b[i]=(b[i]-sum)/a[i][i];
    }
  *iflag=1;
  return;

}
/*
Interpolate using RBF + linear polynomial 
d=number of dimensions
np=number of points
xcloud=cloud of points xcloud[3*np]
P=target point
weights=interpolation weights
*/
inline void d_interprbf(double *xcloud, double *P,scalarList & weights,int np, int d)
{
  int i,j,n,imin;
  double **M, *B;
  double dd,dmin;
  int neqns=np+d+1;
  int iflag;
  
  
  M=(double **)malloc(sizeof(double *)*neqns);
  B=(double *)malloc(sizeof(double) *neqns);
  for(i=0;i<neqns;i++) M[i]=(double *)malloc(sizeof(double)*neqns);
  
  dmin=1.0e15;
  imin=0;

  for(i=0;i<np;i++)
    {
      for(j=0;j<np;j++)
	{
	  dd=0;
	  for(n=0;n<d;n++) dd+=(xcloud[d*i+n]-xcloud[d*j+n])*(xcloud[d*i+n]-xcloud[d*j+n]);
	  M[i][j]=Foam::exp(-dd);
	}
      M[i][np]=M[np][i]=1.0;
      for(j=np+1,n=0;j<np+d+1;j++,n++) M[i][j]=M[j][i]=xcloud[d*i+n];
      dd=0;
      for(n=0;n<d;n++) dd+=(xcloud[d*i+n]-P[n])*(xcloud[d*i+n]-P[n]);
      if (dd < dmin) {
	imin=i;
	dmin=dd;
      }
      B[i]=Foam::exp(-dd);
    }
  for(i=np;i<np+d+1;i++) 
    {
      for(j=np;j<np+d+1;j++)
	M[i][j]=0.0;
      B[i]=(i==np)?1.0:P[i-(np+1)];
    }
  d_solvec(M,B,&iflag,neqns);
  
  if (iflag==0) 
    {
      for(i=0;i<np;i++) weights[i]=0.0;
      weights[imin]=1.0;
    }
  else 
    {
      for(i=0;i<np;i++)weights[i]=B[i];
    }

  for(i=0;i<neqns;i++) free(M[i]);
  free(M);
  free(B);
}

/* get number of terms in a homogeneous polynomial of a given order for d dimensions */
inline void getPolyCount(int *nterms,int order, int d, int ix)
{
 int i;
 if (order==1) { for (i=ix;i<d;i++) (*nterms)++; return;}
 for(i=ix;i<d;i++) getPolyCount(nterms,order-1,d,i);
}

/* get polynomial terms in symbolic notation of a homogeneous polynomial of given order 
 * for d dimensions, note: only up to d=4 is implemented here for symbols */
inline void getPolyTerms(char **termval,int termindx,int *termid,int order, int d, int ix)
{
	int i,j;
	char isave[4];
	char varnames[4]={'x','y','z','t'};
	for(j=0;j<termindx;j++) 
		isave[j]=termval[*termid][j];
	
	for(i=ix;i<d;i++) 
	{
		for(j=0;j<termindx;j++) 
			termval[*termid][j]=isave[j];
		
		termval[*termid][termindx]=varnames[i];
	
		if (order > 1) getPolyTerms(termval,termindx+1,termid,order-1,d,i);
		if (order==1) (*termid)++;
	}
}

/* get values of each polynomial term at a given evaluation point, same as evaluating
 * the symbolic polynomial above */
inline void getPolyValues(double *xp, double *termval, int *termid,int order, int d, int ix)
{
	int i,j;
	double isave;
	isave=termval[*termid];
	for(i=ix;i<d;i++) 
	{
		termval[*termid]=isave;
		termval[*termid]*=xp[i];
		if (order > 1) getPolyValues(xp,termval,termid,order-1,d,i);
		if (order == 1) (*termid)++;
	}
 
	return;
}
      


/*
Interpolate using RBF + arbitrary order polynomial
np=number of points
d=number of dimensions
order=order of polynomial
xcloud=cloud of points xcloud[3*np]
P=target point
weights=interpolation weights
*/
inline void d_interprbf_generic(double *xcloud, double *P, scalarList & weights,int np, int d, int order)
{
	int i,j,n,imin,o;
	double **M, *B;
	double dd,dmin;
	int neqns=np+1;
	int nterms;
	int iflag;

	for(i=1;i<=order;i++)
		getPolyCount(&neqns,i,d,0);
  
	M=(double **)malloc(sizeof(double *)*neqns);
	B=(double *)malloc(sizeof(double) *neqns);
	for(i=0;i<neqns;i++) 
		M[i]=(double *)malloc(sizeof(double)*neqns);
  
	dmin=1.0e15;
	imin=0;

	for(i=0;i<np;i++)
    {
		for(j=0;j<np;j++)
		{
			dd=0;
			for(n=0;n<d;n++) 
				dd+=(xcloud[d*i+n]-xcloud[d*j+n])*(xcloud[d*i+n]-xcloud[d*j+n]);
			
			M[i][j]=Foam::exp(-dd);
		}
      
		M[i][np]=M[np][i]=1.0;
		for(j=np+1;j<neqns;j++) 
			M[i][j]=1.0;
      
		nterms=0;
		for(o=1;o<=order;o++) 
			getPolyValues(&(xcloud[d*i]),&(M[i][np+1]),&nterms,o,d,0);
      
		for(j=np+1;j<neqns;j++) 
			M[j][i]=M[i][j];
      
		dd=0;
		for(n=0;n<d;n++) 
			dd+=(xcloud[d*i+n]-P[n])*(xcloud[d*i+n]-P[n]);
      
		if (dd < dmin) 
		{
			imin=i;
			dmin=dd;
		}
		B[i]=Foam::exp(-dd);
    }
	/* lower corner of matrix is zero */
	for(i=np;i<neqns;i++) 
		for(j=np;j<neqns;j++)
			M[i][j]=0.0;
  /* 
  printf("M=\n");
  for(i=0;i<neqns;i++)
   {
    for(j=0;j<neqns;j++)
      printf("%f ",M[i][j]);
    printf("\n");
  }
  */
  
	/* set up the rest of the RHS */
	for(j=np;j<neqns;j++) 
		B[j]=1.0;
  
	nterms=0;
  
	for(o=1;o<=order;o++)
		getPolyValues(P,&(B[np+1]),&nterms,o,d,0);
  /* 
  printf("B=\n");
  for(i=0;i<neqns;i++)
     printf("%f\n",B[i]);
  */
  
	d_solvec(M,B,&iflag,neqns);
  
	if (iflag==0) 
    {
		for(i=0;i<np;i++) 
			weights[i]=0.0;
		
		weights[imin]=1.0;
    }
	else 
    {
		for(i=0;i<np;i++)
		weights[i]=B[i];
    }

	for(i=0;i<neqns;i++) 
		free(M[i]);
  
	free(M);
	free(B);
}
  

  

    inline void getInterpolationWeightsRBF(const dynamicFvMesh & donorMesh, volScalarField & donorBlanks, label RBFPoints)
    {
        forAll(mesh().C(),i)
        {
            if ( donorCells[i] >=0 )
            {
                iweights[i].resize(donorCellCells[i].size());
                iweights[i]=0.0;
                scalar sumWeights=0.0;

				double *xcloud;
				xcloud=new double[donorCellCells[i].size()*2];
				for(int j=0;j<donorCellCells[i].size();j++) 
				{
					xcloud[2*j  ]=donorMesh.C()[donorCellCells[i][j]].x();
					xcloud[2*j+1]=donorMesh.C()[donorCellCells[i][j]].y();
				}
				double P[2]={mesh().C()[i].x(),mesh().C()[i].y()};
				//d_interprbf(xcloud,P,iweights[i],donorCellCells[i].size(),2);
		
				d_interprbf_generic(xcloud, P, iweights[i], RBFPoints, 2, 2);
		
				delete [] xcloud;
			}
		}
    }
	
	
    inline void finalizeClassification(const dynamicFvMesh & donorMesh, const overlap* donorConnect, const volScalarField & cellmanager)
    {
        // Check if interpolation cells's donors are themselves interpolation cells
        forAll(mesh().C(),i)
        {   
            bool flip=false;
            
            // Is a possible interpolation cell (exclude mandatory interpolation cells)
            if ( donorCells[i] >=0 && cellmanager.primitiveField()[i] > CELLMANAGERINTERPCOMP)
            {
                //Check the primary donor;
                // If donor cell is an interpolation cell
                 for (label j = 0 ; j < donorCellCells[i].size(); j++ )
                 {
                    if ( donorConnect->donorCells[donorCellCells[i][j]] >=0)
                    {
                        //flip=true;
                    }
                 }
                 
                 // Also check if the donor has smaller/higher volume?
                 label primaryDonor = donorCells[i];
                 if (donorMesh.V()[primaryDonor] > mesh().V()[i])
                 {
                    flip = true;
                 }
                 
                 
                 if ( flip)
                 {
                    donorCells[i]=-1;
                    donorCellCells[i]=-1;   
                 }
                 
            }
            
            
            if ( cellmanager.primitiveField()[i] < CELLMANAGERINTERPCOMP && donorCells[i] < 0 )
            {
                Info << " Cell " << i <<" is a mandatory interpolation cell but is classified as a regular cell " << endl;
                
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
        
        nInterpolationCells=0;
        forAll(mesh().C(),i)
        {
            nInterpolationCells+=(cellBlank[i]==0);
        }
        
        oversetInterface=oint.toc();
        
    }
    
    inline label locateSingleOversetInterfaceCell(volScalarField & cellBlank)
    {
        if ( cellBlank[mesh().owner()[oversetInterface[0]]] == 1 )
        {
            return mesh().neighbour()[oversetInterface[0]];
        }
        else
            return mesh().owner()[oversetInterface[0]];
    }
    
    
    
};


class overlapSolver
{
    public:
        
        
    overlapSolver(PtrList<dynamicFvMesh> & componentMesh_, PtrList<volScalarField>& cellManager_, PtrList<volScalarField>& cellBlank_, label alphaMesh_)
    :
    componentMesh(componentMesh_),
    cellManager(cellManager_), 
    cellBlank(cellBlank_),
	alphaMesh(alphaMesh_)
    {
        nRegions = componentMesh_.size();
        connect = new overlap* [nRegions];
        //alphaMesh = 0;
        nCellsTotal = 0;
        
        initialize();
    }
    
    ~overlapSolver()
    {
        delete [] connect;   
    }
    
    inline void initialize()
    {
        // First need to flag the cells adjacent to overset boundaries as mandatory
        for( label mi = 0 ; mi < nRegions ; mi++)
        {
            const dynamicFvMesh & mesh = componentMesh[mi];
            
            connect[mi]= new overlap(componentMesh[mi],mi);
            
            nCellsTotal+=componentMesh[mi].nCells();
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
        
        
        Amat.allocate(nCellsTotal,nCellsTotal);
        source.resize(nCellsTotal);
        solution.resize(nCellsTotal);
        stlSolution.reserve(nCellsTotal);
        
    }
    
    inline void updateMovingMesh()
    {
        for ( label i = 0 ; i < nRegions; i++ )
        {
            componentMesh[i].update();
            forAll(cellBlank[i],ci)
            {
                cellBlank[i][ci]=1.0;
                connect[i]->donorCells[ci]=-1;
                connect[i]->donorCellCells[ci].clear();
                connect[i]->iweights[ci].clear();
            }
        }
        
        Amat.destruct();
        Amat.allocate(nCellsTotal,nCellsTotal);
        stlSolution.clear();
        stlSolution.reserve(nCellsTotal);
        source=0.0;
        solution=0.0;
    }
    
    inline void preProcessAndCompute()
    {
		bool useRBF=false;
		if (RBFPoints > 0)
		{
			useRBF=true;
		}
		
        for ( label i = 0 ; i < nRegions; i++ )
        {
            for ( label j = 0 ; j < nRegions ; j++ )
            {
                if ( i!=j)
                {
                    connect[i]->searchFromMeshAndAssign(componentMesh[j]);
                    connect[i]->initialize();
					if (!useRBF)
					{
						connect[i]->donorSearch(componentMesh[j]);
					}
					else
					{
						connect[i]->donorSearch(componentMesh[j],RBFPoints);
					}
                }
            }
        }
		
        
        
        // Donor-interp finalize Classification
        for ( label i = 0 ; i < nRegions; i++ )
        {        
            for ( label j = 0 ; j < nRegions ; j++ )
            {
                if ( i!=j)
                {
                    connect[i]->finalizeClassification(componentMesh[j],connect[j],cellManager[i]);
                }
            }
            
        }
    
        // Set setCellBlanks
        for ( label i = 0 ; i < nRegions; i++ )
        {  
            connect[i]->setCellBlanks(cellBlank[i]);
        }
    
        // Increment overlap if required (nTimes < 3 is default)
        for ( label i = 0 ; i < nRegions; i++ )
        {  
            for (label nTimes = 0 ; nTimes < 3 ; nTimes++)
            {    
				//Info << " Incrementing width by : " << nTimes << endl;
                connect[i]->incrementOverlap(cellBlank[i],cellManager[i]);
            }
        }
    
        // Set the interfaces
        for ( label i = 0 ; i < nRegions; i++ )
        {  
            connect[i]->setOversetInterfaces(cellBlank[i]);
            Info << "Overset Interfaces on Mesh " << i <<"  = " << connect[i]->oversetInterface.size() << "  Interp Cells = " << connect[i]->nInterpolationCells<<endl;
        }
        
    }
    
    inline void distanceWeightedWeights()
    {
        // Interpolation weights
        for ( label i = 0 ; i < nRegions; i++ )
        {  
            for ( label j = 0 ; j < nRegions ; j++ )
            {
                if ( i!=j)
                {
                    connect[i]->getInterpolationWeights(componentMesh[j],cellBlank[j]);
                }
            }
            
        }
    }
    
    inline void polynomialWeights()
    {
       // Interpolation weights
        for ( label i = 0 ; i < nRegions; i++ )
        {  
            for ( label j = 0 ; j < nRegions ; j++ )
            {
                if ( i!=j)
                {
                    connect[i]->getInterpolationWeightsLLS(componentMesh[j],cellBlank[j]);
                }
            }
            
        } 
    }
	
	inline void rbfWeights()
    {
       // Interpolation weights
        for ( label i = 0 ; i < nRegions; i++ )
        {  
            for ( label j = 0 ; j < nRegions ; j++ )
            {
                if ( i!=j)
                {
                    connect[i]->getInterpolationWeightsRBF(componentMesh[j],cellBlank[j], RBFPoints);
                }
            }
            
        } 
    }
    
    
    // Number of meshes
    label nRegions;   
    
    // Number of RBF support cells
	label RBFPoints;
    
    // Total cell count
    label nCellsTotal;
    
    // The Sparse Matrix holding the full discretization
    SM::SparseMatrix<double> Amat;
    
    // The RHS of the discretized Poisson equation
    scalarField source;
    
    // Solution to the Poisson equation
    std::vector<double> stlSolution;
    
    // Solution to the Poisson equation
    scalarField solution;
    
    // The overlap connect object
    overlap** connect;
    
    // Reference to all meshes
    PtrList<dynamicFvMesh>& componentMesh;
    
    // A flag to determine mandatory overlap cells
    PtrList<volScalarField>& cellManager;
    
    // Cell blanking (0 : interp, 1: solve)
    PtrList<volScalarField>& cellBlank;
	
	// Mesh to fix conservation
    label alphaMesh;
        
    
};
