//---------------------------------------------------------
/**    @file	ConnollySurface.h
*     @brief	ConnollySurface.h is the header for CLASS
*               ConnollySurface.cpp								*/
//---------------------------------------------------------

#ifndef ConnollySurface_h
#define ConnollySurface_h

#include "Surface.h"
#include "SurfaceFactory.h"

#ifdef DBGMEM_CRT
	#define _CRTDBG_MAP_ALLOC
	#define _CRTDBG_MAP_ALLOC_NEW
#endif


//#define DEBUG_CONNOLLY

#define MAX_INCIDENT_PROBES 50

#define DEFAULT_PROBE_RADIUS 1.4 // default probe radius


class ConnollyCell{
public:	
	int patch_type;

	virtual ~ConnollyCell()
	{}
};

class FacetCell: public ConnollyCell
{
public:
	double center[3];
	int id[3]; // atoms triplet indices
	bool isSelfIntersecting;	
	double planes[4][4]; // trimming tethraedron
	int plane_indices[3][2]; // every plane depends on a pair of atoms
	vector<double*> self_intersection_planes; // additional self intersection planes
	FacetCell* mirrorCell;	

	FacetCell()
	{
		mirrorCell=NULL;
	}

	virtual ~FacetCell()
	{
		for (unsigned int i=0;i<self_intersection_planes.size();i++)
			deleteVector<double>(self_intersection_planes[i]);
	}
};

class EdgeCell: public ConnollyCell
{
public:
	double center[3];
	int id[2]; // atom pair indices
	bool isSelfIntersecting;
	double clipping_center[3]; // center of the clipping sphere of the torus
	double clipping_radius;
	double major_radius;
	double u[3],v[3];
	double self_intersection_radius; // radius of the clipping sphere in case of spindle torus
	double Rot[3][3]; // rotation matrix of the torus
	double invrot[3][3]; // inverse rotation matrix
	double cutting_planes[2][4]; // two cutting planes of the torus
	vector<double*> additional_planes; // additional clipping planes due to tori clipping by singular facets
	vector<bool> flags; // acute / non acute flags for additional planes	
	bool acute; // if angle between planes is acute perform "and" side test between planes, otherwise "or" side test is sufficient
	double rcoi;

	virtual ~EdgeCell()
	{
		for (unsigned int i=0;i<additional_planes.size();i++)
			deleteVector<double>(additional_planes[i]);
	}
};


class PointCell: public ConnollyCell
{
public:
	int id; // atom index
	vector<EdgeCell*> neighbours;
	vector<FacetCell*> incidentProbes;
	vector<EdgeCell*> buried_neighbours;

	virtual ~PointCell()
	{
		for (unsigned int i=0;i<buried_neighbours.size();i++)
			delete buried_neighbours[i];
	}

};

// row major
// 3d map
#define GRID_CONNOLLY_CELL_MAP(i,j,k,l,NX,NY,NZ) gridConnollyCellMap[(l)+(MAX_CONNOLLY_CELLS-1)*((k)+(NZ)*((j)+(NY)*(i)))]
// 2d map
#define GRID_CONNOLLY_CELL_MAP_2D(i,j,k,NA,NB) gridConnollyCellMap2D[(k)+(MAX_CONNOLLY_CELLS_2D-1)*((j)+(NB)*(i))]

#define SELF_MAP(i,j,k,l,NX,NY,NZ) gridProbesMap[(l)+(MAX_PROBES-1)*((k)+(NZ)*((j)+(NY)*(i)))]

#ifdef ENABLE_CGAL 
//////////////////////// CGAL ///////////////////////////////////////////////////
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_filtered_traits_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Fixed_alpha_shape_3.h>
#include <CGAL/Fixed_alpha_shape_vertex_base_3.h>
#include <CGAL/Fixed_alpha_shape_cell_base_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Ray_3.h>
#include <CGAL/Line_3.h>
#include <CGAL/intersections.h>
#include <cassert>
#include <vector>
#include <fstream>
////////////////////////////////////////////////////////////////////////////////
#endif

// here singular/regular is in the alpha shape nomenclature

#define POINT_CELL				0
#define REGULAR_EDGE_CELL		1
#define SINGULAR_EDGE_CELL		2
#define REGULAR_FACE_CELL		3
#define SINGULAR_FACE_CELL		4

/** @brief This class builds and converts to a DelPhi suitable representation the Connolly-Richards Surface.
All the gathered info is analytically computed both the intersections and the projections. First the
alpha shape of the set of atoms is computed, then from that the Connolly surface is computed by build all
the patches.

@author Sergio Decherchi 
@date 15/05/2012
*/
class ConnollySurface: public Surface
{

#ifdef ENABLE_CGAL 
private:
	typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
	typedef CGAL::Regular_triangulation_filtered_traits_3<K>	Traits;
	typedef CGAL::Polyhedron_3<K>				                Polyhedron;	
	typedef Polyhedron::Halfedge_around_facet_circulator        HF_circulator;
	typedef CGAL::Vector_3<K>									Vector3;
	typedef CGAL::Plane_3<K>									Plane3;
	typedef CGAL::Ray_3<K>									    Ray3;
	typedef CGAL::Line_3<K>									    Line3;	
	typedef CGAL::Triangulation_vertex_base_with_info_3<PointCell*,Traits> Vertex_with_info;		
	typedef CGAL::Fixed_alpha_shape_vertex_base_3<Traits,Vertex_with_info>		Vb;			
	typedef CGAL::Triangulation_cell_base_with_info_3< FacetCell*[4],Traits> Cell_with_info;
	typedef CGAL::Regular_triangulation_cell_base_3<Traits,Cell_with_info> RT_Cell_with_info;
	typedef CGAL::Fixed_alpha_shape_cell_base_3<Traits,RT_Cell_with_info>         Fb;
	typedef CGAL::Triangulation_data_structure_3<Vb,Fb>			Tds;
	typedef Tds::Cell_circulator								Cell_circulator;
	typedef CGAL::Regular_triangulation_3<Traits,Tds>			Rt;
	typedef CGAL::Fixed_alpha_shape_3<Rt>				        Fixed_alpha_shape_3;
	typedef Fixed_alpha_shape_3::Cell_handle					Alpha_Cell_handle;
	typedef Fixed_alpha_shape_3::Vertex_handle					Alpha_Vertex_handle;
	typedef Fixed_alpha_shape_3::Facet							Alpha_Facet;
	typedef Fixed_alpha_shape_3::Edge							Alpha_Edge;
	typedef Traits::Bare_point                                  Point3;
	typedef Traits::Weighted_point                              Weighted_point;

	typedef Rt::Vertex_iterator                                 Vertex_iterator;
	typedef Rt::Finite_vertices_iterator                        Finite_Vertex_Iterator;
	typedef Rt::Finite_cells_iterator							Finite_Cells_Iterator;
	typedef Rt::Finite_edges_iterator							Finite_Edges_Iterator;
	typedef Rt::Finite_facets_iterator							Finite_Facets_Iterator;
	typedef Rt::Vertex_handle                                   Vertex_handle;
	typedef Rt::Facet											Facet;
	

#endif

private:
	/** number of Connolly cells for type*/
	int type[5];
	/** link each auxuliary grid box to the respective cell list */
	int* gridConnollyCellMap2D;
	int* gridConnollyCellMap;	
	/** auxiliary grid sizes*/
	unsigned int nx,ny,nz;
	double scale;
	double side;
	/** auxiliary grid sizes for the 2d map*/
	unsigned int nx_2d,ny_2d,nz_2d;
	double scale_2d;
	double side_2d;
	double xmin_2d,xmax_2d,ymin_2d,ymax_2d,zmin_2d,zmax_2d;
	unsigned int** ind_2d;

	double*x,*y,*z;
	double s;
	double xmin,xmax,ymin,ymax,zmin,zmax;
	unsigned short*** ind;

	unsigned int AUX_GRID_DIM_CONNOLLY;
	unsigned int MAX_CONNOLLY_CELLS;
	unsigned int AUX_GRID_DIM_CONNOLLY_2D;
	unsigned int MAX_CONNOLLY_CELLS_2D;

	/** for each cell there is a structure that defines the patch*/
	vector<ConnollyCell*> sesComplex; 	
	/** this vector contains NULL if the atom is not exposed and the point cell if the atom contributes to the surface*/
	PointCell** atomPatches;
	/** compute the connolly surface using the CGAL alpha shape module and compute all information
	needed by to ray-trace it*/
	bool buildConnolly();
	/** map each cell to the auxiliary grid*/	
	bool buildAuxiliaryGrid();

	#ifdef ENABLE_CGAL
		/** CGAL build Connolly surface */
		bool buildConnollyCGAL();
	#endif

	bool savePovRay;
	int MAX_PROBES;
	/** self intersections grid perfil. Increase this if you use big probes*/
	double si_perfil;

public:
	/** Default constructor*/
	ConnollySurface();
	/** set DelPhi environment*/
	ConnollySurface(DelPhiShared* ds);			
	/** set configuration and DelPhi environment*/
	ConnollySurface(ConfigFile* cf,DelPhiShared* ds);			

	//////////////////////// INTERFACE MANDATORY METHODS /////////////////////////////////
	/** Compute connolly surface. Call it after load*/
	virtual bool build();
	/** Save it in a simple ASCII format (.ses)*/
	virtual bool save(char* fileName);
	/**Load the surface from a file in .ses format*/
	virtual bool load(char* fileName);
	/** Print number of cells and types*/
	virtual void printSummary();		
	/** Get a projection of a point on the surface. Return projection and normal*/
	virtual bool getProjection(double p[3],double* proj1,double* proj2,
		double* proj3,double* normal1,double* normal2,double* normal3);
	/** Get all the intersections of a ray that goes from P1 to P2 over the surface.
	The interesctions are returned with increasing distance order. 
	the first double in the vector is the t parameter for the intersection of the parametric 
	line and the surface, the double pointer is the normal vector. During ray surface intersection
	the previously built auxiliary grid is used to speed up computations*/
	virtual void getRayIntersection(double p1[3],double p2[3],vector<pair<double,double*> >& intersections,int thdID,bool computeNormals);
	/** function for the constructor without arguments*/
	virtual void init();
	/** functions for the constructor with config file argument*/
	virtual void init(ConfigFile* cf);
	/**function for the denstructor*/
	virtual void clear();
	/** pre-process panel to accelerate ray-tracing*/
	virtual void preProcessPanel();
	virtual void postRayCasting();
	virtual bool preBoundaryProjection();
	/////////////////////////////////////////////////////////////

	void setProbeRadius(double pr) 
	{ 
		double e = 1e-2;
		if (pr>e) 
		{ 
			probe_radius=pr; 
		}
		else
		{
			cout << endl << WARN << "Cannot set " << pr << "<=" << e << ". Setting "<< DEFAULT_PROBE_RADIUS;
			probe_radius = DEFAULT_PROBE_RADIUS;
		}
	}
	
	double getProbeRadius() 
	{
		return probe_radius;
	}

	void setSavePovRay(bool ss) 
	{ 
		savePovRay = ss;
	}
	
	bool getSavePovRay() 
	{
		return savePovRay;
	}

	/**for the 3d grid set the max grid size and the maximal number of patches inside a grid cube*/
	void setAuxGrid(unsigned int dim,unsigned int max)
	{		
		AUX_GRID_DIM_CONNOLLY = dim;
		MAX_CONNOLLY_CELLS = max;
	}

	/**for the 2d grid set the max grid size and the maximal number of patches inside a grid cube.
	The grid cube itself does not exist just a reference, indeed the real quantity is MAX_CONNOLLY_CELLS_2D
	that is the number of patches along the grid tube*/
	void setAuxGrid2D(unsigned int dim,unsigned int max)
	{		
		AUX_GRID_DIM_CONNOLLY_2D = dim;
		MAX_CONNOLLY_CELLS_2D = (max*dim);
	}

	/** set the max number of probes to be tested for self intersections grid.
	The higher the probe radius the higher this parameter should be.*/
	void setMaxProbes(int m)
	{ 	
		if (m<=0)
		{
			cout << endl << WARN << "Cannot set max probes <0";
			return;
		}
		MAX_PROBES = m;
	}

	int getMaxProbes()
	{
		return MAX_PROBES;
	}

	void setSIPerfil(double si)
	{
		si_perfil = si;
	}

	double getSIPerfil()
	{
		return si_perfil;
	}

	virtual ~ConnollySurface();

private:

	bool rayConnollyCellIntersection(double*,double*,ConnollyCell*,double*,double*,double*,double*,int& numInt,int);
	/** gives true if the point is inside the list of planes*/
	bool isFeasible(ConnollyCell* cc,double* point);
	/** project a point to a torus defined by torus_equation*/	
	void projectToTorus(double* y,EdgeCell* ec,double* proj,double* norm,double& dist);
	/** given a point on the torus in EdgeCell, it gives the normal to that point without computing the gradient
	explicitly*/
	void getNormalToTorus(double* y,EdgeCell* ec,double* normal);
	/** given a point y compute the normal on that point. This routine does not check that
	the y point really belongs to the surface, this should be assured by the user. If not
	assured the result is meaningless*/
	void getNormal(double* y,ConnollyCell* cc,double* normal);
	/** save concave sphere patch in Pov-Ray format*/
	void saveConcaveSpherePatch(ofstream& of,FacetCell* fc,int i);

//-------		EXTRA -------------
	// void saveProbeCenters(ofstream&, ofstream&, FacetCell*, const int i) const ;
	void saveProbeCenters(ofstream&, FacetCell*, const int i) const ;
	bool write_regTri(Rt& reg_triang) const;
	bool write_alphaFilt(Fixed_alpha_shape_3& alpha, const int flag = 0) const;
//--------------------------------------------------------


	void saveSphere(ostream& of,double* center,double radius);
	void saveAtomPatch(ofstream& of,PointCell* pc);
	void saveEdgePatch(ofstream& of,EdgeCell* ec,int size,double bigR,double* u,double* v,double* w,bool isComplex);
	/** check the orientation. Assume the planes points toward the visible region of the torus*/
	bool orientation(double* pb_center1,double* pb_center2,double* w1,double* w2);
	/** the aim is to sort probes in clockwise order. As reference the first probe is used*/
	void sortProbes(EdgeCell* ec,FacetCell** fcv,int np,int* sorted);		
	void getCoi(double* torus_center,double rcoi,double** sampledPoints,int numPoints,double* u,double* v);
	/** project a point in 3D to a circle in 3D. Input are the point, the radius,center and the
	plane where circle belongs and the output is the projection and the distance. Assume that
	the normal to the plane is unitary*/
	void projectToCircle(double* point,double radius,double* center,double* plane,double* proj,double& dist);
};


// expand it explicitly because Swig is not able to expand it
static class ConnollySurfaceRegister{ 
	static Surface* createSurface(ConfigFile* conf,DelPhiShared* ds) 
	{ 
		return new ConnollySurface(conf,ds); 
	} 
	public: 
		ConnollySurfaceRegister() 
		{ 
			surfaceFactory().add("ses",createSurface); 
		} 
} ConnollySurfaceRegisterObject;


//static SurfaceRecorder<ConnollySurface> sesRecorder("ses");

#endif