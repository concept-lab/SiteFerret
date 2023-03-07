	
#include "ConnollySurface.h"

void ConnollySurface::init()
{
	gridConnollyCellMap=NULL;
	gridConnollyCellMap2D=NULL;
	atomPatches = NULL;
	// set environment
	ind_2d = NULL;
	ind = NULL;
	x = NULL;
	y = x;
	z = y;	
	savePovRay = false;
	// default grid settings
	AUX_GRID_DIM_CONNOLLY = 100;
	AUX_GRID_DIM_CONNOLLY_2D = 50;
	// max connolly cells for each acc grid cube
	MAX_CONNOLLY_CELLS = 400;	
	MAX_CONNOLLY_CELLS_2D = (400*AUX_GRID_DIM_CONNOLLY_2D);
	surfType = MOLECULAR_SURFACE;
	providesAnalyticalNormals = true;
	MAX_PROBES = 100;
	si_perfil = 1.5;
}

void ConnollySurface::init(ConfigFile* cf)
{
	unsigned int maxSESDim2D = cf->read<unsigned int>( "Max_ses_patches_auxiliary_grid_2d_size", 50 );
	unsigned int maxSESPatches2D = cf->read<unsigned int>( "Max_ses_patches_per_auxiliary_grid_2d_cell", 400 );
	unsigned int maxSESDim = cf->read<unsigned int>( "Max_ses_patches_auxiliary_grid_size", 100 );
	unsigned int maxSESPatches = cf->read<unsigned int>( "Max_ses_patches_per_auxiliary_grid_cell", 400 );	
	bool savePovRay = cf->read<bool>( "Save_PovRay", false );
	int mp = cf->read<int>( "Max_Probes_Self_Intersections",100);
	double si_perfil = cf->read<double>( "Self_Intersections_Grid_Coefficient",1.5);

	setAuxGrid(maxSESDim,maxSESPatches);
	setAuxGrid2D(maxSESDim2D,maxSESPatches2D);
	inside = 5;		
	setMaxProbes(mp);
	setSIPerfil(si_perfil);
	setSavePovRay(savePovRay);	
}

ConnollySurface::ConnollySurface():Surface()
{
	init();
}

ConnollySurface::ConnollySurface(DelPhiShared* ds):Surface()
{
	init();
	// set environment
	delphi = ds;	
}

ConnollySurface::ConnollySurface(ConfigFile* cf,DelPhiShared* ds):Surface(cf)
{
	init();
	// set environment
	delphi = ds;	
	init(cf);
}

bool ConnollySurface::buildAuxiliaryGrid()
{
	if (sesComplex.size() == 0)
	{	
		cout << endl << WARN << "Cannot get surface with an empty SES complex!";
		return false;
	}	

	// Perform pre-processing to speed-up intersections and projections.
	// Compute bounding box for each connoly cell in the surface, and map each
	// bounding box to the proper grid point.
	// This routine uses an auxiliary grid. Delphi grid is not directly used
	// because it may be too memory consuming. To speed up computations
	// a maximal number of connolly cells is allowed in each auxiliary grid cell
	// the macro is MAX_CONNOLLY_CELLS in ConnollySurface.h.
	// The macro AUX_GRID_DIM_CONNOLLY sets the maximally allowed grid size.
	unsigned int igrid=delphi->nx;	
	// cannot have an auxiliary grid smaller than that of delphi.
	// all the rest of the code is based on this assumption

	// auxiliary grid is consistent with delphi grid in order to speed-up tracing
	int gridMul=1;
	while(igrid>AUX_GRID_DIM_CONNOLLY)
	{
		gridMul+=2;
		int digrid = delphi->nx;
		while(1)
		{
			// get nearest odd multiple
			int fixedPoint = (digrid+(gridMul-1))/gridMul;
			igrid = fixedPoint*gridMul;
			if (igrid%2==0)
				digrid=igrid+1;
			else
			{
				igrid = fixedPoint;
				break;
			}
		}
	}
	// auxiliary scale
	scale = delphi->scale/((double)gridMul);

	cout << endl << INFO << "Auxiliary grid is " << igrid;
	
	xmin = delphi->baricenter[0]-(igrid-1)/(2*scale);
	ymin = delphi->baricenter[1]-(igrid-1)/(2*scale);
	zmin = delphi->baricenter[2]-(igrid-1)/(2*scale);

	xmax = delphi->baricenter[0]+(igrid-1)/(2*scale);
	ymax = delphi->baricenter[1]+(igrid-1)/(2*scale);
	zmax = delphi->baricenter[2]+(igrid-1)/(2*scale);

	////////////// Allocate memory for the grid, planes and maps //////////////////////

	if (x!=NULL)
		deleteVector<double>(x);
	if (y!=NULL)
		deleteVector<double>(y);
	if (z!=NULL)
		deleteVector<double>(z);

	x = allocateVector<double>(igrid);
	y = allocateVector<double>(igrid);
	z = allocateVector<double>(igrid);

	// delete before the new dimensions are computed
	if (ind!=NULL)
		deleteMatrix3D<unsigned short>(nx,ny,ind);

	side = 1/scale;
	nx = igrid;
	ny = igrid;
	nz = igrid;

	cout << endl << INFO << "Allocating " << (nx*ny*nz*MAX_CONNOLLY_CELLS)*sizeof(int)/1024.0/1024.0 << " MB" << " for the auxiliary grid...";

	if (gridConnollyCellMap!=NULL)
		deleteVector<int>(gridConnollyCellMap);	
	gridConnollyCellMap = allocateVector<int>(nx*ny*nz*MAX_CONNOLLY_CELLS);

	// look up, it has been deallocated!
	ind = allocateMatrix3D<unsigned short>(nx,ny,nz);

	for (unsigned int i=0;i<nx;i++)
		for (unsigned int j=0;j<ny;j++)
			for (unsigned int k=0;k<nz;k++)
				ind[i][j][k]=0;

	for (unsigned int i=0;i<nx;i++)
		x[i]=xmin+i*side;
	
	for (unsigned int i=0;i<ny;i++)
		y[i]=ymin+i*side;

	for (unsigned int i=0;i<nz;i++)
		z[i]=zmin+i*side;
	
	cout << "ok!";

	// build a bounding box for cell and map it to
	// the auxiliary grid
	int max_t = 0;

	cout << endl << INFO << "Mapping auxiliary grid...";
	
	for (unsigned int it=0;it<sesComplex.size();it++)
	{
		// map connolly cell
		ConnollyCell* cc = sesComplex[it];
		 
		// compute bounding box object
		double downx = INFINITY;
		double downy = -INFINITY;
		double downz = INFINITY;

		double upx = -INFINITY;
		double upy = INFINITY;
		double upz = -INFINITY;
		
		double* sphere_center;
		double radius;

		if (cc->patch_type==REGULAR_FACE_CELL || cc->patch_type==SINGULAR_FACE_CELL)
		{
			FacetCell* fc = (FacetCell*)cc;
			sphere_center = fc->center;		
			radius = probe_radius;			
		}
		else if (cc->patch_type==SINGULAR_EDGE_CELL || cc->patch_type==REGULAR_EDGE_CELL)
		{
			EdgeCell* ec = (EdgeCell*)cc;
			sphere_center = ec->clipping_center;
			radius = ec->clipping_radius;
		}
		else if (cc->patch_type==POINT_CELL)
		{
			PointCell* pc = (PointCell*)cc;
			radius = delphi->atoms[pc->id]->radius;
			sphere_center = delphi->atoms[pc->id]->pos;
		}


		downx = sphere_center[0]-radius;
		downy = sphere_center[1]+radius;
		downz = sphere_center[2]-radius;

		upx = sphere_center[0]+radius;
		upy = sphere_center[1]-radius;
		upz = sphere_center[2]+radius;
	          
		// resolve which are the grid cubes that
		// are cut by the bounding box object. These
		// cubes see the cell bounding box    
		int ix_start = (int)rintp((downx-xmin)/side);
		int iy_start = (int)rintp((downy-ymin)/side);
		int iz_start = (int)rintp((downz-zmin)/side);
	     	     
		int ix_end = (int)rintp((upx-xmin)/side);
		int iy_end = (int)rintp((upy-ymin)/side);
		int iz_end = (int)rintp((upz-zmin)/side);

		if (ix_start<0)
		 ix_start=0;
		if (iy_start<0)
		 iy_start=0;
		if (iz_start<0)
		 iz_start=0;

		if (ix_end<0)
		 ix_end=0;
		if (iy_end<0)
		 iy_end=0;
		if (iz_end<0)
		 iz_end=0;

		if (ix_end>=(int)nx)
		 ix_end=nx-1;
		if (iy_end>=(int)ny)
		 iy_end=ny-1;	          
		if (iz_end>=(int)nz)
		 iz_end=nz-1;	          

		if (ix_start>=(int)nx)
		 ix_start=nx-1;
		if (iy_start>=(int)ny)
		 iy_start=ny-1;	          
		if (iz_start>=(int)nz)
		 iz_start=nz-1;	          

		//printf("\n in ");
		for (int iz=iz_start;iz<=iz_end;iz++)
		 for (int iy=iy_start;iy>=iy_end;iy--)
			 for (int ix=ix_start;ix<=ix_end;ix++)
			 {
				 //printf("-");
				 if (ind[ix][iy][iz]>max_t)
					 max_t =ind[ix][iy][iz];

				 if (ind[ix][iy][iz]>=MAX_CONNOLLY_CELLS)
				 {
					cout << endl << ERR << "Number of connolly cells is superior to maximum allowed, please increase Max_ses_patches_per_auxiliary_grid_2d_cell";						
					exit(-1);
				 }
				 GRID_CONNOLLY_CELL_MAP(ix,iy,iz,(ind[ix][iy][iz]),nx,ny,nz) = it;
				 ind[ix][iy][iz]++;
			 }				 

	 //printf(" out");
	}    	

	cout << "ok!";
	cout << endl << INFO << "Max Connolly cells per auxiliary cell -> " << max_t;

	return true;
}
bool ConnollySurface::build()
{
	bool flag;
	#ifdef ENABLE_CGAL 
		flag = buildConnollyCGAL();
	#else
		cout << endl << WARN << "Connolly surface requires CGAL lib";
		return false;
	#endif

    if (!flag)
	{
		cout << endl << ERR << "Connolly surface construction failed";
		return flag;
	}
	return flag;
}

void ConnollySurface::postRayCasting()
{
	// remove 2d grid for ray casting 
	if (gridConnollyCellMap2D!=NULL)
		deleteVector<int>(gridConnollyCellMap2D);	
	
	if (ind_2d!=NULL)
		deleteMatrix2D<unsigned int>(last_rows_ind,ind_2d);

	// no more needed
	//gridConnollyCellMap2D = NULL;
	//ind_2d = NULL;
}


bool ConnollySurface::preBoundaryProjection()
{
	bool flag;
	// useful to project bgp
	if (projBGP)
	{
		flag = buildAuxiliaryGrid();
		if (!flag)
		{
			cout << endl << ERR << "Cannot build 3D auxiliary grid";
			return flag;
		}
	}	
	return true;
}

bool ConnollySurface::buildConnollyCGAL()
{
	// load atoms
	list<Weighted_point> l;
	char outputFile[BUFLEN];	    
	strcpy(outputFile,"connolly.pov");
	double probe_radius2 = probe_radius*probe_radius;

	double x,y,z,r;
	double max_x=-1e6,min_x=1e6,max_y=-1e6,min_y=1e6,max_z=-1e6,min_z=1e6;
			
	for (int i=0;i<delphi->numAtoms;i++)
	{
		delphi->atoms[i]->pos[0] = delphi->atoms[i]->pos[0]+randDisplacement*(randnum()-0.5);
		delphi->atoms[i]->pos[1] = delphi->atoms[i]->pos[1]+randDisplacement*(randnum()-0.5);
		delphi->atoms[i]->pos[2] = delphi->atoms[i]->pos[2]+randDisplacement*(randnum()-0.5);

		x = delphi->atoms[i]->pos[0];
		y = delphi->atoms[i]->pos[1];
		z = delphi->atoms[i]->pos[2];
		r = delphi->atoms[i]->radius;
		r +=probe_radius;
				
		l.push_front(Weighted_point(Point3(x,y,z),(r*r),i)); 

		r-=probe_radius;
		max_x = max(max_x,x+r);
		max_y = max(max_y,y+r);
		max_z = max(max_z,z+r);

		min_x = std::min(min_x,x-r);
		min_y = std::min(min_y,y-r);
		min_z = std::min(min_z,z-r);	
	}

	ofstream of;
	
	if (savePovRay)
		of.open(outputFile);	

	double mid_x = (max_x+min_x)/2.f;
	double mid_y = (max_y+min_y)/2.f;
	double mid_z = (max_z+min_z)/2.f;

	min_x -= fabs(mid_x-min_x)*2;
	min_y -= fabs(mid_y-min_y)*2;
	min_z -= fabs(mid_z-min_z)*2;

	max_x += fabs(mid_x-max_x)*2;
	max_y += fabs(mid_y-max_y)*2;
	max_z += fabs(mid_z-max_z)*2;
	
	if (savePovRay)
	{
		of << "#include \"shapes.inc\" ";
		of << "\n#include \"colors.inc\" ";
		of << "\nglobal_settings {max_trace_level 3}";
		of << "\nbackground { color Black }"; 
		of << "\ncamera {";
		of << "\n\tlocation  <" << mid_x << ","  << max_y <<"," << mid_z <<">";
		of << "\nlook_at  <" << mid_x << "," << mid_y << "," << mid_z << "> translate-<" << mid_x << "," << mid_y << "," << mid_z << ">" << " rotate<0,0,clock> "<< "translate<" << mid_x << "," << mid_y << "," << mid_z << ">}";		
		of << "\nlight_source {<" << mid_x << ","  << max_y <<"," << mid_z <<">" << " color White "<< "translate-<" << mid_x << "," << mid_y << "," << mid_z << ">" << " rotate<0,0,clock> "<< "translate<" << mid_x << "," << mid_y << "," << mid_z << ">}";		
	}

	// setup self intersections grid of the probe	
	double gside = 2*probe_radius;	
	double gscale = 1./gside;

	unsigned int ggrid = 1;
	
	cout << endl << INFO << "Adjusting self intersection grid ";
	while(1)
	{
		ggrid = (unsigned int)floor(gscale*si_perfil*delphi->rmaxdim);
		//cout << ggrid << ",";
		
		if (ggrid<=100)
		{		
			if (ggrid<=0)
			{
				gside=gside/2.;
				gscale = 1./gside;
			}
			else
				break;
		}
		else
		{
			gside=gside*2;
			gscale = 1./gside;
		}
	}	

	double gxmin = delphi->baricenter[0]-(ggrid-1)/(2*gscale);
	double gymin = delphi->baricenter[1]-(ggrid-1)/(2*gscale);
	double gzmin = delphi->baricenter[2]-(ggrid-1)/(2*gscale);

	cout << endl << INFO << "Allocating self intersection grid....";
	FacetCell** gridProbesMap = allocateVector<FacetCell*>(MAX_PROBES*ggrid*ggrid*ggrid);
	for (unsigned int i=0;i<ggrid;i++)
		for (unsigned int j=0;j<ggrid;j++)
			for (unsigned int k=0;k<ggrid;k++)
				// init the number of probes mapped in that place (the pointer is a number)
				SELF_MAP(i,j,k,0,ggrid,ggrid,ggrid)=0;

	cout << "ok!";

	time_t start, end;
	time(&start);	
	
	

	//Fixed_alpha_shape_3  alpha_shape(l.begin(), l.end(), 0.0);
	cout << endl << INFO << "Computing regular triangulation....";
	Rt reg_triang(l.begin(), l.end());
	cout << "ok!";

	// write_regTri(reg_triang);

	cout << endl << INFO << "Computing alpha shape complex....";
	Fixed_alpha_shape_3 alpha_shape(reg_triang, 0.0);
	cout << "ok!";
	//int flag = 1; //REG  = 0; INT = 1; 
	
	/////
	// write_alphaFilt(alpha_shape,1);

	////////////////////////////







	for (int i=0;i<5;i++)
		type[i]=0;

	if (atomPatches!=NULL)
		deleteVector<PointCell*>(atomPatches);

	atomPatches = allocateVector<PointCell*>(delphi->numAtoms);

	for (int i=0;i<delphi->numAtoms;i++)
		atomPatches[i] = NULL;

	sesComplex.reserve(delphi->numAtoms*4);

	vector<int> exposed;	
	////////////////////////////////// atom cells ///////////////////////////////////////////////
	Finite_Vertex_Iterator fvit = alpha_shape.finite_vertices_begin();
	for (;alpha_shape.finite_vertices_end()!=fvit;fvit++)
	{
		int ctype = alpha_shape.classify(fvit);

		//fvit->info();
		// vertex classification is equal to face classification
		// if you go here, the atom is exposed and this a point cell
		if (ctype == Fixed_alpha_shape_3::REGULAR || ctype==Fixed_alpha_shape_3::SINGULAR)
		{
			type[POINT_CELL]++;
			// retrieve atom index
			//cout << endl << fvit->point().index();
		
			PointCell* pc = new PointCell();
			sesComplex.push_back(pc);
			pc->id = fvit->point().index();
			pc->patch_type = POINT_CELL;
			atomPatches[pc->id]= pc;

			exposed.push_back(pc->id);
						
			// get all buried neighbours
			vector<Alpha_Edge> ll;
			alpha_shape.incident_edges(fvit,std::back_inserter(ll));
			for (unsigned int i=0;i<ll.size();i++)
			{
				int ctype = alpha_shape.classify(ll[i]);

				// connection to buried atoms. Get the hidden clipping sphere
				if (ctype == Fixed_alpha_shape_3::INTERIOR)
				{
					const Vertex_handle& v1 = ll[i].first->vertex(ll[i].second);
					const Vertex_handle& v2 = ll[i].first->vertex(ll[i].third);							
					double diff[3];
					diff[0]=v2->point().x()-v1->point().x();
					diff[1]=v2->point().y()-v1->point().y();
					diff[2]=v2->point().z()-v1->point().z();
			
					double dij2 = DOT(diff,diff);
					double dij = sqrt(dij2);
					double ratio =  (v1->point().weight()+dij2-v2->point().weight())/(2*dij);
					double tt2 = v1->point().weight()-ratio*ratio;
					double bigR = sqrt(tt2);
			
					int ind1,ind2;

					ind1 = v1->point().index();
					ind2 = v2->point().index();

					double r1 = delphi->atoms[ind1]->radius;
					double r2 = delphi->atoms[ind2]->radius;

					EdgeCell* ec = new EdgeCell();
					pc->buried_neighbours.push_back(ec);
						
					ec->id[0]=ind1;
					ec->id[1]=ind2;

					ec->major_radius = bigR;

					double rj2[3],cof;
					cof = (ratio/dij);
					VEC_MUL(rj2,diff,cof)
			
					double torus_center[3];
					torus_center[0]=rj2[0]+v1->point().x();
					torus_center[1]=rj2[1]+v1->point().y();
					torus_center[2]=rj2[2]+v1->point().z();
					
					double torus_center_v1_dist2 = DOT(rj2,rj2);

					// 2) get clipping sphere
					// 2.1) get an arbitrary probe position 
					// get by Pitagora theorem the distance from the torus center to the probe
					double torus_center_probe_dist2 = v1->point().weight()-torus_center_v1_dist2;
					double rcoi = sqrt(torus_center_probe_dist2);
					ec->rcoi = rcoi;
					// get an arbitrary probe position
					// 2.2) get the radius and center of the visibility sphere			
					double w[3],ty;
					NORMALIZE_S_ASSIGN(w,diff,ty)

					double d = -DOT(torus_center,w);

					double va[3],u[3],ttt;
					va[0]=0;
					va[1]=w[2];
					va[2]=-w[1];

					NORMALIZE_S_ASSIGN(u,va,ttt)

					// get a random probe position
					double RandomProbe[3];
					ADD_MUL(RandomProbe,torus_center,u,rcoi)
									
					double temp1[3];
					temp1[0]=RandomProbe[0]-v1->point().x();
					temp1[1]=RandomProbe[1]-v1->point().y();
					temp1[2]=RandomProbe[2]-v1->point().z();
			
					double pai = sqrt(DOT(temp1,temp1));
			
					double temp2[3];
					temp2[0]=RandomProbe[0]-v2->point().x();
					temp2[1]=RandomProbe[1]-v2->point().y();
					temp2[2]=RandomProbe[2]-v2->point().z();
			
					double paj = sqrt(DOT(temp2,temp2));
		
					double xx[3];
					xx[0]=temp1[0]/pai*r1;
					xx[1]=temp1[1]/pai*r1;
					xx[2]=temp1[2]/pai*r1;

					double aiaj[3];
					aiaj[0]=v2->point().x()-v1->point().x();
					aiaj[1]=v2->point().y()-v1->point().y();
					aiaj[2]=v2->point().z()-v1->point().z();
			
					double ctemp[3],coeff;
					coeff = pai/(pai+paj);
					VEC_MUL(ctemp,aiaj,coeff)
			
					double clipping_center[3];
					clipping_center[0] = ctemp[0]+v1->point().x();
					clipping_center[1] = ctemp[1]+v1->point().y();
					clipping_center[2] = ctemp[2]+v1->point().z();
			
					double ttemp[3];
					SUB(ttemp,xx,ctemp)
					double r_clipping = sqrt(DOT(ttemp,ttemp));

					ASSIGN(ec->clipping_center,clipping_center)
					ec->clipping_radius = r_clipping;
				}
			}			
			
			//cout << endl << pc->buried_neighbours.size();			
		}
		//else
		//	cout << endl << "buried";
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////

	FILE* fpExp = fopen("exposed.xyz","w");
	fprintf(fpExp,"%zu\n",exposed.size());
	fprintf(fpExp,"Exposed_atoms\n");	
	for (unsigned int i=0;i<exposed.size();i++)
		fprintf(fpExp,"C %f %f %f\n",delphi->atoms[exposed[i]]->pos[0],delphi->atoms[exposed[i]]->pos[1],delphi->atoms[exposed[i]]->pos[2]);
	fclose(fpExp);

	fpExp = fopen("exposedIndices.txt","w");
	for (unsigned int i=0;i<exposed.size();i++)
		fprintf(fpExp,"%d\n",exposed[i]);
	fclose(fpExp);
	
	set<FacetCell*> checkList;

	////////////////////////// facets cells /////////////////////////////////////////////////////////////////

	Finite_Facets_Iterator ffit = alpha_shape.finite_facets_begin();
	for(;ffit!=alpha_shape.finite_facets_end();ffit++)
	{
		int ctype = alpha_shape.classify((*ffit));

		if (ctype==Fixed_alpha_shape_3::REGULAR || ctype==Fixed_alpha_shape_3::SINGULAR)
		{
			bool isRegular = true;
			// probe touches only one side
			if (ctype==Fixed_alpha_shape_3::REGULAR)
				type[REGULAR_FACE_CELL]++;				
			// probe touches both sides
			else if (ctype==Fixed_alpha_shape_3::SINGULAR)
			{
				isRegular = false;
				type[SINGULAR_FACE_CELL]++;
			}
						
			double P1[3],P2[3],P3[3],P12[3],P13[3],P23[3],n[3];

			int index0,index1,index2,index3;
			double w0,w1,w2,w3;

			// get the mirror facet such that the outwards normal is correct
			// depending on the classification of the cell
			if( alpha_shape.classify(ffit->first) == Fixed_alpha_shape_3::INTERIOR)
			{
				const Alpha_Facet& facetp = alpha_shape.mirror_facet((*ffit)); // which gives the same facet viewed from the other adjacent cells 
				// opposite point
				const Weighted_point& p0 = facetp.first->vertex((facetp.second)&3)->point();
				// get points on facet
				const Weighted_point& p1 = facetp.first->vertex((facetp.second+1)&3)->point();
				const Weighted_point& p2 = facetp.first->vertex((facetp.second+2)&3)->point();
				const Weighted_point& p3 = facetp.first->vertex((facetp.second+3)&3)->point();						
				// get the proper normal			
				const Vector3 nn = ( facetp.second % 2 == 1) ? CGAL::normal(p1, p2, p3) : CGAL::normal(p1, p3, p2);
				n[0] = nn.x();
				n[1] = nn.y();
				n[2] = nn.z();

				P1[0] = p1.x();
				P1[1] = p1.y();
				P1[2] = p1.z();

				P2[0] = p2.x();
				P2[1] = p2.y();
				P2[2] = p2.z();

				P3[0] = p3.x();
				P3[1] = p3.y();
				P3[2] = p3.z();

				index0 = p0.index();
				index1 = p1.index();
				index2 = p2.index();
				index3 = p3.index();

				w0 = p0.weight();
				w1 = p1.weight();
				w2 = p2.weight();
				w3 = p3.weight();
			}
			else
			{
				const Alpha_Facet* facetp = &(*ffit);
				// opposite point
				const Weighted_point& p0 = facetp->first->vertex((facetp->second)&3)->point();
				// get points on facet
				const Weighted_point& p1 = facetp->first->vertex((facetp->second+1)&3)->point();
				const Weighted_point& p2 = facetp->first->vertex((facetp->second+2)&3)->point();
				const Weighted_point& p3 = facetp->first->vertex((facetp->second+3)&3)->point();						
				// get the proper normal			
				const Vector3 nn = ( facetp->second % 2 == 1) ? CGAL::normal(p1, p2, p3) : CGAL::normal(p1, p3, p2);
				n[0] = nn.x();
				n[1] = nn.y();
				n[2] = nn.z();

				P1[0] = p1.x();
				P1[1] = p1.y();
				P1[2] = p1.z();

				P2[0] = p2.x();
				P2[1] = p2.y();
				P2[2] = p2.z();

				P3[0] = p3.x();
				P3[1] = p3.y();
				P3[2] = p3.z();

				index0 = p0.index();
				index1 = p1.index();
				index2 = p2.index();
				index3 = p3.index();

				w0 = p0.weight();
				w1 = p1.weight();
				w2 = p2.weight();
				w3 = p3.weight();
			}
								
			SUB(P12,P2,P1)
			SUB(P13,P3,P1)
			SUB(P23,P3,P2)

			// get probe center
			//const double A = P12.squared_length();
			//const double B = P23.squared_length();
			//const double C = P13.squared_length();
			const double A = DOT(P12,P12);
			const double B = DOT(P23,P23);
			const double C = DOT(P13,P13);
			const double D = w1;
			const double E = w2;
			const double F = w3;
			const double X = D+E-A;
			const double Y = E+F-B;
			const double Z = D+F-C;
			const double V = 1./12.*sqrt(4*D*E*F-F*X*X-D*Y*Y-E*Z*Z+X*Y*Z);
			//const double G = P12*P13;
			const double G = DOT(P12,P13);
			//const double H = n*n;
			const double H = DOT(n,n);
			//const double I = P2*P2-P1*P1;
			const double I = (DOT(P2,P2))-(DOT(P1,P1));
			//const double J = P3*P3-P1*P1;
			const double J = (DOT(P3,P3))-(DOT(P1,P1));
			//const double K = n*P1;
			const double K = DOT(n,P1);
			const double L = I-E+D;
			const double M = J-F+D;
			const double S = C*L-G*M;
			const double T = A*M-G*L;

			//Vector3 Q = (K/H)*n+(S/(2*H))*P12+(T/(2*H))*P13;
			double Q[3];

			Q[0] = (K/H)*n[0]+(S/(2*H))*P12[0]+(T/(2*H))*P13[0];
			Q[1] = (K/H)*n[1]+(S/(2*H))*P12[1]+(T/(2*H))*P13[1];
			Q[2] = (K/H)*n[2]+(S/(2*H))*P12[2]+(T/(2*H))*P13[2];

			//Vector3 Probe = Q + ((6*V)/(H))*n;
			double probe[3],coeff;

			coeff = ((6*V)/(H));

			probe[0] = Q[0] + coeff*n[0];
			probe[1] = Q[1] + coeff*n[1];
			probe[2] = Q[2] + coeff*n[2];

			FacetCell* fc1 = new FacetCell();
			if (isRegular)
				fc1->patch_type = REGULAR_FACE_CELL;
			else
				fc1->patch_type = SINGULAR_FACE_CELL;

			sesComplex.push_back(fc1);

			//fc1->center[0]=Probe.x();
			//fc1->center[1]=Probe.y();
			//fc1->center[2]=Probe.z();
			
			fc1->center[0]=probe[0];
			fc1->center[1]=probe[1];
			fc1->center[2]=probe[2];

			fc1->id[0] = index1;
			fc1->id[1] = index2;
			fc1->id[2] = index3;

			atomPatches[fc1->id[0]]->incidentProbes.push_back(fc1);
			atomPatches[fc1->id[1]]->incidentProbes.push_back(fc1);
			atomPatches[fc1->id[2]]->incidentProbes.push_back(fc1);
			
			//Vector3 PP1,PP2,PP3;
			double PP1[3],PP2[3],PP3[3],YY[3];

			// if it is singular this automatically clips the sphere without adding planes
			if (!isRegular)
			{
				// tethraedron points
				//PP1 = P1;
				ASSIGN(PP1,P1);
				//PP2 = P2;
				ASSIGN(PP2,P2);
				//PP3 = P3;
				ASSIGN(PP3,P3);
				fc1->isSelfIntersecting = true;				
			}
			// avoid unecessary self intersection clipping
			else
			{
				// tethraedron points
				//PP1 = Probe + 50000.*(P1-Probe);
				SUB(YY,P1,probe)
				ADD_MUL(PP1,probe,YY,50000.)
				//PP2 = Probe + 50000.*(P2-Probe);
				SUB(YY,P2,probe)
				ADD_MUL(PP2,probe,YY,50000.)
				//PP3 = Probe + 50000.*(P3-Probe);
				SUB(YY,P3,probe)
				ADD_MUL(PP3,probe,YY,50000.)
				fc1->isSelfIntersecting = false;				
			}
			// 0 meaning that we have no additional planes. 
			// If numSelfIntersections==0 and isSelfIntersecting=0 then this means that 
			// the only self intersections is due to mirror facet
			
			//double t1[3],t2[3],t3[3],probe[3];
			double *t1,*t2,*t3;
			t1 = PP1;
			t2 = PP2;
			t3 = PP3;
			/*
			t1[0]=PP1.x();
			t1[1]=PP1.y();
			t1[2]=PP1.z();

			t2[0]=PP2.x();
			t2[1]=PP2.y();
			t2[2]=PP2.z();

			t3[0]=PP3.x();
			t3[1]=PP3.y();
			t3[2]=PP3.z();

			probe[0]=Probe.x();
			probe[1]=Probe.y();
			probe[2]=Probe.z();
			*/

			plane3points(t1,t2,t3,fc1->planes[0]);
			plane3points(t1,t2,probe,fc1->planes[1]);
			plane3points(t2,t3,probe,fc1->planes[2]);
			plane3points(t1,t3,probe,fc1->planes[3]);

			fc1->plane_indices[0][0]=fc1->id[0];
			fc1->plane_indices[0][1]=fc1->id[1];

			fc1->plane_indices[1][0]=fc1->id[1];
			fc1->plane_indices[1][1]=fc1->id[2];

			fc1->plane_indices[2][0]=fc1->id[0];
			fc1->plane_indices[2][1]=fc1->id[2];

			// fix orientation

			// interior reference point
			//Vector3 inner = ((PP1+PP2+PP3)/3.0+Probe)/2.0;
			double inner2[3];
			/*
			inner2[0]=inner.x();
			inner2[1]=inner.y();
			inner2[2]=inner.z();
			*/
			inner2[0] =  ((PP1[0]+PP2[0]+PP3[0])/3.0+probe[0])/2.0;
			inner2[1] =  ((PP1[1]+PP2[1]+PP3[1])/3.0+probe[1])/2.0;
			inner2[2] =  ((PP1[2]+PP2[2]+PP3[2])/3.0+probe[2])/2.0;

			FacetCell* fc2;

			for (int ff=0;ff<4;ff++)
			{
				double orient = DOT(inner2,fc1->planes[ff]);
				// ok
				if (orient<-(fc1->planes[ff][3]))
				{}
				// swap
				else
				{
					fc1->planes[ff][0]=-fc1->planes[ff][0];
					fc1->planes[ff][1]=-fc1->planes[ff][1];
					fc1->planes[ff][2]=-fc1->planes[ff][2];
					fc1->planes[ff][3]=-fc1->planes[ff][3];
				}
			}
															
			if (ctype==Fixed_alpha_shape_3::SINGULAR)
			{
				//Vector3 Probe = Q - ((6*V)/(H))*n;
				double probe[3];
				probe[0] = Q[0] - coeff*n[0];
				probe[1] = Q[1] - coeff*n[1];
				probe[2] = Q[2] - coeff*n[2];

				fc2 = new FacetCell();
				fc2->mirrorCell = fc1;
				fc1->mirrorCell = fc2;
				fc2->patch_type = SINGULAR_FACE_CELL;

				fc2->isSelfIntersecting = true;

				sesComplex.push_back(fc2);

				//fc2->center[0]=Probe.x();
				//fc2->center[1]=Probe.y();
				//fc2->center[2]=Probe.z();

				fc2->center[0]=probe[0];
				fc2->center[1]=probe[1];
				fc2->center[2]=probe[2];
				
				fc2->id[0]= index1;
				fc2->id[1]= index2;
				fc2->id[2]= index3;

				atomPatches[fc2->id[0]]->incidentProbes.push_back(fc2);
				atomPatches[fc2->id[1]]->incidentProbes.push_back(fc2);
				atomPatches[fc2->id[2]]->incidentProbes.push_back(fc2);
								
				//probe[0]=Probe.x();
				//probe[1]=Probe.y();
				//probe[2]=Probe.z();

				plane3points(t1,t2,t3,fc2->planes[0]);
				plane3points(t1,t2,probe,fc2->planes[1]);
				plane3points(t2,t3,probe,fc2->planes[2]);
				plane3points(t1,t3,probe,fc2->planes[3]);

				fc2->plane_indices[0][0]=fc2->id[0];
				fc2->plane_indices[0][1]=fc2->id[1];

				fc2->plane_indices[1][0]=fc2->id[1];
				fc2->plane_indices[1][1]=fc2->id[2];

				fc2->plane_indices[2][0]=fc2->id[0];
				fc2->plane_indices[2][1]=fc2->id[2];

				// interior reference point
				/*
				Vector3 inner = ((PP1+PP2+PP3)/3.0+Probe)/2.0;				
				inner2[0]=inner.x();
				inner2[1]=inner.y();
				inner2[2]=inner.z();*/
				
				inner2[0] =  ((PP1[0]+PP2[0]+PP3[0])/3.0+probe[0])/2.0;
				inner2[1] =  ((PP1[1]+PP2[1]+PP3[1])/3.0+probe[1])/2.0;
				inner2[2] =  ((PP1[2]+PP2[2]+PP3[2])/3.0+probe[2])/2.0;

				for (int ff=0;ff<4;ff++)
				{
					double orient = DOT(inner2,fc2->planes[ff]);
					// ok
					if (orient<-(fc2->planes[ff][3]))
					{}
					// swap
					else
					{
						fc2->planes[ff][0]=-fc2->planes[ff][0];
						fc2->planes[ff][1]=-fc2->planes[ff][1];
						fc2->planes[ff][2]=-fc2->planes[ff][2];
						fc2->planes[ff][3]=-fc2->planes[ff][3];
					}
				}								
			}
			// no mirror patch
			else
				fc1->mirrorCell = NULL;

			// check if self-intersection must be tested or not
			double plane_dist,proj[3],plane[4];
			//double tt1[3],tt2[3],tt3[3];
			double *tt1,*tt2,*tt3;
			
			tt1 = P1;
			tt2 = P2;
			tt3 = P3;
			
			/*
			tt1[0]=P1.x();
			tt1[1]=P1.y();
			tt1[2]=P1.z();

			tt2[0]=P2.x();
			tt2[1]=P2.y();
			tt2[2]=P2.z();

			tt3[0]=P3.x();
			tt3[1]=P3.y();
			tt3[2]=P3.z();
			*/
			plane3points(tt1,tt2,tt3,plane);
			point2plane(fc1->center,plane,&plane_dist,proj);

			if (plane_dist<=probe_radius && ctype!=Fixed_alpha_shape_3::SINGULAR)
			{
				set<FacetCell*>::iterator it = checkList.find(fc1);
				if (it==checkList.end())
				{
					checkList.insert(fc1);
					int ix = (int)rintp((fc1->center[0]-gxmin)/gside);
					int iy = (int)rintp((fc1->center[1]-gymin)/gside);
					int iz = (int)rintp((fc1->center[2]-gzmin)/gside);

					if (ix>=(int)ggrid || iy>=(int)ggrid || iz>=(int)ggrid || ix<=0 || iy<=0 || iz<=0)
					{
						cout << endl << ERR << "Too big probe, to support this probe increase Self_Intersections_Grid_Coefficient\n";
						exit(-1);
					}

					size_t max_ind = (size_t)SELF_MAP(ix,iy,iz,0,ggrid,ggrid,ggrid);
					max_ind++;
					if ((int)max_ind>=MAX_PROBES)
					{
						cout << endl << ERR << "Increase MAX_PROBES";
						exit(-1);
					}
					SELF_MAP(ix,iy,iz,0,ggrid,ggrid,ggrid)=(FacetCell*)max_ind;
					SELF_MAP(ix,iy,iz,max_ind,ggrid,ggrid,ggrid) = fc1;
				}
			}
		}
	}
	
	////////////////////////// edge cells /////////////////////////////////////////////////////////////////
	Finite_Edges_Iterator feit = alpha_shape.finite_edges_begin();
	for(;feit!=alpha_shape.finite_edges_end();feit++)
	{
		int ctype = alpha_shape.classify((*feit));
	
		// exposed edge
		if (ctype==Fixed_alpha_shape_3::REGULAR || ctype==Fixed_alpha_shape_3::SINGULAR)
		{
			// compute major torus radius
			double bigR = 0;
			double smallR = probe_radius;
								
			// up and down edge vertices
			const Vertex_handle& v1 = feit->first->vertex(feit->second);
			const Vertex_handle& v2 = feit->first->vertex(feit->third);

			//Vector3 diff  = v2->point()-v1->point();			
			double diff[3];
			diff[0]=v2->point().x()-v1->point().x();
			diff[1]=v2->point().y()-v1->point().y();
			diff[2]=v2->point().z()-v1->point().z();

			//double dij2 = diff.x()*diff.x()+diff.y()*diff.y()+diff.z()*diff.z();
			double dij2 = DOT(diff,diff);

			double dij = sqrt(dij2);
			double ratio =  (v1->point().weight()+dij2-v2->point().weight())/(2*dij);
			double tt2 = v1->point().weight()-ratio*ratio;
			bigR = sqrt(tt2);
			
			int ind1,ind2;

			ind1 = v1->point().index();
			ind2 = v2->point().index();

			double r1 = delphi->atoms[ind1]->radius;
			double r2 = delphi->atoms[ind2]->radius;

			double ratio2 =  (v2->point().weight()+dij2-v1->point().weight())/(2*dij);
			double bigR2 = sqrt(v2->point().weight()-ratio2*ratio2);
			
			if (fabs(bigR-bigR2)>1e-6)
				cout << endl << WARN << "Non precise torus inner radius " << bigR << " "<< bigR2;

			EdgeCell* ec = new EdgeCell();
			sesComplex.push_back(ec);
						
			ec->id[0]=ind1;
			ec->id[1]=ind2;

			ec->major_radius = bigR;

			// classify the type of the edge
			if ((ctype==Fixed_alpha_shape_3::REGULAR))
			{
				type[REGULAR_EDGE_CELL]++;
				ec->patch_type = REGULAR_EDGE_CELL;
			}
			else
			{				
				type[SINGULAR_EDGE_CELL]++;
				ec->patch_type = SINGULAR_EDGE_CELL;
			}
			
			// compute torus clipping sphere:
			// 1) get torus center			
			double rj2[3],cof;
			cof = (ratio/dij);
			VEC_MUL(rj2,diff,cof)

			double torus_center[3];
			torus_center[0]=rj2[0]+v1->point().x();
			torus_center[1]=rj2[1]+v1->point().y();
			torus_center[2]=rj2[2]+v1->point().z();

			double torus_center_v1_dist2 = DOT(rj2,rj2);

			// 2) get clipping sphere
			// 2.1) get an arbitrary probe position 
			// get by Pitagora theorem the distance from the torus center to the probe
			double torus_center_probe_dist2 = v1->point().weight()-torus_center_v1_dist2;
			double rcoi = sqrt(torus_center_probe_dist2);
			ec->rcoi = rcoi;
			// get an arbitrary probe position
			// 2.2) get the radius and center of the visibility sphere
			double w[3],ty;
			NORMALIZE_S_ASSIGN(w,diff,ty)

			double d = -DOT(torus_center,w);

			// plane of the COI is w*x+d=0. Torus center lies on this plane
			// get a random point on the COI (Circle Of Intersection). The COI can be analytically computed.			
			double cross_temp[3],va[3],vb[3],u[3],v[3],ttt;
			va[0]=0;
			va[1]=w[2];
			va[2]=-w[1];
			vb[0]=w[2];
			vb[1]=0;
			vb[2]=-w[0];

			NORMALIZE_S_ASSIGN(u,va,ttt)
			CROSS(cross_temp,va,vb)
			CROSS(v,cross_temp,va)
			NORMALIZE_S(v,ttt)

			ec->u[0]=u[0];
			ec->u[1]=u[1];
			ec->u[2]=u[2];

			ec->v[0]=v[0];
			ec->v[1]=v[1];
			ec->v[2]=v[2];
			
			// sample coi for debug
			/*double tt = 0;
			for(float tcoi=0.0;tcoi<6.28;tcoi+=0.1)
			{
				// sample the equation of the coi
				Vector3 P = torus_center+rcoi*(cos(tcoi)*U+sin(tcoi)*V);
				char buff2[BUFLEN];				
				sprintf(buff2,"\n\n sphere { \n <%lf,%lf,%lf>,%lf",(P.x()),(P.y()),(P.z()),0.1);	  
				of << buff2;
				of << "\n pigment{ color Yellow}}";	  
			}
			*/

			// get a random probe position
			double RandomProbe[3];
			ADD_MUL(RandomProbe,torus_center,u,rcoi)
						
			double temp1[3];
			temp1[0]=RandomProbe[0]-v1->point().x();
			temp1[1]=RandomProbe[1]-v1->point().y();
			temp1[2]=RandomProbe[2]-v1->point().z();

			double pai = sqrt(DOT(temp1,temp1));

			double temp2[3];
			temp2[0]=RandomProbe[0]-v2->point().x();
			temp2[1]=RandomProbe[1]-v2->point().y();
			temp2[2]=RandomProbe[2]-v2->point().z();

			double paj = sqrt(DOT(temp2,temp2));

			double xx[3];
			xx[0]=temp1[0]/pai*r1;
			xx[1]=temp1[1]/pai*r1;
			xx[2]=temp1[2]/pai*r1;

			double aiaj[3];
			aiaj[0]=v2->point().x()-v1->point().x();
			aiaj[1]=v2->point().y()-v1->point().y();
			aiaj[2]=v2->point().z()-v1->point().z();

			double ctemp[3],coeff;
			coeff = pai/(pai+paj);
			VEC_MUL(ctemp,aiaj,coeff)
			
			double clipping_center[3];
			clipping_center[0] = ctemp[0]+v1->point().x();
			clipping_center[1] = ctemp[1]+v1->point().y();
			clipping_center[2] = ctemp[2]+v1->point().z();


			double ttemp[3];
			SUB(ttemp,xx,ctemp)
			double r_clipping = sqrt(DOT(ttemp,ttemp));
						
			ASSIGN(ec->center,torus_center)

			//  set of incident probes to a given edge
			FacetCell* fcv[MAX_INCIDENT_PROBES];

			// if it is singular no clipping plane is needed
			bool isComplex=false;
			int index = 0;

			// if it is regular needs clipping
			// now get the two triangle patches if the edge is regular
			if (ec->patch_type == REGULAR_EDGE_CELL)
			{
				vector<FacetCell*>& f1 = atomPatches[ec->id[0]]->incidentProbes;
								
				int ref1 = ec->id[0];
				int ref2 = ec->id[1];				
				unsigned int np = (unsigned int)f1.size();
											
				// get all the incident probes. Most of the cases will be 2.
				for (unsigned int indf=0;indf<np;indf++)
				{				
					int i1 = (f1)[indf]->id[0];
					int i2 = (f1)[indf]->id[1];
					int i3 = (f1)[indf]->id[2];

					if	(i1==ref1 && i2==ref2)
					{
						fcv[index]=(f1)[indf];
						index++;
					}
					else if	(i2==ref1 && i1==ref2)
					{
						fcv[index]=(f1)[indf];
						index++;
					}
					else if	(i1==ref1 && i3==ref2)
					{
						fcv[index]=(f1)[indf];
						index++;
					}
					else if	(i3==ref1 && i1==ref2)
					{
						fcv[index]=(f1)[indf];
						index++;
					}
					else if	(i2==ref1 && i3==ref2)
					{
						fcv[index]=(f1)[indf];
						index++;
					}
					else if	(i3==ref1 && i2==ref2)
					{
						fcv[index]=(f1)[indf];
						index++;
					}					
				}
			
				if (fcv[0]==NULL || fcv[1]==NULL)
				{
					cout << endl << ERR << "Error at atoms " << ec->id[0] << "," << ec->id[1];
					cout << endl << ERR << "Regular edge with no probe stations?";
					exit(-1);
				}
				
				// now get the two planes that clip the torus and decide
				// if the clipping is in "AND" or "OR" mode. That is the angle is respectively  <pi or >pi between planes
			
				// two border probes: this is the most common situation
				// recover the two reference planes and deduce if the "AND" or "OR" configuaration is needed				
				if (index==2)
				{
					bool found = false;
					// get first plane
					for (int ii=0;ii<3;ii++)
					{
						int i1 = fcv[0]->plane_indices[ii][0];
						int i2 = fcv[0]->plane_indices[ii][1];

						//printf("\n(1) %d,%d (%d,%d)",i1,i2,ind1,ind2);

						if (((i1==ind1) && (i2==ind2)) || ((i1==ind2) && (i2==ind1)))
						{
							found = true;
							ec->cutting_planes[0][0] = -fcv[0]->planes[ii+1][0];
							ec->cutting_planes[0][1] = -fcv[0]->planes[ii+1][1];
							ec->cutting_planes[0][2] = -fcv[0]->planes[ii+1][2];
							ec->cutting_planes[0][3] = -fcv[0]->planes[ii+1][3];
							//cout << endl << ec->cutting_planes[0][0] << "," << ec->cutting_planes[0][1] << "," << ec->cutting_planes[0][2] << "," << ec->cutting_planes[0][3];
							break;
						}
					}
					
					if (!found)
					{
						cout << endl << ERR << "Cannot detect correct plane!";
						exit(-1);
					}

					found = false;

					// get second plane
					for (int ii=0;ii<3;ii++)
					{
						int i1 = fcv[1]->plane_indices[ii][0];
						int i2 = fcv[1]->plane_indices[ii][1];

						//printf("\n(2) %d,%d (%d,%d)",i1,i2,ind1,ind2);

						if (((i1==ind1) && (i2==ind2)) || ((i1==ind2) && (i2==ind1)))
						{
							found = true;
							ec->cutting_planes[1][0] = -fcv[1]->planes[ii+1][0];
							ec->cutting_planes[1][1] = -fcv[1]->planes[ii+1][1];
							ec->cutting_planes[1][2] = -fcv[1]->planes[ii+1][2];
							ec->cutting_planes[1][3] = -fcv[1]->planes[ii+1][3];
							//cout << endl << ec->cutting_planes[1][0] << "," << ec->cutting_planes[1][1] << "," << ec->cutting_planes[1][2] << "," << ec->cutting_planes[1][3];
							break;
						}
					}
					
					if (!found)
					{
						cout << endl << ERR << "Cannot detect correct plane!";
						exit(-1);
					}

					// get the two reference vectors of the planes.
					double vref1[3],vref2[3];

					SUB(vref1,fcv[0]->center,ec->center)
					SUB(vref2,fcv[1]->center,ec->center)

					// decide AND or TEST based on orientation
					ec->acute = orientation(fcv[0]->center,fcv[1]->center,ec->cutting_planes[0],ec->cutting_planes[1]);					
				}
				// more than 2 incident probes: complicate situation...
				//else 
				if (index>2)
				{										
					// indices of the sorted set of probes
					int sorted[MAX_INCIDENT_PROBES];					

					if ((index%2)!=0)
					{
						cout << endl << ERR << "Torus circulator gives an odd number of probes!";
						exit(-1);
					}

					// manage this situation by explictly trasversing the set of stations of the probes along the current edge
					// sort the probe positions along the COI in order to get a pair of planes every time.
					sortProbes(ec,fcv,index,sorted);

					// as first index find a singular face. Choose the direction by 
					// picking a non singular face
					int start_index=0;
					int direction=0;
					
					// very probably at least one facet is singular: do fast detection of the first probe
					for (int u=0;u<index;u++)
					{
						if (fcv[sorted[u]]->patch_type==SINGULAR_FACE_CELL)
						{
							start_index = u;
							int next = (u+1)%index;
							int prev = (index+u-1)%index;

							if (fcv[sorted[u]]->mirrorCell==fcv[sorted[next]])
								direction=-2;
							else if (fcv[sorted[u]]->mirrorCell==fcv[sorted[prev]])
								direction=+2;
							else
							{
								cout << endl << ERR << "Incosistency in torus circulator";
								exit(-1);
							}

							break;
						}
					}

					// no singular facets, use other, slower criterion						
					if (direction==0)
					{
						isComplex = true;
						double* c1 = fcv[sorted[0]]->center;
						double* c2 = fcv[sorted[1]]->center;
						double mid2[3],plane[3];
						MID(mid2,c1,c2)
						SUB(mid2,mid2,ec->center)

						bool found = false;
						// get refence plane for the test
						for (int ii=0;ii<3;ii++)
						{
							int i1 = fcv[sorted[0]]->plane_indices[ii][0];
							int i2 = fcv[sorted[0]]->plane_indices[ii][1];
							
							if (((i1==ind1) && (i2==ind2)) || ((i1==ind2) && (i2==ind1)))
							{
								found = true;																								
								plane[0] = -fcv[sorted[0]]->planes[ii+1][0];
								plane[1] = -fcv[sorted[0]]->planes[ii+1][1];
								plane[2] = -fcv[sorted[0]]->planes[ii+1][2];								
								//cout << endl << ec->cutting_planes[0][0] << "," << ec->cutting_planes[0][1] << "," << ec->cutting_planes[0][2] << "," << ec->cutting_planes[0][3];
								break;
							}
						}
						if (!found)
						{
							cout << endl <<ERR<< "Cannot identify plane!";
							exit(-1);
						}
					
						double test = -DOT(mid2,plane);
						start_index = 0;
						if (test>0)
							direction=+2;
						else
							direction=-2;
					}

					int currentProbe1=start_index;
					int currentProbe2=(index+start_index+direction/2)%index;

					//cout << endl << "Choosen direction "<<direction;
										
					// get all needed planes
					while(1)
					{						
						FacetCell* fc1 = fcv[sorted[currentProbe1]];
						FacetCell* fc2 = fcv[sorted[currentProbe2]];
						
						bool found = false;
						double* plane1,*plane2;

						// get first plane
						for (int ii=0;ii<3;ii++)
						{
							int i1 = fc1->plane_indices[ii][0];
							int i2 = fc1->plane_indices[ii][1];

							//printf("\n(1) %d,%d (%d,%d)",i1,i2,ind1,ind2);

							if (((i1==ind1) && (i2==ind2)) || ((i1==ind2) && (i2==ind1)))
							{
								found = true;
								double* plane = allocateVector<double>(4);
								plane1  = plane;
								ec->additional_planes.push_back(plane);
								plane[0] = -fc1->planes[ii+1][0];
								plane[1] = -fc1->planes[ii+1][1];
								plane[2] = -fc1->planes[ii+1][2];
								plane[3] = -fc1->planes[ii+1][3];
								//cout << endl << ec->cutting_planes[0][0] << "," << ec->cutting_planes[0][1] << "," << ec->cutting_planes[0][2] << "," << ec->cutting_planes[0][3];
								break;
							}
						}
					
						if (!found)
						{
							cout << endl << ERR << "Cannot detect correct plane!";
							exit(-1);
						}

						found = false;

						// get second plane
						for (int ii=0;ii<3;ii++)
						{
							int i1 = fc2->plane_indices[ii][0];
							int i2 = fc2->plane_indices[ii][1];

							//printf("\n(2) %d,%d (%d,%d)",i1,i2,ind1,ind2);

							if (((i1==ind1) && (i2==ind2)) || ((i1==ind2) && (i2==ind1)))
							{
								found = true;
								double* plane = allocateVector<double>(4);
								plane2 = plane;
								ec->additional_planes.push_back(plane);
								plane[0] = -fc2->planes[ii+1][0];
								plane[1] = -fc2->planes[ii+1][1];
								plane[2] = -fc2->planes[ii+1][2];
								plane[3] = -fc2->planes[ii+1][3];
								//cout << endl << ec->cutting_planes[1][0] << "," << ec->cutting_planes[1][1] << "," << ec->cutting_planes[1][2] << "," << ec->cutting_planes[1][3];
								break;
							}
						}
					
						if (!found)
						{
							cout << endl << ERR << "Cannot detect correct plane!";
							exit(-1);
						}

						// get the two reference vectors of the planes.
						double vref1[3],vref2[3];

						SUB(vref1,fc1->center,ec->center)
						SUB(vref2,fc2->center,ec->center)

						// decide AND or TEST based on orientation
						//bool flag = orientation(vref1,vref2,ec,fc1->center,fc2->center,plane1,plane2);
						bool flag = orientation(fc1->center,fc2->center,plane1,plane2);
						//cout << endl << "Direction " << direction;
						//cout << endl << "Orientation " << flag;						
						ec->flags.push_back(flag);

						// go the next pair using the correct direction
						currentProbe1+=direction;
						currentProbe1=(index+currentProbe1)%index;
						currentProbe2=(index+currentProbe1+direction/2)%index;
						// cycle is restarting, stop
						if (currentProbe1==start_index)
							break;
					}
				}
			}
			// formal "OR". Not executed in practice
			else
				ec->acute=false;

			if ((bigR-probe_radius)>=0)
				ec->isSelfIntersecting = false;
			else
			{
				ec->isSelfIntersecting = true;
				ec->self_intersection_radius = sqrt(-bigR*bigR+probe_radius*probe_radius);

			// if the edge is singular and self intersecting now probes pair must be collected	
			}
			
			ec->clipping_center[0]=clipping_center[0];
			ec->clipping_center[1]=clipping_center[1];
			ec->clipping_center[2]=clipping_center[2];

			ec->clipping_radius = r_clipping;
			
			ec->Rot[0][0]=v[0];
			ec->Rot[0][1]=u[0];
			ec->Rot[0][2]=w[0];
			ec->Rot[1][0]=v[1];
			ec->Rot[1][1]=u[1];
			ec->Rot[1][2]=w[1];
			ec->Rot[2][0]=v[2];
			ec->Rot[2][1]=u[2];
			ec->Rot[2][2]=w[2];

			double det2;
			INVERT_3X3(ec->invrot,det2,ec->Rot)
			
			atomPatches[ec->id[0]]->neighbours.push_back(ec);
			atomPatches[ec->id[1]]->neighbours.push_back(ec);
					
			if (savePovRay)
				saveEdgePatch(of,ec,(int)sesComplex.size(),bigR,u,v,w,isComplex);					
		}		
	}
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	cout << endl << INFO << "Checking " << checkList.size() << " probes for self intersections...";

	// remove self intersections
	for (unsigned int i=0;i<sesComplex.size();i++)
	{
		ConnollyCell* cc = sesComplex[i];
		int type = cc->patch_type;
		if (type==REGULAR_FACE_CELL || type==SINGULAR_FACE_CELL)
		{
			FacetCell* fc1 = (FacetCell*)cc;
			int ix = (int)rintp((fc1->center[0]-gxmin)/gside);
			int iy = (int)rintp((fc1->center[1]-gymin)/gside);
			int iz = (int)rintp((fc1->center[2]-gzmin)/gside);

			//cout << endl << ix << " " << iy << " " << iz;

			for (int k=0;k<27;k++)
			{
				unsigned int cx = ix+shift_map[k][0];
				unsigned int cy = iy+shift_map[k][1];
				unsigned int cz = iz+shift_map[k][2];

				if (cx>(ggrid-1) || cy>(ggrid-1) || cz>(ggrid-1) || cx<0 || cy<0 || cz<0)
					continue;

				size_t max_ind = (size_t)SELF_MAP(cx,cy,cz,0,ggrid,ggrid,ggrid);

				//cout << " mi " << max_ind << endl;

				for (size_t j=0;j<max_ind;j++)
				{
					FacetCell* fc2 = SELF_MAP(cx,cy,cz,j+1,ggrid,ggrid,ggrid);					
					 
					if (fc1==fc2)
						continue;

					double* c1 = fc1->center;
					double* c2 = fc2->center;			
					double dist2;
					DIST2(dist2,c1,c2)
					//cout << endl << dist2 << " " << 4*probe_radius2;
					// there is a self intersection, must clip
					if (dist2<(4*probe_radius2))
					{						
						// a reference point is the mid point
						double ref_point[3],dir[3],temp;
						ADD(ref_point,c2,c1)										
						VEC_MUL(ref_point,ref_point,0.5)					

						fc1->isSelfIntersecting=true;
						fc2->isSelfIntersecting=true;
						
						// get plane direction 
						SUB(dir,c1,c2)
						NORMALIZE_S(dir,temp)
					
						double bias = -DOT(dir,ref_point);					
	
						double* plane = new double[4];
						plane[0] = dir[0]; 
						plane[1] = dir[1];
						plane[2] = dir[2];
						plane[3] = bias;

						double orient = DOT(fc1->center,plane);
						// ok
						if (orient<-(plane[3]))
						{}	
						// swap
						else
						{
							plane[0]=-dir[0];
							plane[1]=-dir[1];
							plane[2]=-dir[2];							
							plane[3]=-bias;
						}
					
						fc1->self_intersection_planes.push_back(plane);
						double* plane2 = new double[4];
						plane2[0]=-plane[0];
						plane2[1]=-plane[1];
						plane2[2]=-plane[2];							
						plane2[3]=-plane[3];
						fc2->self_intersection_planes.push_back(plane2);
					}
				}				
			}
		}
	}
	cout << "ok!";
	
	if (savePovRay)
	{
		cout << endl << INFO << "Saving surface in Pov-Ray in connolly.pov...";
		cout.flush();

		
		for (unsigned int i=0;i<sesComplex.size();i++)
		{
			ConnollyCell* cp = sesComplex[i];
			if (cp->patch_type==REGULAR_FACE_CELL || cp->patch_type==SINGULAR_FACE_CELL)
			{
				FacetCell* fc = (FacetCell*)cp;
				saveConcaveSpherePatch(of,fc,i);
			}
		}

		
	}


	for (unsigned int i=0;i<sesComplex.size();i++)
		{
			ConnollyCell* cp = sesComplex[i];
			if (cp->patch_type==POINT_CELL)
			{
				PointCell* pc = (PointCell*)cp;
				saveAtomPatch(of,pc);
			}
		}
		cout << "ok!";
	

	


		/////////////////////////////////////
//		NEW METHOD FOR TAKING Facet Cell desired information

		std ::ofstream oftxt,ofpqr,oftxtSing,ofpqrSing;
		char outputFile2[1024];
		
		strcpy(outputFile2,"centers.txt");
		oftxt.open(outputFile2);
		// strcpy(outputFile2,"centers.pqr");
		// ofpqr.open(outputFile2);
		// strcpy(outputFile2,"centersSingular.txt");
		// oftxtSing.open(outputFile2);
		// strcpy(outputFile2,"centersSingular.pqr");
		// ofpqrSing.open(outputFile2);

		cout << endl << INFO << "Saving centers of probes related to facet cells";
		//oftxt <<"#Probe radius = "<< probe_radius << endl;
		oftxt << "#Coordinates[3], probe radius \t IDs[3] of neighbouring atoms" << endl;
		for (unsigned int i=0;i<sesComplex.size();i++)
		{
			ConnollyCell* cp = sesComplex[i];

		if (cp->patch_type==REGULAR_FACE_CELL){
			FacetCell* fc = (FacetCell*)cp;
			// saveProbeCenters(oftxt,ofpqr,fc,i);
			saveProbeCenters(oftxt,fc,i);
			}
		// else if (cp->patch_type==SINGULAR_FACE_CELL){
		// 	FacetCell* fc = (FacetCell*)cp;
		// 	saveProbeCenters(oftxtSing,ofpqrSing,fc,i);
		// 	}
		}
		oftxt.close();
		// ofpqr.close();
		// oftxtSing.close();
		// ofpqrSing.close();

		////////////////////////////////



if (savePovRay)
		of.close();

	time(&end);	
	double duration = difftime(end, start);
	cout << endl << INFO << "Surface build-up time.. " << duration << " [s]"; 	

	printSummary();
	// free memory
	deleteVector<FacetCell*>(gridProbesMap);









	return true;

}

bool ConnollySurface::save(char* fileName)
{
	return false;
}

bool ConnollySurface::load(char* fileName)
{
	return false;
}

void ConnollySurface::printSummary()
{
	cout << endl << INFO << "Probe Radius value " << getProbeRadius();
	{
		cout << endl << INFO << "Number of ses cells -> " << sesComplex.size();
		cout << endl << INFO << "Number of del_point cells -> " << type[POINT_CELL];
		cout << endl << INFO << "Number of regular del_edge cells -> " << type[REGULAR_EDGE_CELL];
		cout << endl << INFO << "Number of singular del_edge cells -> " << type[SINGULAR_EDGE_CELL];
		cout << endl << INFO << "Number of regular del_facet cells -> " << type[REGULAR_FACE_CELL];		
		cout << endl << INFO << "Number of singular del_facet cells -> " << type[SINGULAR_FACE_CELL];		

		if (internals!=NULL)
		{
			*internals << endl << "cells " << sesComplex.size();
			*internals << endl << "del_point " << type[POINT_CELL];
			*internals << endl << "rdel_edge " << type[REGULAR_EDGE_CELL];
			*internals << endl << "sdel_edge " << type[SINGULAR_EDGE_CELL];
			*internals << endl << "rdel_facet " << type[REGULAR_FACE_CELL];
			*internals << endl << "sdel_facet " << type[SINGULAR_FACE_CELL];
		}
	}
}

bool ConnollySurface::getProjection(double p[3],double* proj1,double* proj2,
		double* proj3,double* normal1,double* normal2,double* normal3)
{	
	// get the cells that are associated to this grid point
	// by querying the auxiliary grid
	double dist;
	set<int>cells;
	double hside = delphi->hside;

	double aposx = p[0];
	double aposy = p[1];
	double aposz = p[2];

	// move from delphi grid to auxiliary grid
	int irefx = (int)rintp((aposx-xmin)/side);
	int irefy = (int)rintp((aposy-ymin)/side);
	int irefz = (int)rintp((aposz-zmin)/side);

	int ixaux = irefx;
	int iyaux = irefy;
	int izaux = irefz;

	for (int i=0;i<ind[ixaux][iyaux][izaux];i++)
		cells.insert((GRID_CONNOLLY_CELL_MAP(ixaux,iyaux,izaux,i,nx,ny,nz)));

	// keep the nearest patch
	set<int>::iterator it;
	double locProj[3], minDist=INFINITY, locNorm[3];
	bool ff=false;

	dist = 0;
	locProj[0]=0;
	locProj[1]=0;
	locProj[2]=0;

	locNorm[0]=0;
	locNorm[1]=0;
	locNorm[2]=0;

	bool isFeasy = false;
	int ptype;

	for (it=cells.begin();it!=cells.end();it++)
	{	
		ConnollyCell* cc = sesComplex[(*it)];

		if (cc->patch_type==POINT_CELL)
		{
			PointCell* pc = (PointCell*)cc;
			projectToSphere(p,delphi->atoms[pc->id]->pos,delphi->atoms[pc->id]->radius,locProj,dist);
		}
		else if (cc->patch_type==REGULAR_FACE_CELL || cc->patch_type==SINGULAR_FACE_CELL)
		{
			FacetCell* fc = (FacetCell*)cc;
			projectToSphere(p,fc->center,probe_radius,locProj,dist);
		}
		else
		{
			EdgeCell* ec = (EdgeCell*)cc;
			projectToTorus(p,ec,locProj,locNorm,dist);
		}
		
		bool flag = isFeasible(cc,locProj);

		if (!flag)
			continue;

		if (dist<minDist)
		{
			minDist = dist;
			(*proj1)=locProj[0];
			(*proj2)=locProj[1];
			(*proj3)=locProj[2];

			(*normal1)=locNorm[0];
			(*normal2)=locNorm[1];
			(*normal3)=locNorm[2];
			isFeasy = flag;
			ptype = cc->patch_type;
		}
	}

	if (cells.size()==0)
	{
		cout << endl << WARN << "Empty cell in getProjection!";
		return false;
	}

	// Connolly surface is not that nice. Projections can't be analytical everywhere.
	// We are trying to project a point into a singularity? Old DelPhi style.
	if (minDist == INFINITY)
	{
		bool fixed = false;

		double** sampledPoints;									
		int coiNum = 500;
		sampledPoints = allocateMatrix2D<double>(coiNum,3);

		for (it=cells.begin();it!=cells.end();it++)			
		{	
			ConnollyCell* cc = sesComplex[(*it)];					
			if (cc->patch_type==SINGULAR_EDGE_CELL ||  cc->patch_type==REGULAR_EDGE_CELL)
			{
				EdgeCell* ec = (EdgeCell*)cc;
				
				// try to manage singularity. Explicitly project to probes and keep the nearest one
				if (ec->isSelfIntersecting)
				{
					// identify the nearest the common probes and check the feasibility of the projection
					// for at least one of them.
					
					//  set of incident probes to a given edge
					FacetCell* fcv[MAX_INCIDENT_PROBES];					
					int index = 0;

					// manage SINGULAR EDGE?
					if (ec->patch_type == SINGULAR_EDGE_CELL)
					{
						{
						#ifdef ENABLE_BOOST_THREADS
							boost::mutex::scoped_lock scopedLock(mutex);
						#endif
						(*errorStream) << endl << WARN << "Singular edge projection";
						}
					}

					// get all the sourrounding probe stations
					if (ec->patch_type == REGULAR_EDGE_CELL)
					{
						vector<FacetCell*>& f1 = atomPatches[ec->id[0]]->incidentProbes;
										
						int ref1 = ec->id[0];
						int ref2 = ec->id[1];				
						size_t np = (int)f1.size();
																			
						// get all the incident probes. Most of the cases will be 2.
						for (unsigned int indf=0;indf<np;indf++)
						{				
							int i1 = (f1)[indf]->id[0];
							int i2 = (f1)[indf]->id[1];
							int i3 = (f1)[indf]->id[2];

							if	(i1==ref1 && i2==ref2)
							{
								fcv[index]=(f1)[indf];
								index++;
							}
							else if	(i2==ref1 && i1==ref2)
							{
								fcv[index]=(f1)[indf];
								index++;
							}
							else if	(i1==ref1 && i3==ref2)
							{
								fcv[index]=(f1)[indf];
								index++;
							}
							else if	(i3==ref1 && i1==ref2)
							{
								fcv[index]=(f1)[indf];
								index++;
							}
							else if	(i2==ref1 && i3==ref2)
							{
								fcv[index]=(f1)[indf];
								index++;
							}
							else if	(i3==ref1 && i2==ref2)
							{
								fcv[index]=(f1)[indf];
								index++;
							}					
						}
					}	
				
					getCoi(ec->center,ec->rcoi,sampledPoints,coiNum,ec->u,ec->v);
					
					// for all the points get projection. Keep the nearest feasible
					for (int kk=0;kk<coiNum;kk++)
					{
						projectToSphere(p,sampledPoints[kk],probe_radius,locProj,dist);
						for (int ll=0;ll<index;ll++)
						{						
							bool flag = isFeasible(fcv[ll],locProj);
							
							if (!flag)
								continue;

							if (dist<minDist)
							{
								minDist = dist;
								(*proj1)=locProj[0];
								(*proj2)=locProj[1];
								(*proj3)=locProj[2];
								(*normal1)=0;
								(*normal2)=0;
								(*normal3)=0;
								fixed = true;
							}
						}
					}					
				}				
			}			
		} // end singular projection
		deleteMatrix2D<double>(coiNum,sampledPoints);
		
		if (!fixed)
		{
			(*proj1)=p[0];
			(*proj2)=p[1];
			(*proj3)=p[2];

			{  
			#ifdef ENABLE_BOOST_THREADS
				boost::mutex::scoped_lock scopedLock(mutex);
			#endif			
			(*errorStream) << endl << WARN << "Approximating bgp with grid point "; 					
			}
		}
		else
		{
			// non tangent continuity has been properly managed. 
			// the obtained projection is the best thing allowed by Connolly definition.
		}	
	}
	return true;	
}

void ConnollySurface::preProcessPanel()
{
	unsigned int igrid=delphi->nx;	
	// cannot have an auxiliary grid smaller than that of delphi.
	// all the rest of the code is based on this assumption

	// auxiliary grid is consistent with delphi grid in order to speed-up tracing
	int gridMul=1;
	while(igrid>AUX_GRID_DIM_CONNOLLY_2D)
	{
		gridMul+=2;
		int digrid = delphi->nx;
		while(1)
		{
			// get nearest odd multiple
			int fixedPoint = (digrid+(gridMul-1))/gridMul;
			igrid = fixedPoint*gridMul;
			if (igrid%2==0)
				digrid=igrid+1;
			else
			{
				igrid = fixedPoint;
				break;
			}
		}
	}
	// auxiliary scale
	scale_2d = delphi->scale/((double)gridMul);

	//cout << endl << INFO << "Auxiliary grid is " << igrid;
	
	xmin_2d = delphi->baricenter[0]-(igrid-1)/(2*scale_2d);
	ymin_2d = delphi->baricenter[1]-(igrid-1)/(2*scale_2d);
	zmin_2d = delphi->baricenter[2]-(igrid-1)/(2*scale_2d);

	xmax_2d = delphi->baricenter[0]+(igrid-1)/(2*scale_2d);
	ymax_2d = delphi->baricenter[1]+(igrid-1)/(2*scale_2d);
	zmax_2d = delphi->baricenter[2]+(igrid-1)/(2*scale_2d);

	side_2d = 1/scale_2d;
	nx_2d = igrid;
	ny_2d = igrid;
	nz_2d = igrid;

	if (gridConnollyCellMap2D!=NULL)
		deleteVector<int>(gridConnollyCellMap2D);	
	
	if (ind_2d!=NULL)
		deleteMatrix2D<unsigned int>(last_rows_ind,ind_2d);

	if (panel==0)
	{
		last_rows_ind = ny_2d;
		last_cols_ind = nz_2d;
	}
	else if (panel==1)
	{
		last_rows_ind = nx_2d;
		last_cols_ind = ny_2d;
	}
	else
	{
		last_rows_ind = nx_2d;
		last_cols_ind = nz_2d;	
	}
	
	gridConnollyCellMap2D = allocateVector<int>(last_rows_ind*last_cols_ind*MAX_CONNOLLY_CELLS_2D);

	// look up, it has been deallocated
	ind_2d = allocateMatrix2D<unsigned int>(last_rows_ind,last_cols_ind);
	for (int i=0;i<last_rows_ind;i++)
		for (int j=0;j<last_cols_ind;j++)			
			ind_2d[i][j]=0;

	
	// build a bounding box for cell and map it to
	// the auxiliary grid
	unsigned int max_t = 0;

	//cout << endl << INFO << "Mapping auxiliary grid...";
	
	for (unsigned int it=0;it<sesComplex.size();it++)
	{
		// map connolly cell
		ConnollyCell* cc = sesComplex[it];
		 
		// compute bounding box object
		double downx = INFINITY;
		double downy = -INFINITY;
		double downz = INFINITY;

		double upx = -INFINITY;
		double upy = INFINITY;
		double upz = -INFINITY;
		
		double* sphere_center;
		double radius;

		if (cc->patch_type==REGULAR_FACE_CELL || cc->patch_type==SINGULAR_FACE_CELL)
		{
			FacetCell* fc = (FacetCell*)cc;
			sphere_center = fc->center;		
			radius = probe_radius;			
		}
		else if (cc->patch_type==SINGULAR_EDGE_CELL || cc->patch_type==REGULAR_EDGE_CELL)
		{
			EdgeCell* ec = (EdgeCell*)cc;
			sphere_center = ec->clipping_center;
			radius = ec->clipping_radius;
		}
		else if (cc->patch_type==POINT_CELL)
		{
			PointCell* pc = (PointCell*)cc;
			radius = delphi->atoms[pc->id]->radius;
			sphere_center = delphi->atoms[pc->id]->pos;
		}


		downx = sphere_center[0]-radius;
		downy = sphere_center[1]+radius;
		downz = sphere_center[2]-radius;

		upx = sphere_center[0]+radius;
		upy = sphere_center[1]-radius;
		upz = sphere_center[2]+radius;
	          
		// resolve which are the grid cubes that
		// are cut by the bounding box object. These
		// cubes see the cell bounding box    
		int ix_start = (int)rintp((downx-xmin_2d)/side_2d);
		int iy_start = (int)rintp((downy-ymin_2d)/side_2d);
		int iz_start = (int)rintp((downz-zmin_2d)/side_2d);
	     	     
		int ix_end = (int)rintp((upx-xmin_2d)/side_2d);
		int iy_end = (int)rintp((upy-ymin_2d)/side_2d);
		int iz_end = (int)rintp((upz-zmin_2d)/side_2d);

		if (ix_start<0)
		 ix_start=0;
		if (iy_start<0)
		 iy_start=0;
		if (iz_start<0)
		 iz_start=0;

		if (ix_end<0)
		 ix_end=0;
		if (iy_end<0)
		 iy_end=0;
		if (iz_end<0)
		 iz_end=0;

		if (ix_end>=(int)nx_2d)
		 ix_end=nx_2d-1;
		if (iy_end>=(int)ny_2d)
		 iy_end=ny_2d-1;	          
		if (iz_end>=(int)nz_2d)
		 iz_end=nz_2d-1;	          

		if (ix_start>=(int)nx_2d)
		 ix_start=nx_2d-1;
		if (iy_start>=(int)ny_2d)
		 iy_start=ny_2d-1;	          
		if (iz_start>=(int)nz_2d)
		 iz_start=nz_2d-1;	          

		// map the points into the 2D matrices of each plane
		
		if (panel==0)
		{
			// plane yz
			for (int iz=iz_start;iz<=iz_end;iz++)
			{
				for (int iy=iy_start;iy>=iy_end;iy--)			 
				{
					if (ind_2d[iy][iz]>max_t)
						max_t =ind_2d[iy][iz];			
				
					if (ind_2d[iy][iz]>=MAX_CONNOLLY_CELLS_2D)
					{
						cout << endl << ERR << "Number of connolly cells is superior to maximum allowed, please increase Max_ses_patches_per_auxiliary_grid_2d_cell";						
						exit(-1);
					}
					GRID_CONNOLLY_CELL_MAP_2D(iy,iz,(ind_2d[iy][iz]),ny_2d,nz_2d) = it;
					ind_2d[iy][iz]++;
				}
			}
		}

		else if (panel==1)
		{
			// plane xy
			for (int ix=ix_start;ix<=ix_end;ix++)
			{
				for (int iy=iy_start;iy>=iy_end;iy--)			 
				{
					if (ind_2d[ix][iy]>max_t)
						max_t =ind_2d[ix][iy];			
				
					if (ind_2d[ix][iy]>=MAX_CONNOLLY_CELLS_2D)
					{
						cout << endl << ERR << "Number of connolly cells is superior to maximum allowed, please increase Max_ses_patches_per_auxiliary_grid_2d_cell";						
						exit(-1);
					}
					GRID_CONNOLLY_CELL_MAP_2D(ix,iy,(ind_2d[ix][iy]),nx_2d,ny_2d) = it;
					ind_2d[ix][iy]++;
				}
			}
		}
		else
		{
			// plane xz
			for (int ix=ix_start;ix<=ix_end;ix++)
			{
				for (int iz=iz_start;iz<=iz_end;iz++)
				{
					if (ind_2d[ix][iz]>max_t)
						max_t =ind_2d[ix][iz];			
				
					if (ind_2d[ix][iz]>=MAX_CONNOLLY_CELLS_2D)
					{
						cout << endl << ERR << "Number of connolly cells is superior to maximum allowed, please increase Max_ses_patches_per_auxiliary_grid_cell";						
						exit(-1);
					}
					GRID_CONNOLLY_CELL_MAP_2D(ix,iz,(ind_2d[ix][iz]),nx_2d,nz_2d) = it;
					ind_2d[ix][iz]++;
				}
			}
		}
	}    	

	/*
	if (gridLoad!=NULL)
		deleteVector<int>(gridLoad);

	//printf("\n Load balancer...");
	// the load is splitted along the dimension defined by surface class
	// threading strategy thus: ny,ny,nz
	// query the auxiliary on how many are patches encoutered by that ray
	if (panel==0)
	{
		totalLoad = 0;
		gridLoad = allocateVector<int>(delphi->ny);
		for (int i=0;i<delphi->ny;i++)
		{
			gridLoad[i]=0;
			// along z
			for (int j=0;j<delphi->nz;j++)
			{
				int t1 = (int)rintp((delphi->y[i]-ymin_2d)/side_2d);
				int t2 = (int)rintp((delphi->z[j]-zmin_2d)/side_2d);
				gridLoad[i]+=ind_2d[t1][t2];
			}
			totalLoad+=gridLoad[i];
		}
	}
	else if (panel == 1)
	{
		totalLoad = 0;
		gridLoad = allocateVector<int>(delphi->ny);
		for (int i=0;i<delphi->ny;i++)
		{
			gridLoad[i]=0;
			// along x
			for (int j=0;j<delphi->nx;j++)
			{
				int t1 = (int)rintp((delphi->y[i]-ymin_2d)/side_2d);
				int t2 = (int)rintp((delphi->x[j]-xmin_2d)/side_2d);
				gridLoad[i]+=ind_2d[t2][t1];
			}
			totalLoad+=gridLoad[i];
		}
	}
	else 
	{
		totalLoad = 0;
		gridLoad = allocateVector<int>(delphi->nz);
		for (int i=0;i<delphi->nz;i++)
		{
			gridLoad[i]=0;
			// along x
			for (int j=0;j<delphi->nx;j++)
			{
				int t1 = (int)rintp((delphi->z[i]-zmin_2d)/side_2d);
				int t2 = (int)rintp((delphi->x[j]-xmin_2d)/side_2d);
				gridLoad[i]+=ind_2d[t2][t1];
			}
			totalLoad+=gridLoad[i];
		}
	}	
	*/
}

void ConnollySurface::getRayIntersection(double pa[3],double pb[3],vector<pair<double,double*> >& intersections,int thdID,bool computeNormals)
{
	double dir[3];
	double aposx = pa[0],aposy = pa[1],aposz = pa[2];
	dir[0]=pb[0]-aposx;
	dir[1]=pb[1]-aposy;
	dir[2]=pb[2]-aposz;
	double apos[3];
	apos[0]=aposx;
	apos[1]=aposy;
	apos[2]=aposz;		
	int iyaux,izaux,ixaux;
		
	intersections.reserve(2000);

	int numInt = 0;

	if (panel==0)
	{
		int t1 = (int)rintp((pa[1]-ymin_2d)/side_2d);
		int t2 = (int)rintp((pa[2]-zmin_2d)/side_2d);

		unsigned int numcells = ind_2d[t1][t2];		

		iyaux = t1;
		izaux = t2;

		for (unsigned int iter = 0;iter<numcells;iter++)
		{				
			int it = GRID_CONNOLLY_CELL_MAP_2D(iyaux,izaux,iter,ny_2d,nz_2d);
			double t1,t2,t3,t4;
			int ni;
			double intPoint[3];

			bool ff = rayConnollyCellIntersection(pa,dir,sesComplex[it],&t1,&t2,&t3,&t4,ni,thdID);				
			numInt++;
				
			// no intersection
			if (!ff)
				continue;

			intPoint[0]=apos[0]+dir[0]*t1;
			intPoint[1]=apos[1];
			intPoint[2]=apos[2];
			
			// feasibility test: check the intersection is inside the cell
			bool testFeasibility;

			testFeasibility = isFeasible(sesComplex[it],intPoint);

			if(testFeasibility)
			{
				// compute normal								
				if (computeNormals)
				{
					double* n = allocateVector<double>(3);
					getNormal(intPoint,sesComplex[it],n);
					intersections.push_back(pair<double,double*>(t1,n));
				}
				else
					intersections.push_back(pair<double,double*>(t1,(double*)NULL));
			}

			if (ni==1)
				continue;

			intPoint[0]=apos[0]+dir[0]*t2;
			intPoint[1]=apos[1];
			intPoint[2]=apos[2];

			testFeasibility = isFeasible(sesComplex[it],intPoint);
			if(testFeasibility)
			{
				if (computeNormals)
				{
					double* n = allocateVector<double>(3);
					getNormal(intPoint,sesComplex[it],n);
					intersections.push_back(pair<double,double*>(t2,n));
				}
				else
					intersections.push_back(pair<double,double*>(t2,(double*)NULL));
			}
				
			if (ni<3)
				continue;

			intPoint[0]=apos[0]+dir[0]*t3;
			intPoint[1]=apos[1];
			intPoint[2]=apos[2];

			testFeasibility = isFeasible(sesComplex[it],intPoint);
			
			if(testFeasibility)
			{
				// compute normal								
				if (computeNormals)
				{
					double* n = allocateVector<double>(3);
					getNormal(intPoint,sesComplex[it],n);
					intersections.push_back(pair<double,double*>(t3,n));
				}
				else
					intersections.push_back(pair<double,double*>(t3,(double*)NULL));
			}

			if (ni<4)
				continue;

			intPoint[0]=apos[0]+dir[0]*t4;
			intPoint[1]=apos[1];
			intPoint[2]=apos[2];

			testFeasibility = isFeasible(sesComplex[it],intPoint);
			if(testFeasibility)
			{
				// compute normal								
				if (computeNormals)
				{
					double* n = allocateVector<double>(3);
					getNormal(intPoint,sesComplex[it],n);
					intersections.push_back(pair<double,double*>(t4,n));
				}
				else
					intersections.push_back(pair<double,double*>(t4,(double*)NULL));			
			}
		}
	}
	else if(panel==1)
	{
		int t1 = (int)rintp((pa[1]-ymin_2d)/side_2d);
		int t2 = (int)rintp((pa[0]-xmin_2d)/side_2d);

		iyaux = t1;
		ixaux = t2;
				
		unsigned int numcells = ind_2d[ixaux][iyaux];

		for (unsigned int iter = 0;iter<numcells;iter++)
		{				
			int it = GRID_CONNOLLY_CELL_MAP_2D(ixaux,iyaux,iter,nx_2d,ny_2d);
			double t1,t2,t3,t4;
			int ni;
			double intPoint[3];

			// ray intersection
			bool ff;

			ff = rayConnollyCellIntersection(pa,dir,sesComplex[it],&t1,&t2,&t3,&t4,ni,thdID);				

			// no intersection
			if (!ff)
				continue;

			intPoint[0]=apos[0];
			intPoint[1]=apos[1];
			intPoint[2]=apos[2]+dir[2]*t1;

			// feasibility test: check the intersection is inside the cell				
			bool testFeasibility = isFeasible(sesComplex[it],intPoint);

			if(testFeasibility)
			{
				// compute normal								
				if (computeNormals)
				{
					double* n = allocateVector<double>(3);
					getNormal(intPoint,sesComplex[it],n);
					intersections.push_back(pair<double,double*>(t1,n));
				}
				else
					intersections.push_back(pair<double,double*>(t1,(double*)NULL));
			}

			if (ni==1)
				continue;

			intPoint[0]=apos[0];
			intPoint[1]=apos[1];
			intPoint[2]=apos[2]+dir[2]*t2;

			testFeasibility = isFeasible(sesComplex[it],intPoint);
			if(testFeasibility)
			{
				// compute normal								
				if (computeNormals)
				{
					double* n = allocateVector<double>(3);
					getNormal(intPoint,sesComplex[it],n);
					intersections.push_back(pair<double,double*>(t2,n));
				}
				else
					intersections.push_back(pair<double,double*>(t2,(double*)NULL));
			}

			if (ni<3)
				continue;
				
			intPoint[0]=apos[0];
			intPoint[1]=apos[1];
			intPoint[2]=apos[2]+dir[2]*t3;

			testFeasibility = isFeasible(sesComplex[it],intPoint);
			
			if(testFeasibility)
			{
				// compute normal
				if (computeNormals)
				{
					double* n = allocateVector<double>(3);
					getNormal(intPoint,sesComplex[it],n);
					intersections.push_back(pair<double,double*>(t3,n));
				}
				else
					intersections.push_back(pair<double,double*>(t3,(double*)NULL));
			}

			if (ni<4)
				continue;
				
			intPoint[0]=apos[0];
			intPoint[1]=apos[1];
			intPoint[2]=apos[2]+dir[2]*t4;

			testFeasibility = isFeasible(sesComplex[it],intPoint);

			if(testFeasibility)
			{
				// compute normal
				if (computeNormals)
				{
					double* n = allocateVector<double>(3);
					getNormal(intPoint,sesComplex[it],n);
					intersections.push_back(pair<double,double*>(t4,n));
				}
				else
					intersections.push_back(pair<double,double*>(t4,(double*)NULL));
			}
		}	
	}
	else
	{

		int t1 = (int)rintp((pa[0]-xmin_2d)/side_2d);
		int t2 = (int)rintp((pa[2]-zmin_2d)/side_2d);

		ixaux = t1;
		izaux = t2;
		
		unsigned int numcells = ind_2d[ixaux][izaux];

		for (unsigned int iter = 0;iter<numcells;iter++)
		{				
			int it = GRID_CONNOLLY_CELL_MAP_2D(ixaux,izaux,iter,nx_2d,nz_2d);
			double t1,t2,t3,t4;
			int ni;
			double intPoint[3];

			// ray intersection
			bool ff;

			ff = rayConnollyCellIntersection(pa,dir,sesComplex[it],&t1,&t2,&t3,&t4,ni,thdID);				

			numInt++;

			// no intersection
			if (!ff)
				continue;

			intPoint[0]=apos[0];
			intPoint[1]=apos[1]+dir[1]*t1;
			intPoint[2]=apos[2];

			// feasibility test: check the intersection is inside the cell								
			bool testFeasibility;
			testFeasibility = isFeasible(sesComplex[it],intPoint);
				
			// accept the intersection if the intersection is in the cube
			// and the intersection exists 							
			if(testFeasibility)
			{
				// compute normal
				if (computeNormals)
				{
					double* n = allocateVector<double>(3);
					getNormal(intPoint,sesComplex[it],n);
					intersections.push_back(pair<double,double*>(t1,n));
				}
				else
					intersections.push_back(pair<double,double*>(t1,(double*)NULL));
			}

			if (ni==1)
				continue;

			intPoint[0]=apos[0];
			intPoint[1]=apos[1]+dir[1]*t2;
			intPoint[2]=apos[2];

			testFeasibility = isFeasible(sesComplex[it],intPoint);

			if(testFeasibility)
			{
				// compute normal
				if (computeNormals)
				{
					double* n = allocateVector<double>(3);
					getNormal(intPoint,sesComplex[it],n);
					intersections.push_back(pair<double,double*>(t2,n));
				}
				else
					intersections.push_back(pair<double,double*>(t2,(double*)NULL));
			}
				
			if (ni<3)
				continue;

			intPoint[0]=apos[0];
			intPoint[1]=apos[1]+dir[1]*t3;
			intPoint[2]=apos[2];

			testFeasibility = isFeasible(sesComplex[it],intPoint);
			if(testFeasibility)
			{
				// compute normal
				if (computeNormals)
				{
					double* n = allocateVector<double>(3);
					getNormal(intPoint,sesComplex[it],n);
					intersections.push_back(pair<double,double*>(t3,n));
				}
				else
					intersections.push_back(pair<double,double*>(t3,(double*)NULL));
			}
				
			if (ni<4)
				continue;

			intPoint[0]=apos[0];
			intPoint[1]=apos[1]+dir[1]*t4;
			intPoint[2]=apos[2];

			testFeasibility = isFeasible(sesComplex[it],intPoint);
			
			if(testFeasibility)
			{
				// compute normal
				if (computeNormals)
				{
					double* n = allocateVector<double>(3);
					getNormal(intPoint,sesComplex[it],n);
					intersections.push_back(pair<double,double*>(t4,n));
				}
				else
					intersections.push_back(pair<double,double*>(t4,(double*)NULL));
			}
		}	
	}
	
	if (intersections.size()>0)
		sort(intersections.begin(), intersections.end(),compKeepIndex);	
}

void ConnollySurface::clear()
{
	// deallocate all patches
	for (unsigned int i=0;i<sesComplex.size();i++)
			delete sesComplex[i];

	// deallocate acceleration grid 2D
	if (gridConnollyCellMap2D!=NULL)
		deleteVector<int>(gridConnollyCellMap2D);
	
	if (ind_2d!=NULL)
		deleteMatrix2D<unsigned int>(last_cols_ind,ind_2d);

	// deallocate acceleration grid 3D
	if (gridConnollyCellMap!=NULL)
		deleteVector<int>(gridConnollyCellMap);

	// deallocate vector
	if (ind!=NULL)
		deleteMatrix3D<unsigned short>(nx,ny,ind);

	if (x!=NULL)
		deleteVector<double>(x);

	if (y!=NULL)
		deleteVector<double>(y);

	if (z!=NULL)
		deleteVector<double>(z);

	// deallocate exposed atoms pointers
	if (atomPatches!=NULL)
		deleteVector<PointCell*>(atomPatches);

	sesComplex.clear();	
}
ConnollySurface::~ConnollySurface()
{
	clear();
}

bool ConnollySurface::rayConnollyCellIntersection(double* orig,double* dir,ConnollyCell* cc,double* t1,double* t2,double* t3,double* t4,int& numInt,int thdID)
{
	double* sphere_center;
	double* torus_center;
	EdgeCell* ec;
	double sphere_radius;
	int patch_type = cc->patch_type;

	bool doTorus=false;

	if (patch_type==REGULAR_FACE_CELL || patch_type==SINGULAR_FACE_CELL)
	{
		FacetCell* fc = (FacetCell*)cc;		
		sphere_center = fc->center;
		sphere_radius = probe_radius;
	}
	else if (patch_type==SINGULAR_EDGE_CELL || patch_type==REGULAR_EDGE_CELL)
	{
		ec = (EdgeCell*)cc;		
		sphere_center = ec->clipping_center;
		sphere_radius = ec->clipping_radius;
		torus_center = ec->center;
		doTorus=true;
	}
	else if (patch_type==POINT_CELL)
	{
		PointCell* pc = (PointCell*)cc;
		sphere_center = delphi->atoms[pc->id]->pos;
		sphere_radius = delphi->atoms[pc->id]->radius;		
	}
	else
		cout << endl << ERR << "Cannot get patch type in intersection test";

	bool det = raySphere(orig,dir,sphere_center,sphere_radius,t1,t2);

	if (!det)
		return false;

	numInt = 2;

	// if it is a sphere we have finished
	if (!doTorus)
	{		
		return true;
	}

	int index = 0;
	double dir_norm;

	double torus_roots[4],pp[3],ddir[3],tdir[3],temp[3];

	double orig1[3],orig2[3];
	// start the ray from the clipping sphere intersection
	ADD_MUL(orig1,orig,dir,(*t1))
	ADD_MUL(orig2,orig,dir,(*t2))
	SUB(tdir,orig2,orig1)

	// build torus roots equation
	SUB(temp,orig1,torus_center)
		
	// roto-translate origin point
	for (int i=0;i<3;i++)
	{
		pp[i]=0;
		for (int j=0;j<3;j++)			
			pp[i]+=ec->invrot[i][j]*temp[j];
	}

	double tt;
	// rotate ray direction
	for (int i=0;i<3;i++)
	{
		ddir[i]=0;
		for (int j=0;j<3;j++)
			ddir[i]+=ec->invrot[i][j]*tdir[j];
	}

	NORMALIZE_S(ddir,tt)		
	dir_norm = tt;

	double pd,p2,coeff,r2,R2;
	r2 = probe_radius*probe_radius;
	R2 = ec->major_radius*ec->major_radius;	
	pd = DOT(pp,ddir);
	p2 = DOT(pp,pp);
	coeff = (p2-r2-R2);

	double A = 1;
	double B = 4*pd;
	double C = 2*(p2-(R2+r2)+2*pd*pd+2*R2*ddir[2]*ddir[2]);
	double D = 4*(pd*coeff+2*R2*pp[2]*ddir[2]);
	double E = coeff*coeff-4*R2*(r2-pp[2]*pp[2]);

	int numrootss = 4;
	double polyy[5];

	polyy[0]=E;
	polyy[1]=D;
	polyy[2]=C;
	polyy[3]=B;
	polyy[4]=A;

	double 	roots[4];
	int	 order=4;
	
	getRealRootsSturm(polyy,order,roots,numrootss);
	
	if (numrootss==0)
		return false;

	for (int i = 0; i != numrootss; i++)
		torus_roots[i]=roots[i];

	vector<double> ll;

	for (int i=0;i<numrootss;i++)
	{
		if (torus_roots[i]<0)
			continue;
		
		torus_roots[i] = torus_roots[i]/dir_norm;		
		double dd,point[3];
		ADD_MUL(point,orig1,tdir,torus_roots[i])	

		DIST2(dd,sphere_center,point)
		// acceptable
		if (dd<(sphere_radius*sphere_radius))
		{		
			SUB(point,point,orig)
			for (int k=0;k<3;k++)
			{
				if (dir[k]!=0)
				{
					ll.push_back(point[k]/dir[k]);
					break;
				}
			}
		}
	}

	if (ll.size()==0)
		return false;

	numInt = (int)ll.size();
	for (unsigned int i=0;i<ll.size();i++)
	{
		if (i==0)
			(*t1)=ll[i];
		else if (i==1)
			(*t2)=ll[i];
		else if (i==2)
			(*t3)=ll[i];
		else
			(*t4)=ll[i];
	}
	return true;
}

bool ConnollySurface::isFeasible(ConnollyCell* cc,double* point)
{	
	if (cc->patch_type==REGULAR_FACE_CELL || cc->patch_type==SINGULAR_FACE_CELL)
	{
		FacetCell* fc = (FacetCell*)cc;
		int start;

		if (cc->patch_type==SINGULAR_FACE_CELL)
			start = 0;
		else
			start = 1;

		for (int i=start;i<4;i++)
		{
			double tst = DOT(fc->planes[i],point)+fc->planes[i][3];
			if (tst>0)
				return false;
		}

		for (unsigned int i=0;i<fc->self_intersection_planes.size();i++)
		{
			double tst = DOT(fc->self_intersection_planes[i],point)+fc->self_intersection_planes[i][3];
			if (tst>0)
				return false;
		}

		return true;
	}
	else if (cc->patch_type==SINGULAR_EDGE_CELL)
	{
		EdgeCell* ec = (EdgeCell*)cc;
		
		if (ec->isSelfIntersecting)
		{			
			double r = ec->self_intersection_radius;
			double dd;
			DIST2(dd,ec->center,point)			
			if (dd<(r*r))
				return false;
		}
		return true;
	}
	else if (cc->patch_type==REGULAR_EDGE_CELL)
	{
		EdgeCell* ec = (EdgeCell*)cc;

		if (ec->isSelfIntersecting)
		{			
			// check for self intersection
			double r = ec->self_intersection_radius;
			double dd;
			DIST2(dd,ec->center,point)
			if (dd<(r*r))
				return false;
		}
		
		// it is inside clipping sphere, check the clipping planes
		if (ec->additional_planes.size()!=0)
		{
			int plane_index = 0;
			for (unsigned int i = 0;i<ec->flags.size();i++,plane_index+=2)
			{
				bool acute = ec->flags[i];

				if (acute)
				{
					double tst = DOT(ec->additional_planes[plane_index],point)+ec->additional_planes[plane_index][3];
					if (tst>0)
						continue;

					tst = DOT(ec->additional_planes[plane_index+1],point)+ec->additional_planes[plane_index+1][3];
					if (tst>0)
						continue;

					return true;
				}
				else
				{
					double tst = DOT(ec->additional_planes[plane_index],point)+ec->additional_planes[plane_index][3];
					if (tst<0)
						return true;

					tst = DOT(ec->additional_planes[plane_index+1],point)+ec->additional_planes[plane_index+1][3];
					if (tst<0)
						return true;

					continue;			
				}
			}
			// no planes pair gives a positive result
			return false;
		}
		// fast path
		else
		{
			if (ec->acute)
			{
				double tst = DOT(ec->cutting_planes[0],point)+ec->cutting_planes[0][3];
				if (tst>0)
					return false;

				tst = DOT(ec->cutting_planes[1],point)+ec->cutting_planes[1][3];
				if (tst>0)
					return false;

				return true;
			}
			else
			{
				double tst = DOT(ec->cutting_planes[0],point)+ec->cutting_planes[0][3];
				if (tst<0)
					return true;

				tst = DOT(ec->cutting_planes[1],point)+ec->cutting_planes[1][3];
				if (tst<0)
					return true;

				return false;			
			}
		}

		return false;
	}
	else if (cc->patch_type==POINT_CELL)
	{
		PointCell* pc = (PointCell*)cc;
		
		// filter on tori clipping spheres (shifted voronoi planes of exposed atoms)
		for (unsigned int i=0;i<pc->neighbours.size();i++)
		{
			double* center = pc->neighbours[i]->clipping_center;
			double radius = pc->neighbours[i]->clipping_radius;
			double dd;
			DIST2(dd,center,point)
			dd-=radius*radius;			
			if (dd<0)
				return false;
		}

		// filter on buried clipping tori (shifted voronoi planes of buried atoms)
		for (unsigned int i=0;i<pc->buried_neighbours.size();i++)
		{
			double* center = pc->buried_neighbours[i]->clipping_center;
			double radius = pc->buried_neighbours[i]->clipping_radius;
			double dd;
			DIST2(dd,center,point)
			dd-=radius*radius;
			if (dd<0)
				return false;
		}
		
		return true;
	}
		
	cout << endl << ERR << "Cannot get an answer in feasibility test";
	return false;
}


void ConnollySurface::projectToCircle(double* point,double radius,double* center,double* plane,double* proj,double& dist)
{
	double pdist;
	point2plane(point,plane,&pdist,proj);
	SUB(proj,proj,center)
	double n = sqrt((DOT(proj,proj)));
	radius = radius/n;
	ADD_MUL(proj,center,proj,radius)
	DIST(dist,point,proj)
}

void ConnollySurface::projectToTorus(double* y,EdgeCell* ec,double* proj,double* norm,double& dist)
{	
	double temp[3],pp[3],torus_plane[4]={0,0,+1,0},proj_major[3],minor_plane[4],torus_center[3]={0,0,0};
			
	// roto-translate the point into the torus reference 	
	SUB(temp,y,ec->center)
	for (int i=0;i<3;i++)
	{
		pp[i]=0;
		for (int j=0;j<3;j++)			
			pp[i]+=ec->invrot[i][j]*temp[j];
	}

	// first project to the big circle
	projectToCircle(pp,ec->major_radius,torus_center,torus_plane,proj_major,dist);

	// need 2 more points to get the plane of the previously identified circle. 
	// Two good points are along the axis together with center of the circle.
	double p2[3]={0,0,+1},p3[3]={0,0,-1};
	plane3points(proj_major,p2,p3,minor_plane);
	
	// second project to the small circle
	projectToCircle(pp,probe_radius,proj_major,minor_plane,proj,dist);

	// deroto-translate				
	for (int i=0;i<3;i++)
	{
		pp[i]=0;
		for (int j=0;j<3;j++)			
			pp[i]+=ec->Rot[i][j]*proj[j];
	}

	ADD(proj,pp,ec->center)
	DIST(dist,proj,y)
}

void ConnollySurface::getNormal(double* y,ConnollyCell* cc,double* normal)
{
	int patch_type = cc->patch_type;
	double *sphere_center,sphere_radius;

	if (patch_type==REGULAR_FACE_CELL || patch_type==SINGULAR_FACE_CELL)
	{
		FacetCell* fc = (FacetCell*)cc;
		sphere_center = fc->center;
		sphere_radius = probe_radius;
		SUB(normal,sphere_center,y)
		normal[0]/=sphere_radius;
		normal[1]/=sphere_radius;
		normal[2]/=sphere_radius;
	}	
	else if (patch_type==SINGULAR_EDGE_CELL || patch_type==REGULAR_EDGE_CELL)
	{
		EdgeCell* ec = (EdgeCell*)cc;		
		getNormalToTorus(y,ec,normal);		
	}	
	else if (patch_type==POINT_CELL)
	{
		PointCell* pc = (PointCell*)cc;
		sphere_center = delphi->atoms[pc->id]->pos;
		sphere_radius = delphi->atoms[pc->id]->radius;
		SUB(normal,y,sphere_center)
		normal[0]/=sphere_radius;
		normal[1]/=sphere_radius;
		normal[2]/=sphere_radius;
	}
}
void ConnollySurface::getNormalToTorus(double* y,EdgeCell* ec,double* normal)
{	
	double temp[3],pp[3],torus_plane[4]={0,0,+1,0},torus_center[3]={0,0,0},dist,C[3],aa;
			
	// roto-translate the point into the torus reference 	
	SUB(temp,y,ec->center)
	for (int i=0;i<3;i++)
	{
		pp[i]=0;
		for (int j=0;j<3;j++)			
			pp[i]+=ec->invrot[i][j]*temp[j];
	}

	projectToCircle(pp,ec->major_radius,torus_center,torus_plane,C,dist);
	SUB(C,C,pp)
	NORMALIZE(C,aa)

	// de-rotate the normal vector
	for (int i=0;i<3;i++)
	{
		normal[i]=0;
		for (int j=0;j<3;j++)			
			normal[i]+=ec->Rot[i][j]*C[j];
	}
}



void ConnollySurface::saveConcaveSpherePatch(ofstream& of,FacetCell* fc,int i)
{	
	char buff2[BUFLEN],buff[BUFLEN],temp[BUFLEN];
	sprintf(buff,"Mesh_CS%d",i);
	of << "\n\n #declare " << buff << "=" ;

	sprintf(buff2,"\n\n intersection { \n ");
	of << endl << buff2;

	for (int l=0;l<4;l++)
	{
		sprintf(temp,"plane{<%f,%f,%f>,%f}",fc->planes[l][0],fc->planes[l][1],fc->planes[l][2],-(fc->planes[l][3])/(sqrt((DOT(fc->planes[l],fc->planes[l])))));				
		of << endl << temp;
	}

	//cout << endl << "Num " << fc->numSelfIntersections;
	if (fc->isSelfIntersecting)
	{
		of << endl << "// self intersection planes" << endl;
		for (unsigned int l=0;l<fc->self_intersection_planes.size();l++)
		{
			//cout << endl << fc->self_intersection_planes[l][0] << " " << fc->self_intersection_planes[l][1] <<" "<< fc->self_intersection_planes[l][2]
			sprintf(temp,"plane{<%f,%f,%f>,%f}",fc->self_intersection_planes[l][0],fc->self_intersection_planes[l][1],fc->self_intersection_planes[l][2],-(fc->self_intersection_planes[l][3])/(sqrt((DOT(fc->self_intersection_planes[l],fc->self_intersection_planes[l])))));				
			of << endl << temp;				
		}
	}
	//cout << endl << "-------";

	of << "\n}";
	// save the sphere with its radius
	sprintf(buff2,"\n\n sphere { \n <%lf,%lf,%lf>,%lf",fc->center[0],fc->center[1],fc->center[2],probe_radius);	  
	of << buff2;
//	if (fc->patch_type==SINGULAR_FACE_CELL)
//		of << "\n pigment{ color Orange}";	  
//	else
		of << "\n pigment{ color Red}";	  
	of << "\n bounded_by { " << buff << " } \n clipped_by{ bounded_by}}";			
}


void ConnollySurface::saveSphere(ostream& of,double* center,double radius)
{
	char buff2[BUFLEN];
	sprintf(buff2,"\n\n sphere { \n <%lf,%lf,%lf>,%lf  pigment{color Green}}",center[0],center[1],center[2],radius);	  				
	of << buff2;
}

void ConnollySurface::saveAtomPatch(ofstream& of,PointCell* pc)
{
	char buff2[BUFLEN],buff[BUFLEN];

	of << "\n// -------------- ATOM " << pc->id << "-------------------";
	sprintf(buff,"AtomClip%d",pc->id);
	// add 
	of << "\n\n #declare " << buff << " =" ;	
	of << "\n union {";
	double* sphere;
	double radius;
	for (unsigned int k=0;k<pc->neighbours.size();k++)
	{
		sphere = pc->neighbours[k]->clipping_center;
		radius = pc->neighbours[k]->clipping_radius;
		sprintf(buff2,"\n\n sphere { \n <%lf,%lf,%lf>,%lf}",sphere[0],sphere[1],sphere[2],radius);	  				
		of << buff2;
	}
	for (unsigned int k=0;k<pc->buried_neighbours.size();k++)
	{
		sphere = pc->buried_neighbours[k]->clipping_center;
		radius = pc->buried_neighbours[k]->clipping_radius;
		sprintf(buff2,"\n\n sphere { \n <%lf,%lf,%lf>,%lf}",sphere[0],sphere[1],sphere[2],radius);	  				
		of << buff2;
	}

	of << "\n}";

	of << "\ndifference{";
	sprintf(buff2,"\n\n sphere { \n <%lf,%lf,%lf>,%lf  pigment{color Green}}",delphi->atoms[pc->id]->pos[0],delphi->atoms[pc->id]->pos[1],delphi->atoms[pc->id]->pos[2],delphi->atoms[pc->id]->radius);	  				
	of << buff2;
	of << "\nobject{" << buff << "}}";	
}


void ConnollySurface::saveEdgePatch(ofstream& of,EdgeCell* ec,int size,double bigR,double* u,double* v,double* w,bool isComplex)
{
	char buff2[BUFLEN],buff[BUFLEN],buff3[BUFLEN],clipplanes[BUFLEN],clip[BUFLEN];											

	of << "\n// ------------- torus -------------//";
	of << "\n// pair " << ec->id[0] << "," << ec->id[1];
	of << "\n\n";
	sprintf(clipplanes,"Clipping_Planes%d",size);

	// singular edge cell has not clipping planes
	if (ec->patch_type!=SINGULAR_EDGE_CELL)				
	{
		// add clipping planes
		of << "\n\n #declare " << clipplanes << " =" ;						

		if (ec->additional_planes.size()==0)
		{
			if (ec->acute)
				sprintf(buff3,"\n\n // clipping planes \n\n intersection { \n ");
			else
				sprintf(buff3,"\n\n // clipping planes \n\n union { \n ");

			of << buff3;
			for (int l=0;l<2;l++)
			{
				sprintf(buff3,"\nplane{<%f,%f,%f>,%f}",ec->cutting_planes[l][0],ec->cutting_planes[l][1],ec->cutting_planes[l][2],-(ec->cutting_planes[l][3])/(sqrt((DOT(ec->cutting_planes[l],ec->cutting_planes[l])))));				
				of << buff3;
			}
			of << "\n}";
		}
		else
		{
			of << "\nunion{\n";
			int jj=0;
			for (unsigned int i=0;i<ec->additional_planes.size();i+=2)
			{
				if (ec->flags[jj])
					sprintf(buff3,"\n\n // clipping planes \n\n intersection { \n ");
				else
					sprintf(buff3,"\n\n // clipping planes \n\n union { \n ");

				jj++;

				of << buff3;
				sprintf(buff3,"\nplane{<%f,%f,%f>,%f}",ec->additional_planes[i][0],ec->additional_planes[i][1],ec->additional_planes[i][2],-(ec->additional_planes[i][3])/(sqrt((DOT(ec->additional_planes[i],ec->additional_planes[i])))));				
				of << buff3;
				sprintf(buff3,"\nplane{<%f,%f,%f>,%f}",ec->additional_planes[i+1][0],ec->additional_planes[i+1][1],ec->additional_planes[i+1][2],-(ec->additional_planes[i+1][3])/(sqrt((DOT(ec->additional_planes[i+1],ec->additional_planes[i+1])))));				
				of << buff3;

				of << "\n}";

			}
			of << "\n}";
		}					
	}				

	sprintf(clip,"Clip%d",size);
	of << "\n\n #declare " << clip << " =" ;	

	// manage self intersection if necessary 
	if ((bigR-probe_radius)<0)
	{					
		sprintf(buff2,"\n\n difference { \n "
			"sphere {  <%lf,%lf,%lf>,%lf }  "
			"\n sphere {  <%lf,%lf,%lf>,%lf } }",
			//(clipping_center.x()),(clipping_center.y()),(clipping_center.z()),r_clipping,
			//torus_center.x(),torus_center.y(),torus_center.z(),sqrt(-bigR*bigR+probe_radius*probe_radius));
			(ec->clipping_center[0]),(ec->clipping_center[1]),(ec->clipping_center[2]),ec->clipping_radius,
			ec->center[0],ec->center[1],ec->center[2],sqrt(-bigR*bigR+probe_radius*probe_radius));
		of << buff2;
	}
	else
	{
		sprintf(buff2,"\nsphere {  <%lf,%lf,%lf>,%lf }  ",
			//(clipping_center.x()),(clipping_center.y()),(clipping_center.z()),r_clipping);
			(ec->clipping_center[0]),(ec->clipping_center[1]),(ec->clipping_center[2]),ec->clipping_radius);
		of << buff2;				
	}


	// global clipping object 
	sprintf(buff,"Clipping_Object%d",(int)sesComplex.size());

	// add all
	of << "\n\n #declare " << buff << " =" ;	
	// torus,clipping spheres, and also clipping planes
	if (ec->patch_type!=SINGULAR_EDGE_CELL)
		sprintf(buff2,"\n\n intersection { \n"
		"object{%s}\nobject{%s}}",clipplanes,clip);
	// only torus with clipping spheres
	else
		sprintf(buff2,"\n\nobject{%s}",clip);

	of << buff2;


	sprintf(buff2,"\n\n torus { \n %lf,%lf",bigR,probe_radius);	  
	of << buff2;
//	if (!isComplex)
		of << "\n pigment{ color White}";
//	else
//		of << "\n pigment{ color Yellow}";

	//sprintf(buff2,"\n matrix <%lf,%lf,%lf,\n%lf,%lf,%lf,\n%lf,%lf,%lf,\n0,0,0>",	U.x(),U.y(),U.z(),
	//																				w.x(),w.y(),w.z(),
	//																				V.x(),V.y(),V.z());

	sprintf(buff2,"\n matrix <%lf,%lf,%lf,\n%lf,%lf,%lf,\n%lf,%lf,%lf,\n0,0,0>",	u[0],u[1],u[2],
		w[0],w[1],w[2],
		v[0],v[1],v[2]);
	of << buff2;
	//of << "\n translate"<< "<" << torus_center.x() << "," << torus_center.y() << "," << torus_center.z() << ">";
	of << "\n translate"<< "<" << ec->center[0] << "," << ec->center[1] << "," << ec->center[2] << ">";
	of << "\nsturm";
	of << "\n bounded_by { " << buff << " } \n clipped_by{ bounded_by}}"; 

	// to draw the clipping sphere
	//sprintf(buff2,"\nsphere {  <%lf,%lf,%lf>,%lf   \n pigment{ color Yellow filter 0.7}}",
	//(ec->clipping_center[0]),(ec->clipping_center[1]),(ec->clipping_center[2]),ec->clipping_radius);
	//of << buff2;					
}

/** check the orientation. Assume the planes points toward the visible region of the torus*/
bool ConnollySurface::orientation(double* pb_center1,double* pb_center2,double* w1,double* w2)
{
	// intersect the lines starting at the probes position and following the normal orientations
	//double a=DOT(w1,w1);
	double b=DOT(w1,w2);
	double c=DOT(w2,w2);
	double w0[3];
	SUB(w0,pb_center1,pb_center2)
		double d=DOT(w1,w0);
	double e=DOT(w2,w0);

	//double sc = (b*e-c*d)/(a*c-b*b);
	// only sign
	double sc = (b*e-c*d);

	//cout << endl << "sc param " << sc;
	if (sc<0)
		return true;
	else 
		return false;
}

void ConnollySurface::sortProbes(EdgeCell* ec,FacetCell** fcv,int np,int* sorted)
{
	vector<indexed_double> angles;
	double ref[3];
	double tt;
	SUB(ref,fcv[0]->center,ec->center)
		NORMALIZE_S(ref,tt);
	double ws[3],temp[3],vs[3];
	SUB(temp,fcv[1]->center,ec->center)
		NORMALIZE_S(temp,tt);
	CROSS(vs,ref,temp)
		CROSS(ws,vs,ref)

		for (int i=1;i<np;i++)
		{
			double temp[3];
			SUB(temp,fcv[i]->center,ec->center)
				NORMALIZE_S(temp,tt)
				double a1 = DOT(ws,temp);
			double a2 = DOT(ref,temp);
			double angle = atan2(a1,a2);
			// [0,2pi]	
			if (angle<0)
				angle = TWO_PI+angle;
			//cout << endl << "Angle " << angle;
			angles.push_back(indexed_double(angle,i));
		}		
		// remember indexing and sort
		sort(angles.begin(),angles.end(),index_double_comparator);
		// save
		sorted[0]=0;
		//cout << endl << "Ordered";
		for (int i=1;i<np;i++)
		{
			sorted[i] = angles[i-1].second;
			//cout << endl << "\t" << sorted[i];			
		}
		return;
}



void ConnollySurface::getCoi(double* torus_center,double rcoi,double** sampledPoints,int numPoints,double* u,double* v)
{
	//FILE* fp;
	//fp = fopen("tempp.txt","a");

	double step = 6.28/numPoints;
	double tcoi =0;
	for(int i=0;i<numPoints;i++)
	{			
		sampledPoints[i][0]=torus_center[0]+rcoi*(cos(tcoi)*u[0]+sin(tcoi)*v[0]);
		sampledPoints[i][1]=torus_center[1]+rcoi*(cos(tcoi)*u[1]+sin(tcoi)*v[1]);
		sampledPoints[i][2]=torus_center[2]+rcoi*(cos(tcoi)*u[2]+sin(tcoi)*v[2]);			
		tcoi += step;

		//fprintf(fp,"\nsphere {"); 
		//fprintf(fp,"\n<%f,%f,%f>,0.1\n",sampledPoints[i][0],sampledPoints[i][1],sampledPoints[i][2]);
		//fprintf(fp,"\npigment{ color White}}"); 			
	}			
	//fclose(fp);
}



///////////////////////////////////////

// EXTRA NEW
// void ConnollySurface :: saveProbeCenters(ofstream& oftxt,ofstream& ofpqr,FacetCell* fc,const int i) const {
void ConnollySurface :: saveProbeCenters(ofstream& oftxt,FacetCell* fc,const int i) const {
	char temp1[BUFLEN],temp2[BUFLEN];

	sprintf(temp1, "%f %f %f %f \t %d %d %d\t %d", fc->center[0], fc->center[1], fc->center[2], probe_radius, fc->id[0], fc->id[1], fc->id[2], fc->isSelfIntersecting);
	oftxt << temp1 << std:: endl;
	//   ** PQR FORMAT ** 
	//i atom index needed bt pqr format, 0 is the charge of the atom
	// sprintf(temp2, "ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f%8.4f%8.4f",i%100000,"C","NSR",'A',1,fc->center[0],fc->center[1],fc->center[2],0.,probe_radius);
	// ofpqr << temp2 << std:: endl;
}

//

bool ConnollySurface :: write_regTri(Rt& reg_triang) const{
	const bool FAIL =0;
	const bool SUCCESS =1;
	char temp[BUFLEN];

	std :: cout << "\n";
	std :: cout << endl << INFO<< "Fetching Regular Triangulation" ;
	FILE* fp = fopen("regTriang.OFF","w");
	sprintf(temp,"# Regular triangulation with atomic radii augmented with probe radius= %f",probe_radius);
	std::string title = temp;
	std :: string buff(title.size() + 15, ' ');
	fprintf(fp,"%s",buff.c_str());

	
	std :: map <Weighted_point,int> vert_map;
	// 3 vertices per triangle
	Weighted_point v1,v2,v3;

	Finite_Vertex_Iterator fvit = reg_triang.finite_vertices_begin();
	int i = 0;
	std :: cout << endl << INFO << "Building vertex list..";
	std :: cout << endl << INFO << "Number of vertices regular triangulation: "<< reg_triang.number_of_vertices();
	for (;reg_triang.finite_vertices_end()!=fvit;fvit++){
		vert_map.insert(pair<Weighted_point,int> (fvit->point(),i));
		i++;
		fprintf(fp,"\n %f %f %f",fvit->point().x(),fvit->point().y(),fvit->point().z());
	}
	Finite_Facets_Iterator ffit = reg_triang.finite_facets_begin();
	int k=0;
	std :: cout << endl << INFO << "Iterating through facets..";
	for(;ffit!=reg_triang.finite_facets_end();ffit++){
		//const Facet& facetp = reg_triang.mirror_facet((*ffit));
		const Facet& facetp = (*ffit);
		v1= facetp.first->vertex((facetp.second+1)&3)->point();
		v2= facetp.first->vertex((facetp.second+2)&3)->point();
		v3= facetp.first->vertex((facetp.second+3)&3)->point();
		fprintf(fp,"\n3 %d %d %d",vert_map[v1],vert_map[v2], vert_map[v3] );
		if ( (vert_map.find(v1) == vert_map.end())||(vert_map.find(v2) == vert_map.end())||(vert_map.find(v3) == vert_map.end()) )
		{
			std :: cout << endl<< WARN << "ERROR, vertex not mapped\n";
			return FAIL;
		}
		k++;
	}

	std :: cout << endl << INFO << "Done: Number of facets: "<< k << std ::endl;
	//std :: cout << reg_triang.number_of_facets() << std ::endl; here counts also those connected to vertex to inf I think.. whilst k iterates through finite facets
	rewind(fp);
	fprintf(fp,"%s\nOFF %d %d 0\n",title.c_str(),i,k);
	fclose(fp);


	return SUCCESS;
} 

bool ConnollySurface :: write_alphaFilt(Fixed_alpha_shape_3& alpha, const int flag) const{
	/*info: 
		CGAL classifiers:
			REGULAR = not singular AND on the boundary: that is alpha-exposed
			EXTERIOR = I think, not alpha-exposed and not interior.
	 		SINGULAR = Exterior and ONLY connected to other exterior elements
			INTERIOR = Within the boundary (REGULAR), not alpha exposed
			OTHER OPTIONS THAN REGULAR AND INTERIOR NOT SUPPORTED YET!
		Rationale:
			The function creates an OFF file. All possible vertex coordinates (independently from the CGAL classification)
			are printed and a map associating a unique label to each vertex is created.
			Then the facets corresponding to the chosen classification are printed on the file
			by listing the associated vertices labels (3 per facet, one facet per line)
	*/
	const int REG = 0;
	const int INT = 1;
	const int EXT = 2;
	const int SING =3;
	const bool FAIL =0;
	const bool SUCCESS =1;
	char temp[BUFLEN];

	std :: cout << endl << INFO<< "Fetching alpha-simplex" ;
	if(flag==REG){
		 std :: cout << endl << INFO<< "tag = REGULAR ";
		 sprintf(temp,"# Boundary of the 0-alpha shape (alpha-exposed) with atomic radii augmented with probe radius= %f",probe_radius);
		 }
	else if(flag==INT) {
		std :: cout << endl << INFO<< "tag = INTERNAL";
		sprintf(temp,"# Interior of the 0-alpha shape (alpha-exposed) with atomic radii augmented with probe radius= %f",probe_radius);
		}
	// else if(flag==EXT) std :: cout << endl << INFO<< "tag = EXTERNAL";
	FILE* fp = fopen("alphaTriang.OFF","w");
	//FILE* fpe = fopen("alphaTriang_discarded.OFF","w");
	


	std::string title = temp;
	std :: string buff(title.size() + 15, ' ');

	fprintf(fp,"%s",buff.c_str());

	// here store unique map for vertices
	std :: map <Weighted_point,int> vert_map;
	// 3 vertices per triangle
	Weighted_point v1,v2,v3;

	Finite_Vertex_Iterator fvit = alpha.finite_vertices_begin();
	int i = 0;
	
	std :: cout << endl << INFO << "Building vertex list..";
	//creating "verbose vertex list" some won't be used in the OFF file
	// with some changes, should not be too hard to just print vertices that will be used..
	for (;alpha.finite_vertices_end()!=fvit;fvit++){
		int vtype = alpha.classify(fvit);
				vert_map.insert(pair<Weighted_point,int> (fvit->point(),i));
				i++;
				//save top half of OFF file, WHAT ABOUAT WEIGHTS?
				fprintf(fp,"\n %f %f %f",fvit->point().x(),fvit->point().y(),fvit->point().z());
	}
	std :: cout << endl << INFO << "Total Number of vertices: "<< i;

	// assign vertex to correct facet(s) 
	std :: cout << endl << INFO << "Iterating through facets..";
	Finite_Facets_Iterator ffit = alpha.finite_facets_begin();
	int k=0;
	for(;ffit!=alpha.finite_facets_end();ffit++){
		int ttype = alpha.classify(*ffit);
		if(flag == REG){
			if(ttype == Fixed_alpha_shape_3:: REGULAR){
				const Alpha_Facet& facetp = (*ffit);
				v1= facetp.first->vertex((facetp.second+1)&3)->point();
				v2= facetp.first->vertex((facetp.second+2)&3)->point();
				v3= facetp.first->vertex((facetp.second+3)&3)->point();

				if ( (vert_map.find(v1) == vert_map.end())||(vert_map.find(v2) == vert_map.end())||(vert_map.find(v3) == vert_map.end()) )
				{
					std :: cout << endl<< WARN << "ERROR, vertex not mapped while retrieving alpha shape\n";
					return FAIL;
				}
				fprintf(fp,"\n3 %d %d %d",vert_map.at(v1),vert_map.at(v2), vert_map.at(v3) );
				k++;
			}
		}
		else if (flag==INT){
			if(ttype == Fixed_alpha_shape_3:: INTERIOR)
			{	
				const Alpha_Facet& facetp = alpha.mirror_facet((*ffit)); //for consistency to what done when building patches
				//const Alpha_Facet& facetp = (*ffit);
				v1= facetp.first->vertex((facetp.second+1)&3)->point();
				v2= facetp.first->vertex((facetp.second+2)&3)->point();
				v3= facetp.first->vertex((facetp.second+3)&3)->point();

				if ( (vert_map.find(v1) == vert_map.end())||(vert_map.find(v2) == vert_map.end())||(vert_map.find(v3) == vert_map.end()) )
				{
					std :: cout << endl << WARN << "ERROR, vertex not mapped while retrieving alpha shape\n";
					return FAIL;
				}
				fprintf(fp,"\n3 %d %d %d",vert_map.at(v1),vert_map.at(v2), vert_map.at(v3) );	
				k++;	
			}
		}
		else{
			std :: cout << endl << WARN<<"Invalid option for printing alpha shape.\n"<<std::endl;
			return FAIL;
		}
	}
	//move to the top and add extra line with number of vertices and facets, as required by OFF format
	rewind(fp);
	//fprintf(fp,"OFF %d %d 0\n",i,k);
	fprintf(fp,"%s\nOFF %d %d 0\n",title.c_str(),i,k);
	//fp << "OFF" << i<< k << std::endl;
	fclose(fp);
	//print off file
	std :: cout << endl << INFO << "Done: Number of facets with indicated attributes: "<< k << "\n "<< std ::endl;
	return SUCCESS;


}