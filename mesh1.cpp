#include "mesh.h"
#include "dg.h"
#include <cmath>

std::vector<std::vector<double>> reader::readv(std::string s)
{
	std::vector<std::vector<double>> v;
	std::ifstream input_file(s);
	if (input_file.is_open())
	{
		std::string line;
		while (std::getline(input_file, line))
		{
			std::stringstream line_stream;
			line_stream << line;
			std::vector<double> tempvec;
			double val;
			std::string word;
			while (line_stream >> word)
			{
				if (std::stringstream(word) >> val)
				{
					tempvec.push_back(val);
				}
			}
			v.push_back(tempvec);
		}
	}
	else
	{
		std::cout << "file " << s << " not found" << std::endl;
	}
	return v;
}

// construtor for mesh class generates intpoel,coords and bface data sructures from mesh file
grid::mesh::mesh(std::string s1, std::string s2)
{
  
	func_count = 1; //start function counter to maintain order
	std::vector<std::vector<double>> v = reader::readv(s1);
	std::vector<std::vector<double>> c = reader::readv(s2);
	// get all the int members of the class
	ndimn = v[5][0];
	ntype = v[5][1];
	nelem = v[7][0];
	npoin = v[7][1];
	nbface = v[7][2];
	neqns = ndimn + 2;
	ndegr = (c[0][0] + 1) * (c[0][0] + 2) * 0.5;
  ngauss_boun = c[1][0];
  ngauss_domn = c[2][0];
  domweight = c[3][0];
  bounweight = c[4][0];
	U_infty.init(neqns, 1); //intial condition
	U_infty(0,0) = c[5][0]; //reference density
	U_infty(1,0) = cos(c[6][0]*const_properties::pi/180);//X velocity
	U_infty(2,0) = sin(c[6][0]*const_properties::pi/180); //Y velocity
	U_infty(3,0) = 0.5 + 1/(const_properties::gamma*(const_properties::gamma-1)*c[7][0]*c[7][0]);
	// for loop populate intpoel matrix
	inpoel.init(nelem, ntype); // init and give size to inpoel
	for (int i = 0; i < nelem; i++)
	{
		for (int j = 0; j < ntype; j++)
		{
			inpoel(i, j) = v[i + 9][j + 1];
		}
	}
	coords.init(npoin, ndimn); // init and give size to coords matrix
	for (int i = 0; i < npoin; i++)
	{
		for (int j = 0; j < ndimn; j++)
		{
			coords(i, j) = v[i + 10 + nelem][j + 1];
		}
	}
	// bface format
	//  ith row (ip1 ip2 host el ghost el face_flag)
	bface.init(nbface, 5); 
	for (int i = 0; i < nbface; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			bface(i, j) = v[i + nelem + 2 * npoin + 12][j + 1];
		}
		bface(i, 4) = v[i + nelem + 2 * npoin + 12][3];
	}
}
// end constructor

// sub to generate elemenst surrounding points linked list
void grid::pre_proc::set_esup(grid::mesh &mesh1)
{
	assert(mesh1.func_count == 1); // function must run after mesh constructor
	mesh1.func_count = 2;		   // update func_counter
	mesh1.esup2.init(mesh1.npoin + 1, 1); // initialize 
	// this pass counts number of elements surrounding every point and stores that number in +1 index corresponding to that particular point
	// for eg if point 9 has 6 elements around it, then
	//  the number stored in index 10(for c++ it will be index 9) will be the number 6 after the first element pass
	for (int i = 0; i < mesh1.nelem; i++) // loop over all elements
	{
		for (int j = 0; j < mesh1.ntype; j++) // loop over all the points of a node
		{
			mesh1.esup2(mesh1.inpoel(i, j), 0)++;
		}
	}
	// reshuffle pass for generating storage requirement for the esup1 array
	for (int i = 1; i < mesh1.npoin + 1; i++)
	{
		mesh1.esup2(i, 0) = mesh1.esup2(i, 0) + mesh1.esup2(i - 1, 0); // pushing back one index for each point
	}

	int storage_req = mesh1.esup2(mesh1.npoin, 0);
	// mesh1.esup1.resize(storage_req, 1); // init esup1 based on storage req
	mesh1.esup1.init(storage_req, 1);
	// element pass 2, storage
	for (int i = 0; i < mesh1.nelem; i++) // loop over all elements
	{
		for (int j = 0; j < mesh1.ntype; j++) // loop over all nodes of the element
		{
			int ipoin = mesh1.inpoel(i, j);	   // local point of element from connectivity matrix
			int istor = mesh1.esup2(ipoin, 0); // get storage number from reshuffled esup2
			mesh1.esup2(ipoin, 0) = istor - 1; // shift storage counter back one index to accept new element push back
			mesh1.esup1(istor - 1, 0) = i + 1; // push element number in to esup1
		}
	}
	// reshuffling pass 2 store in the starting index of the elementsupoin1 array for each element
	for (int i = 1; i < mesh1.npoin; i++)
	{
		mesh1.esup2(i, 0) = mesh1.esup2(i + 1, 0);
	}
	mesh1.esup2(mesh1.npoin, 0) = storage_req;
}
// endsub

// sub to generate elements surrounding elements data
void grid::pre_proc::set_esuel(grid::mesh &mesh1)
{
	assert(mesh1.func_count == 2); // function must run after setting up esup data structure
	mesh1.func_count = 3;		   // update func_count
	// init elsuel data structure
	// mesh1.esuel.resize(mesh1.nelem, mesh1.ntype); // number of elements x side for each element;
	mesh1.esuel.init(mesh1.nelem, mesh1.ntype);
	// loop over all elementsf
	for (int i = 0; i < mesh1.nelem; i++)
	{
		int ip1 = mesh1.inpoel(i, 0);
		int ip2 = mesh1.inpoel(i, 1);
		int ip3 = mesh1.inpoel(i, 2);
		int nesup1 = mesh1.esup2(ip1, 0) - mesh1.esup2(ip1 - 1, 0);
		int nesup2 = mesh1.esup2(ip2, 0) - mesh1.esup2(ip2 - 1, 0);
		int nesup3 = mesh1.esup2(ip3, 0) - mesh1.esup2(ip3 - 1, 0);
		int temp1 = 0;					 // flag to tell if sucessful push has happened
		for (int j = 0; j < nesup1; j++) // looping over all elements surrounding point p1 of element in
		{
			// get starting index for point p1 in the esup1 array from esup2 array
			int s_index = mesh1.esup2(ip1 - 1, 0);
			int e = mesh1.esup1(s_index + j, 0); // get element tag from esup1
			// std::cout<<e<<std::endl;
			if (e - 1 != i) // check whether the element is itself
			{
				for (int k = 0; k < mesh1.ntype; k++) // loop over all points of the element
				{
					if (mesh1.inpoel(e - 1, k) == ip2)
					{
						mesh1.esuel(i, 0) = e;
						temp1 = 1;
						// std::cout << inpoel(e - 1, k) << " " << ip2 << std::endl;
					}
				}
			}
		}
		if (temp1 != 1)
		{
			mesh1.esuel(i, 0) = 0;
		}
		int temp2 = 0;
		for (int j = 0; j < nesup2; j++) // looping over all elements surrounding point p1 of element in
		{
			// get starting index for point p1 in the esup1 array from esup2 array
			int s_index = mesh1.esup2(ip2 - 1, 0);
			int e = mesh1.esup1(s_index + j, 0); // get element tag from esup1
			if (e - 1 != i)						 // check whether the element is itself
			{
				for (int k = 0; k < mesh1.ntype; k++) // loop over all points of the element
				{
					if (mesh1.inpoel(e - 1, k) == ip3)
					{
						mesh1.esuel(i, 1) = e;
						temp2 = 1;
					}
				}
			}
		}
		if (temp2 != 1)
		{
			mesh1.esuel(i, 1) = 0;
		}
		int temp3 = 0;
		for (int j = 0; j < nesup3; j++) // looping over all elements surrounding point p3 of element in
		{
			// get starting index for point p3 in the esup1 array from esup2 array
			int s_index = mesh1.esup2(ip3 - 1, 0);
			int e = mesh1.esup1(s_index + j, 0); // get element tag from esup1
			if (e - 1 != i)						 // check whether the element is itself
			{
				for (int k = 0; k < mesh1.ntype; k++) // loop over all points of the element
				{
					if (mesh1.inpoel(e - 1, k) == ip1)
					{
						mesh1.esuel(i, 2) = e;
						temp3 = 1;
					}
				}
			}
		}
		if (temp3 != 1)
		{
			mesh1.esuel(i, 2) = 0;
		}
	}
}
// end sub

// subroutine to update and complete the bface data structuren initiated in the mesh constructor
void grid::pre_proc::set_bface(grid::mesh &mesh1)
{
	assert(mesh1.func_count == 3);		  // function must run after setting up elsuel
	mesh1.func_count = 4;				  // update function counter
	int bf = 0;							  // counter for boundary face number;
	for (int i = 0; i < mesh1.nelem; i++) // loop through all the elements
	{
		for (int j = 0; j < mesh1.ntype; j++) // loop over all the faces of the element
		{
			int E = mesh1.esuel(i, j);
			if (E == 0 and E - 1 < i) // this is a boundary cell and has boundar face
			{
				mesh1.bface(bf, 3) = bf + 1 + mesh1.nelem;
				int &p1 = mesh1.bface(bf, 0);
				int &p2 = mesh1.bface(bf, 1);							  // get local points
				int nesup1 = mesh1.esup2(p1, 0) - mesh1.esup2(p1 - 1, 0); // number of elements  surrounding point p1
				for (int r = 0; r < nesup1; r++)
				{
					int &s_index = mesh1.esup2(p1 - 1, 0); // get starting index for esup1 array
					int &e = mesh1.esup1(s_index + r, 0);
					for (int s = 0; s < mesh1.ntype; s++)
					{
						if (mesh1.inpoel(e - 1, s) == p2)
							mesh1.bface(bf, 2) = e;
					}
				}
				bf++;
			}
		}
	}
}
// endsub

// format of ith row of the intface structure (p1, p2, left element, right element)
// sub to compute and generate the interface connectivity matrix
void grid::pre_proc::set_intface(grid::mesh &mesh1)
{
	assert(mesh1.func_count == 4); // function must run after bface has been set
	mesh1.func_count = 5;		   // update function counter
	mesh1.nmaxface = (3 * mesh1.nelem + mesh1.nbface) / 2;
	mesh1.nintface = mesh1.nmaxface - mesh1.nbface;
	mesh1.intface.init(mesh1.nintface, 4); // init intface matrix container
	int m = 0;							   // face counter
	for (int i = 0; i < mesh1.nelem; i++)  // loop through all the elements
	{
		for (int j = 0; j < mesh1.ntype; j++) // loop over all faces of the element
		{
			int &e = mesh1.esuel(i, j); // get element surrouding current element
			if (e != 0 and e < i + 1)	// this is an internal cell
			{
				if (j != mesh1.ntype - 1)
				{
					mesh1.intface(m, 0) = mesh1.inpoel(i, j);
					mesh1.intface(m, 1) = mesh1.inpoel(i, j + 1);
					mesh1.intface(m, 2) = e;	 // left cell  i_e
					mesh1.intface(m, 3) = i + 1; // right cell j_e
				}
				else
				{
					mesh1.intface(m, 0) = mesh1.inpoel(i, 2);
					mesh1.intface(m, 1) = mesh1.inpoel(i, 0);
					mesh1.intface(m, 2) = e;
					mesh1.intface(m, 3) = i + 1;
				}
				m++;
 			}
		}
	}
}
// endsub


// sub to compute element data 
void grid::pre_proc::set_geoel(mesh &mesh1)
{
  assert(mesh1.func_count == 5);
  mesh1.func_count = 6;
	mesh1.geoel.init(mesh1.nelem, 15);
	for (int i = 0; i < mesh1.nelem; i++)
	{
		int &p1 = mesh1.inpoel(i, 0);
		int &p2 = mesh1.inpoel(i, 1);
		int &p3 = mesh1.inpoel(i, 2);
		double &x1 = mesh1.coords(p1 - 1, 0);
		double &y1 = mesh1.coords(p1 - 1, 1);
		double &x2 = mesh1.coords(p2 - 1, 0);
		double &y2 = mesh1.coords(p2 - 1, 1);
		double &x3 = mesh1.coords(p3 - 1, 0);
		double &y3 = mesh1.coords(p3 - 1, 1);
		mesh1.geoel(i, 0) = grid::pre_proc::el_jacobian(x1, x2, x3, y1, y2, y3)(0,0);
		mesh1.geoel(i, 1) = (x1 + x2 + x3) / 3; //xc  
		mesh1.geoel(i, 2) = (y1 + y2 + y3) / 3; //yc
		mesh1.geoel(i,3) = std::max({x1,x2,x3}) - std::min({x1,x2,x3}); // delta_x
		mesh1.geoel(i,4) = std::max({y1,y2,y3}) - std::min({y1,y2,y3}); //delta_y 
		mesh1.geoel(i,5) = mesh1.coords(p1-1,0)*0.5 + mesh1.coords(p2-1,0)*0.5; //guass points for domain integral 
		mesh1.geoel(i,6) = mesh1.coords(p1-1,1)*0.5 + mesh1.coords(p2-1,1)*0.5;
		mesh1.geoel(i,7) = mesh1.coords(p2-1,0)*0.5 + mesh1.coords(p3-1,0)*0.5;
		mesh1.geoel(i,8) = mesh1.coords(p2-1,1)*0.5 + mesh1.coords(p3-1,1)*0.5;
		mesh1.geoel(i,9) = mesh1.coords(p3-1,0)*0.5 + mesh1.coords(p1-1,0)*0.5;
		mesh1.geoel(i,10) = mesh1.coords(p3-1,1)*0.5 + mesh1.coords(p1-1,1)*0.5;
	}
}
// endsub

// sub tom compute the mass matrix componenets and store in initiated geoel array for DGP(1)
void grid::pre_proc::set_massMat(grid::mesh &mesh1)
{
	assert(mesh1.func_count == 6);
  	mesh1.func_count = 7;     // make sure geoel is inited
	for (int i = 0; i < mesh1.nelem; i++)//loop over all elements
	{
		// get cell data from current mesh object
		int &p1 = mesh1.inpoel(i, 0);
		int &p2 = mesh1.inpoel(i, 1);
		int &p3 = mesh1.inpoel(i, 2);
		double &x1 = mesh1.coords(p1 - 1, 0);
		double &x2 = mesh1.coords(p2 - 1, 0);
		double &x3 = mesh1.coords(p3 - 1, 0);
		double &y1 = mesh1.coords(p1 - 1, 1);
		double &y2 = mesh1.coords(p2 - 1, 1);
		double &y3 = mesh1.coords(p3 - 1, 1);
		// get centroids of the triangle
		double &xc = mesh1.geoel(i,1);
		double &yc = mesh1.geoel(i,2);
		// get delta_x and delta_y
		double &delta_x = mesh1.geoel(i,3);
		double &delta_y = mesh1.geoel(i,4);
		// get area of the triangle
		double &A = mesh1.geoel(i,0);
		mesh1.geoel(i,11) = 	2*A/(delta_x*delta_y)*(x1*x1/12 + x1*x2/12 + x1*x3/12 - x1*xc/3 + x2*x2/12 + x2*x3/12 - x2*xc/3 + x3*x3/12 - x3*xc/3 + xc*xc/2);
		mesh1.geoel(i,12) =  2*A/(delta_x*delta_y)*(x1*y1/12 + x1*y2/24 + x2*y1/24 + x1*y3/24 + x2*y2/12 + x3*y1/24 + x2*y3/24 + x3*y2/24 + x3*y3/12 - x1*yc/6- xc*y1/6 - x2*yc/6 - xc*y2/6 - x3*yc/6 - xc*y3/6 + xc*yc/2);
		mesh1.geoel(i,13) =  2*A/(delta_y*delta_y)*(y1*y1/12 + y1*y2/12 + y1*y3/12 - y1*yc/3 + y2*y2/12 + y2*y3/12 - y2*yc/3 + y3*y3/12 - y3*yc/3 + yc*yc/2);
	}
}
// endsub

// sub to generate internal face data needed for DG formulation
void grid::pre_proc::set_int_geoface(mesh &mesh1)
{
  assert(mesh1.func_count == 7);
  mesh1.func_count = 8;
	mesh1.int_geoface.init(mesh1.nintface, 7);
	for (int i = 0; i < mesh1.nintface; i++)
	{
		int &ip1 = mesh1.intface(i, 0); // 1st point on the face
		int &ip2 = mesh1.intface(i, 1); // second point on the face
		double &p1x = mesh1.coords(ip1 - 1, 0);
		double &p2x = mesh1.coords(ip2 - 1, 0);
		double &p1y = mesh1.coords(ip1 - 1, 1);
		double &p2y = mesh1.coords(ip2 - 1, 1);
		double mag = sqrt(std::pow((p2y-p1y),2) + std::pow((p2x-p1x),2));// push the components of Area weighted normal vectors to the geoface matrix
		mesh1.int_geoface(i, 0) = (p2y - p1y)/mag; //unit normal vector components x and y 
		mesh1.int_geoface(i, 1) = -1*(p2x-p1x)/mag;
		// compute the coords of the gauss points
		double E1 = -1 / sqrt(3);
		double E2 = 1 / sqrt(3);
		double gx1 = 0.5 * (1 - E1) * p1x + 0.5 * (1 + E1) * p2x;
		double gx2 = 0.5 * (1 - E2) * p1x + 0.5 * (1 + E2) * p2x;
		double gy1 = 0.5 * (1 - E1) * p1y + 0.5 * (1 + E1) * p2y;
		double gy2 = 0.5 * (1 - E2) * p1y + 0.5 * (1 + E2) * p2y;
		mesh1.int_geoface(i, 2) = gx1;
		mesh1.int_geoface(i, 3) = gy1;
		mesh1.int_geoface(i, 4) = gx2;
		mesh1.int_geoface(i, 5) = gy2;
		mesh1.int_geoface(i, 6) = len(p1x, p2x, p1y, p2y);
	}
}
// endsub

// sub to compute face data for boundary faces needed for rhsboun
void grid::pre_proc::set_boun_geoface(grid::mesh &mesh1)
{
  assert(mesh1.func_count == 8);
  mesh1.boun_geoface.init(mesh1.nbface,7);
  for(int i=0;i<mesh1.nbface;i++) //loop through all the boundary faces
  {
    int &p1 = mesh1.bface(i,0);
    int &p2 = mesh1.bface(i,1);
    double &p1x = mesh1.coords(p1-1,0);
    double &p1y = mesh1.coords(p1-1,1);
    double &p2x = mesh1.coords(p2-1,0);
    double &p2y = mesh1.coords(p2-1,1);
    double mag = sqrt(std::pow((p2y-p1y),2) + std::pow((p2x-p1x),2));   
    mesh1.boun_geoface(i,0) = p2y - p1y/mag; //unit normal vectors components x and yu
    mesh1.boun_geoface(i,1) = p1x - p2x/mag;
		double E1 = -1 / sqrt(3);
		double E2 = 1 / sqrt(3);
		double gx1 = 0.5 * (1 - E1) * p1x + 0.5 * (1 + E1) * p2x;
		double gx2 = 0.5 * (1 - E2) * p1x + 0.5 * (1 + E2) * p2x;
		double gy1 = 0.5 * (1 - E1) * p1y + 0.5 * (1 + E1) * p2y;
		double gy2 = 0.5 * (1 - E2) * p1y + 0.5 * (1 + E2) * p2y;
		mesh1.boun_geoface(i, 2) = gx1;
		mesh1.boun_geoface(i, 3) = gy1;
		mesh1.boun_geoface(i, 4) = gx2;
		mesh1.boun_geoface(i, 5) = gy2;
		mesh1.boun_geoface(i, 6) = len(p1x, p2x, p1y, p2y);
  }
}
// endsub




// sub to  to write out mesh in .vtk format for paraview to read and visualize
void grid::post_proc::writevtk_mesh(grid::mesh &mesh1, std::string file_name)
{
	std::ofstream file1(file_name);
	if (file1.is_open())
	{
		file1 << "# vtk DataFile Version 2.0\nfeflo\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS " << mesh1.npoin << " double  "
			  << "\n";
		for (int i = 0; i < mesh1.npoin; i++)
		{
			file1 << mesh1.coords(i, 0) << " " << mesh1.coords(i, 1) << " "
				  << "0 \n";
		}
		file1 << "CELLS " << mesh1.nelem << " " << mesh1.ntype * mesh1.nelem + mesh1.nelem << "\n";
		for (int k = 0; k < mesh1.nelem; k++)
		{
			file1 << mesh1.ntype << " " << mesh1.inpoel(k, 0) - 1 << " " << mesh1.inpoel(k, 1) - 1 << " " << mesh1.inpoel(k, 2) - 1 << "\n";
		}
		file1 << "CELL_TYPES " << mesh1.nelem << "\n";
		for (int k = 0; k < mesh1.nelem; k++)
		{
			file1 << mesh1.ntype + 2 << "\n"; // VTK uses flag 5 to identify a triangle surface
		}
	}
}
// endsub

// method to computes the pressure given values of properties
double EOS::perf_gas(matrix2d &cons_var)
{
	double pressure;
	pressure = (const_properties::gamma - 1) * cons_var(0, 0) * (cons_var(0, 0) - 0.5 / cons_var(0, 0) * (cons_var(1, 0) * cons_var(1, 0) + cons_var(2, 0) * cons_var(2, 0)));
	return pressure;
}

// this method computes the jacobian for a triangular element
matrix2d grid::pre_proc::el_jacobian(double &x1, double &x2, double &x3, double &y1, double &y2, double &y3)
{
  matrix2d el(2,1);
	el(0,0) = 0.5 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
  double x0 = 1*((x1*x1 + y1*y1)*(y2-y3)-(x2*x2 + y2*y2)*(y1-y3) + (x3*x3 + y3*y3)*(y1-y2))/(4*el(0,0));
  double y0 = -1*((x1*x1 + y1*y1)*(x2-x3)-(x2*x2 + y2*y2)*(x1-x3) + (x3*x3 + y3*y3)*(x1-x2))/(4*el(0,0));
  el(1,0) = len(x2,x0,y2,y0); 
  return el;
}

// this method computes the length of the 2d face
double len(double &x1, double &x2, double &y1, double &y2)
{
	return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
} 


//assemble the mesh
void grid::construct(grid::mesh &mesh1)
{
  grid::pre_proc::set_esup(mesh1); 
  grid::pre_proc::set_esuel(mesh1);
  grid::pre_proc::set_bface(mesh1);
  grid::pre_proc::set_intface(mesh1);
  grid::pre_proc::set_geoel(mesh1);
  grid::pre_proc::set_massMat(mesh1);
  grid::pre_proc::set_int_geoface(mesh1);
  grid::pre_proc::set_boun_geoface(mesh1);
}



