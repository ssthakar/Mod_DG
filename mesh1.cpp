#include "mesh.h"
#include "dg.h"
#include <cmath>
#include <sstream>

// TODO think of a better implementation using mmatrix class
// stores data in a buffer of vectors of vectors
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
				if (std::stringstream(word) >> val) // convers strings to doubles as well, TODO find better method to do this
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

	std::vector<std::vector<double>> v = reader::readv(s1);
	std::vector<std::vector<double>> c = reader::readv(s2);
	// get all the int members of the class
	ndimn = v[5][0];
	ntype = v[5][1];
	nelem = v[7][0];
	npoin = v[7][1];
	nbface = v[7][2];
	neqns = ndimn+2;
	ndegr = (c[0][0] + 1)*(c[0][0]+2)*0.5;
	// for loop populate intpoel matrix

	inpoel.resize(nelem, ntype); // init and give size to inpoel
	for (int i = 0; i < nelem; i++)
	{
		for (int j = 0; j < ntype; j++)
		{
			inpoel(i, j) = v[i + 9][j + 1];
		}
	}
	coords.resize(npoin, ndimn); // init and give size to coords matrix
	for (int i = 0; i < npoin; i++)
	{
		for (int j = 0; j < ndimn; j++)
		{
			coords(i, j) =v[i + 10 + nelem][j + 1];
		}
	}
	bface.resize(nbface,3);//init bface matrix and give size 
	for(int i=0;i<nbface;i++)
	{
		for(int j=0;j<3;j++)
		{
			bface(i,j) = v[i+nelem+2*npoin+12][j+1];
		}
	}
}

//sub to generate elemenst surrounding points linked list
void grid::set_esup(mesh &mesh1)
{
	mesh1.esup2.resize(mesh1.npoin + 1, 1); // init matrix to store in index 
	// this pass counts number of elements surrounding every point and stores that number in +1 index corresponding to that particular point
	// for eg if point 9 has 6 elements around it, then
	//  the number stored in index 10(for c++ it will be index 9) will be the number 6 after the first element pass
	//
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
	mesh1.esup1.resize(storage_req, 1); // init esup1 based on storage req
	// element pass 2, storage
	for (int i = 0; i < mesh1.nelem; i++) // loop over all elements
	{
		for (int j = 0; j < mesh1.ntype; j++) // loop over all nodes of the element
		{
			int ipoin = mesh1.inpoel(i, j);	 // local point of element from connectivity matrix
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

//sub to generate elements surrounding elements data
void grid::set_esuel(mesh &mesh1)
{
	// subroutine to generalte element surrounding element data structure
	// init elsuel data structure
	mesh1.esuel.resize(mesh1.nelem, mesh1.ntype); // number of elements x side for each element;
	// loop over all elementsf
	for (int i = 0; i < mesh1.nelem; i++)
	{
		int ip1 = mesh1.inpoel(i, 0);
		int ip2 = mesh1.inpoel(i, 1);
		int ip3 = mesh1.inpoel(i, 2);
		int nesup1 = mesh1.esup2(ip1, 0) - mesh1.esup2(ip1 - 1, 0);
		int nesup2 = mesh1.esup2(ip2, 0) - mesh1.esup2(ip2 - 1, 0);
		int nesup3 = mesh1.esup2(ip3, 0) - mesh1.esup2(ip3 - 1, 0);
		int temp1 = 0; //flag to tell if sucessful push has happened
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
						//std::cout << inpoel(e - 1, k) << " " << ip2 << std::endl;
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
			if (e - 1 != i)				   // check whether the element is itself
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
			if (e - 1 != i)				   // check whether the element is itself
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

// subroutine to generate the inter face connectivity matrix 

	//FORMAT FOR INTFACE
	// (ip1 | ip2 |cell on left | cell on right| face flag)
	// 0 for internal face
	// 2 for inlet or outlet face 
	// 4 for wall face  
	//Note on convetion
	/*

	_______p2
  \     /\
	 \ L /  \
		\	/ R  \
		 p1_____\

	*/ 



void grid::set_intfafce(mesh &mesh1)
{
	mesh1.nmaxface = (3 * mesh1.nelem + mesh1.nbface) / 2;
	std::cout << " max faces is " << mesh1.nmaxface << std::endl;
	mesh1.intface.resize(mesh1.nmaxface, 5); // init intface matrix
	int nface = 0;				 // counter for periodic boundary
	int m = 0; // face counter
	int bf = 0;;
	// subroutinea to populate intface matrix, using periodic boundary condition
	for (int i = 0; i < mesh1.nelem; i++) //loop through all the elements 
	{
		for (int j = 0; j < mesh1.ntype; j++) //loop over all faces of the element 
		{
			int e = mesh1.esuel(i, j);
			if(e == 0 and e-1<i) //this is a boundary cell
			{
				nface++;
				mesh1.intface(m,0) = mesh1.bface(bf,0); //push boundary face data in ip1 
				mesh1.intface(m,1) = mesh1.bface(bf,1);// ip2 
				mesh1.intface(m,4) = mesh1.bface(bf,2);// boundary type flag
				mesh1.intface(m,3) = mesh1.nelem+nface;
				int ip1 = mesh1.bface(bf,0);
				int ip2 = mesh1.bface(bf,1);
				int nesup1 = mesh1.esup2(ip1, 0) - mesh1.esup2(ip1 - 1, 0);
				for(int r=0;r<nesup1;r++) //loop over all elements surrounding point p1
				{
					int s_index = mesh1.esup2(ip1-1,0); //get starting index
					int e = mesh1.esup1(s_index+r,0); //get element tag
					for(int s=0;s<mesh1.ntype;s++) //loop through all the points in the element 
					{
						if(mesh1.inpoel(e-1,s) == ip2)
						mesh1.intface(m,2) = e;
					}
				}
				m++;
				bf++;
			}
			if (e != 0 && e - 1 < i) // this is an internal cell
			{
				if (j != 2)
				{
					mesh1.intface(m, 0) = mesh1.inpoel(i, j);
					mesh1.intface(m, 1) = mesh1.inpoel(i, j + 1);
					mesh1.intface(m, 3) = e;
					mesh1.intface(m, 2) = i + 1;
					mesh1.intface(m,4) = 0; //flag for internal face 
					m++;
				}
				else
				{
					mesh1.intface(m, 0) = mesh1.inpoel(i, 2);
					mesh1.intface(m, 1) = mesh1.inpoel(i, 1);
					mesh1.intface(m, 3) = e;
					mesh1.intface(m, 2) = i + 1;
					mesh1.intface(m,4) = 0; //flag for internal face 
					m++;
				}
			}
		}
	} 
}


//sub to generate face data needed for DG formulation
void grid::set_geoface(mesh &mesh1)
{
	//format for geoface data structure (nx ny G1x G1y G2x G2y)

	// generate geoface data structure using intface, intpoel, coords 
	mesh1.geoface.resize(mesh1.nmaxface,6);
	for(int i=0;i<mesh1.nmaxface;i++)
	{
		double nx,ny;// components of the outward area normal vector
		int ip1 = mesh1.intface(i,0); //1st point on the face
		int ip2 = mesh1.intface(i,1); //second point on the face 
		double p1x = mesh1.coords(ip1-1,0);
		double p2x = mesh1.coords(ip2-1,0);
		double p1y = mesh1.coords(ip1-1,1);
		double p2y = mesh1.coords(ip2-1,1);
		nx = p2y - p1y;
		ny = -1*(p2x-p1x);
		// push the components of Area weighted normal vectors to the geoface matrix 
		mesh1.geoface(i,0) = nx;
		mesh1.geoface(i,1) = ny;
		

		// Note for two points the weights are 1 
		//calculate gauss quadrature points needed for rhs integral, for now i am using 2, but later can change it to user defined
		//transformation to get coords for the gauss points 
		double E1  = -1/sqrt(3);
		double E2 = 1/sqrt(3);
		double gx1 = 0.5*(1-E1)*p1x + 0.5*(1+E1)*p2x;
		double gx2 = 0.5*(1-E2)*p1x + 0.5*(1+E2)*p2x;
		double gy1 = 0.5*(1-E1)*p1y + 0.5*(1+E1)*p2y;
		double gy2 = 0.5*(1-E2)*p1y + 0.5*(1+E2)*p2y;
		mesh1.geoface(i,2) = gx1;
		mesh1.geoface(i,3) = gy1;
		mesh1.geoface(i,4) = gx2;
		mesh1.geoface(i,5) = gy2;
	}
}


// sub to compute element data like shape functions etc
void grid::set_geoel(mesh &mesh1)
{

}


//computes the pressure given values of properties
double EOS::perf_gas(double rho, double E, double u, double v)
{
	double pressure;
	pressure = (const_properties::gamma-1)*rho*(E - 0.5*sqrt(u*u + v*v));
	return pressure;
}



