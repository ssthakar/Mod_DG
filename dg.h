#ifndef DG_H
#define DG_H 
#include "mesh.h"

namespace DG
{
	matrix2d get_state(grid::mesh &mesh1,int &i); //get conservative state vars of a particular elem
	void init_field(grid::mesh &mesh1);
	void rhsboun_bface(grid::mesh &mesh1);// calculate the contribution of the boundary integral from boundary faces
	void rhsboun_iface(grid::mesh &mesh1); // calculate the contribution from boundary integral for internal faces
	void rhsdomn(grid::mesh &mesh1); //calculate the contribution of the boundary integral to rhsel
};	

namespace FDS
{
	class RoeFlux
	{
		public:
			double nx;
			double ny;
			matrix2d dissFlux; //disspative flux substracted from the CD
			RoeFlux(matrix2d &Ul, matrix2d &Ur, double &nx, double &ny);//constructor computes the averaged vars, stores in nx and ny  
			matrix2d avg; // store in Roes averaged vars, computed with object instantiation through the constructor
			matrix2d lamda; //store in the eigen vals 
			matrix2d W_amp; //store in the wave amplitude
			void set_lamda(); //compute and store in eign values of the problem
			void set_W_amp(); //compute and store in wave strength
			void flux_out(); //finally compute and store in dissflux 
	};

}
namespace ddt  //time marching methods for local cells 
{
	double calc_deltaT(grid::mesh &mesh1); //method to calculate the allowed delta T for all cells (local time step)
	namespace explct
	{
		matrix2d fwd_euler(grid::mesh &mesh1, double &delta);
		matrix2d RK3(grid::mesh &mesh1,double &delta_t); // TVD runge Kutta 3rd order takes reference to mesh object and timestep 
	};
};

#endif 

