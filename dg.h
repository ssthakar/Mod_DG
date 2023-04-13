#ifndef DG_H
#define DG_H
#include "mesh.h"

namespace DG
{
  double vel_sound(matrix2d &U); //get the speed of sound given the conservative variables vector 
	matrix2d Fx(matrix2d &state);// get x direction flux from conservative variables  (standard flux nothing fancy)
	matrix2d Fy(matrix2d &state);
	matrix2d U_at_poin(grid::mesh &mesh1, double &gx, double &gy, int i); // get solution at point  given the point co-ordinates
	matrix2d get_state(grid::mesh &mesh1, int &i);						  // get conservative state vars of a particular elem
	void init_field(grid::mesh &mesh1);
	void rhsboun_bface(grid::mesh &mesh1); // calculate the contribution of the boundary integral from boundary faces
	void rhsboun_iface(grid::mesh &mesh1); // calculate the contribution from boundary integral for internal faces
	void rhsdomn(grid::mesh &mesh1);	   // calculate the contribution of the boundary integral to rhsel
};

namespace FDS
{

	// Roes averaged flux difference split with Harten's conrrection for the eigen values
	class RoeFlux
	{

	private:
		double Nx;
		double Ny;
		double delta_rho;
		double deltaVx;
		double deltaVy;
		double deltaP;
		double deltaVn;
		matrix2d dissFlux; // disspative flux substracted from the CD
		matrix2d avg;	   // store in Roes averaged vars, computed with object instantiation through the constructor
										 // computes the central averaged flux in the normal direction
		matrix2d lamda;	   // st  ore in the eigen vals
		matrix2d W_amp;// store in the wave amplitude
		matrix2d avg_flux; // central averaged flux in cell normal direction
	public:
		matrix2d intface_flux;
	private:	
		void set_lamda();  // compute and store in eign values of the problem
		void set_W_amp();  // compute and store in wave strength
		void set_diss_flux();
	public:
		RoeFlux();
		void compute_req(matrix2d &Ul, matrix2d &Ur, double &nx, double &ny); // constructor computes the averaged vars, stores in nx and ny
		void compute_flux();
		matrix2d flux_intface();
	};
}
namespace FVS
{
	// Van Leers quadtratic function vector split
	class VanLeer
	{
	};
}



namespace ddt // time marching methods for local cells
{
	void local_ts(grid::mesh &mesh1,int &i); // method to calculate the allowed delta T for all cells (local time step)
	namespace explct
	{
		matrix2d fwd_euler(grid::mesh &mesh1, double &delta);
		matrix2d RK3(grid::mesh &mesh1, double &delta_t); // TVD runge Kutta 3rd order takes reference to mesh object and timestep
	};
};

#endif
