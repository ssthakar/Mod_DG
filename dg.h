#ifndef DG_H
#define DG_H
#include "mesh.h"
#include <stdexcept>

class soln//instantiate at the beginning of the simulation 
{
	public:
	double abstol;
	int max_iter;
	soln(grid::mesh &mesh1,std::string s1); //read in the control file and intialize the flowfield for the mesh
};

namespace DG
{
  void cons(grid::mesh &mesh1); //check if scheme is conservativeh
	void delta_T(grid::mesh &mesh1);
	void residual(grid::mesh &mesh1);
	bool isSolnConverged(grid::mesh &mesh1,soln &soln1);
	double det(matrix2d &coords); //function to get dererminant of from any three points;
	double circum_rad(matrix2d &coords); //function to get the radius of the circum-circle of a triangle
  double vel_sound(matrix2d &U); //get the speed of sound given the conservative variables vector 
	matrix2d Fx(matrix2d &state);// get x direction flux from conservative variables  (standard flux nothing fancy)
	matrix2d Fy(matrix2d &state);
	matrix2d U_at_poin(grid::mesh &mesh1, double &gx, double &gy, int i); // get solution at point  given the point co-ordinates
	matrix2d U_at_poin(grid::mesh &mesh1, double &gx, double &gy, int i,int e,double &nx,double &ny); //overload function to get state when we don't know if the cell is host or ghost
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
		matrix2d W_amp;// store in the wave amplitude
		matrix2d avg_flux; // central averaged flux in cell normal direction
	public:
		matrix2d lamda;
		matrix2d intface_flux;
	private:	
		void set_lamda();  // compute and store in eign values of the problem
		void set_W_amp();  // compute and store in wave strength
		void set_diss_flux();
	public:
		void compute_req(matrix2d &Ul, matrix2d &Ur, double &nx, double &ny); // constructor computes the averaged vars, stores in nx and ny
		void compute_flux();
		matrix2d flux_intface();
	};
}


namespace ddt // time marching methods for local cells
{
	double local_ts1(grid::mesh &mesh1, int &i);
	double local_ts(grid::mesh &mesh1,int &i); // method to calculate the allowed delta T for all cells (local time step)
	namespace RK3 //TVD Runge Kutta 3 stage 
	{
		void RK3_outer(grid::mesh &mesh1, soln &soln1);// outer function 
		void RK_s1(grid::mesh &mesh1); //first stage, also forward Euler
		void RK_s2(grid::mesh &mesh1); //second stage 
		void RK_s3(grid::mesh &mesh1); //third stage
	};
};

#endif
