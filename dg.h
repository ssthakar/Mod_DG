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
  matrix2d fv_state(grid::mesh &mesh1, int i); //gives the finite volume state for any cell i
  void fv_rhsboun_iface(grid::mesh &mesh1); //sub to compute the bounday contribution for finite volume only internal faces 
  matrix2d fv_state(grid::mesh &mesh1, int i,int e, double &nx, double &ny); //gives the finite volume state for any cell i
  void fv_rhsboun_bface(grid::mesh &mesh1); //sub to compute the boundary contribution for finitie volume only boundary faces 


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


namespace FVS
{
  matrix2d Van_leer(grid::mesh &mesh1,matrix2d &Ul, matrix2d &Ur, double &nx, double &ny); //Van Leer FLux vector split
}


namespace ddt // time marching methods for local cells
{
  double calc_ts(grid::mesh &mesh1); //sub to calculate global timestep
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
