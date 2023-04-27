#include "mesh.h"

namespace FDS
{

	// Roes averaged flux difference split with Harten's conrrection for the eigen values
	class RoeFlux
	{

	private:
    matrix2d left;
    matrix2d right;
    //matrix2d tengential;
    matrix2d state;
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
	  matrix2d Roe_matx; //matrix to store in right eigen vectors
    double neg_eig_vals; //number of left moving waves
    double ul,vl,ur,vr,Hl,Hr,pl,pr,rhoR,rhoL,re,le;
  public:
    double fl;
		matrix2d lamda; //eigen values of the PDE
		matrix2d intface_flux; //return the interface flux at any face 
	private:
    matrix2d harten_corrector(); //adds the Harten correction
		void set_lamda();  // compute and store in eign values of the problem
		void set_vec(); //computes the row matrix from Roes averages and other data
    void set_W_amp();  // compute and store in wave strength
		void set_diss_flux();
	public:
		void compute_req(matrix2d &Ul, matrix2d &Ur, double &nx, double &ny); // constructor computes the averaged vars, stores in nx and ny
		void compute_flux();
		matrix2d flux_intface();
	};
}


/*---------------------------------------------------------------------------------------------------------------------------------------------------------------------------*/
//  imppmentation of Roes averaged flux, constructs flux object and has a method that returns the dissipative flux vector to be used in the rhsboun_iface subroutine
// instantiante object outside loop and then compute flux for every loop iteration
void FDS::RoeFlux::compute_req(matrix2d &Ul, matrix2d &Ur, double &nx, double &ny) // compute function, computes the required data for other setter functions/kinda like a constructor
{

  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT); 
  //std::cout<<"flux req called"<<std::endl;
  assert(Ul.rows() == 4 and Ur.rows() == 4); // this implementation only works for 2d problems where the conservative vectors size is 4
  // compute and store in roe-averaged values //stores in other data that is needed for the caluclation of the dissipative flux
  left = Ul;
  right = Ur;
  Nx = nx;
  Ny = ny; // store in normal vectors for further usage
  rhoL = Ul(0,0);
  rhoR = Ur(0,0);
  double invRhoL = 1/Ul(0,0);
  double invRhoR = 1/Ur(0,0);
  assert(Ul.cols() == 1 and Ur.cols() == 1);
  double Pl = EOS::perf_gas(Ul); // compute left and right pressure states using EOS
  double Pr = EOS::perf_gas(Ur);
  deltaP = Pr - Pl; // compute delta values needed to compute wave amp
  deltaVx = Ur(1, 0)/Ur(0,0) - Ul(1, 0)/Ul(0,0);
  deltaVy = Ur(2, 0)/Ur(0,0) - Ul(2, 0)/Ul(0,0);
  delta_rho = Ur(0, 0) - Ul(0, 0);
  deltaVn = (Ur(1, 0)/Ur(0,0) * Nx + Ur(2, 0)/Ur(0,0) * ny) - (Ul(1, 0)/Ul(0,0) * Nx + Ul(2, 0)/Ul(0,0) * Ny);
  avg.init(5, 1);      // init container , stores in roe agveraged vars
  avg_flux.init(4, 1); // init container to store in central ave flux
  //for left cell
  ul = (Ul(1,0)*nx + Ul(2,0)*ny)*invRhoL; //normal velocity, this is the x direction split equivalent of x 
  vl = (-Ul(2,0)*ny + Ul(1,0)*nx)*invRhoL; //tangential velocity, this is the y vel for a x direction like split but in the normal direction
  pl = (const_properties::gamma-1)*(Ul(3,0)-0.5*(ul*ul + vl*vl)); //get the pressure for the left side
  Hl = invRhoL*const_properties::cp*pl + 0.5*(ul*ul + vl*vl);
  le = Ul(3,0);
  re = Ur(3,0);

  //for right cell 
  ur = (Ur(1,0)*nx + Ur(2,0)*ny)*invRhoR; //normal velocity, this is the x direction split equivalent of x 
  vr = (-Ul(2,0)*ny + Ul(1,0)*nx)*invRhoR; //tangential velocity, this is the y vel for a x direction like split but in the normal direction
  pr = (const_properties::gamma-1)*(Ul(3,0)-0.5*(ur*ur + vr*vr)); //get the pressure for the left side
  Hr = invRhoR*const_properties::cp*pr + 0.5*(ur*ur + vr*vr);

  //std::cout<<"break_b4 delta"<<std::endl;
  //print2Term(Ul);
  //std::cout<<"break_a4_detla"<<std::endl;
  double Rij = sqrt(Ur(0, 0) / Ul(0, 0));

  std::cout<<"befor Roes average computation"<<std::endl;
  print2Term(Ul);
  avg(0, 0) = Rij * Ul(0, 0);                                                                                            // rhoij
  avg(1, 0) = (Rij*ur + ul)/(1.0 + Rij);                                             // Vxij
  avg(2, 0) = (Rij*vr + vl )/(1.0 + Rij);                                             // Vyij
  avg(3, 0) = (Rij*Hr + Hl)/(Rij + 1.0);                               // Hij
  std::cout<<"break"<<(avg(3,0) - 0.5*(avg(1,0)*avg(1,0) + avg(2,0)*avg(2,0)))<<std::endl;
  avg(4, 0) = sqrt((const_properties::gamma - 1.0)*(avg(3,0) - 0.5*(avg(1,0)*avg(1,0) + avg(2,0)*avg(2,0)))); // Cij speed of sound
  std::cout<<"break2"<<std::endl;
  //print2Term(avg);
  //std::cout<<"before averaga flux computation"<<std::endl; 
  // compute averaged flux from right and left states at athe cell interface
  avg_flux(0, 0) = (Ul(1, 0) * nx + Ul(2, 0) * ny + Ur(1, 0) * nx + Ur(2, 0) * ny);
  avg_flux(1, 0) = ((Ul(1, 0) * Ul(1, 0) + Pl) * nx + Ul(1, 0) * Ul(2, 0) / Ul(0, 0) * ny + (Ur(1, 0) * Ur(1, 0) + Pr) * nx + Ur(1, 0) * Ur(2, 0) / Ur(0, 0) * ny);
  avg_flux(2, 0) = (Ul(1, 0) * Ul(2, 0) / Ul(0, 0) * nx + (Ul(2, 0) * Ul(2, 0) / Ul(0, 0) + Pl) * ny + Ur(1, 0) * Ur(2, 0) / Ur(0, 0) * nx + (Ur(2, 0) * Ur(2, 0) / Ur(0, 0) + Pr) * ny);
  avg_flux(3, 0) = ((Ul(3,0) + Pl)*Ul(1,0)/Ul(0, 0)*nx + (Ul(3, 0) + Pl) * Ul(2, 0) / Ul(0, 0) * ny + (Ur(3, 0) + Pr) * Ur(1, 0) / Ur(0, 0) * nx + (Ur(3, 0) + Pr) * Ur(2, 0) / Ur(0, 0) * ny);
  //std::cout<<"flux req compute complete"<<std::endl;
}

// compute the eigen values of the A matrix
void FDS::RoeFlux::set_lamda()
{
  //std::cout<<"set lam called"<<std::endl;
  // first two lamdas are the same
  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT); 
  lamda.init(4, 1); // init lamdada
  lamda(0, 0) = avg(1, 0) - avg(4,0);
  lamda(1,0) =  avg(1, 0);
  lamda(2, 0) = avg(1,0);
  lamda(3, 0) = avg(1,0) + avg(4, 0);
  neg_eig_vals = 0;
  for(int i=0;i<lamda.rows();i++)
  {
    if(lamda(i,0)<0)
      neg_eig_vals++; //count the number of right moving eigen vectors
  }
  //std::cout<<"set lam complete"<<std::endl;
}

void FDS::RoeFlux::set_vec()
{
  Roe_matx.init(4,4);
  Roe_matx(0,0) = Roe_matx(0,1) = Roe_matx(0,3) = 1;
  Roe_matx(0,2) = 0;
  Roe_matx(1,0) = avg(1,0) - avg(4,0);
  Roe_matx(1,1)  = avg(1,0);
  Roe_matx(1,2) = 0;
  Roe_matx(1,3) = avg(1,0)+avg(4,0);
  Roe_matx(2,0) = avg(2,0);
  Roe_matx(2,1) = avg(2,0);
  Roe_matx(2,2) = 1;
  Roe_matx(2,3) = avg(2,0);
  Roe_matx(3,0) = avg(3,0) - avg(1,0)*avg(4,0);
  Roe_matx(3,1) = 0.5*(std::pow(avg(1,0),2)+std::pow(avg(2,0),2));
  Roe_matx(3,2) = avg(2,0);
  Roe_matx(3,3) = avg(3,0)+avg(4,0)*avg(1,0);
}


// compute the wave amplitudes and store them
void FDS::RoeFlux::set_W_amp()
{
  double dq_rho = delta_rho;
  double dq_u = rhoR*ur-rhoL*ul;
  double dq_v = rhoR*vr-rhoL*vl;
  double dq_e = rhoR*re - rhoL*le;
  //std::cout<<"wamp called"<<std::endl;
  W_amp.init(4, 1); // inti container and compute store in wave amplitudes
  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT); 
  W_amp(2, 0) = dq_v - avg(2,0)*dq_rho;
  W_amp(1,0) = (const_properties::gamma-1)*(dq_rho * (avg(3,0) - avg(1,0)*avg(1,0)) + avg(1,0)*dq_u - dq_e + (dq_v - avg(2,0)*dq_rho*avg(2,0)))/(avg(4,0)*avg(4,0));
  W_amp(0, 0) = (dq_rho * (avg(1,0) + avg(4,0)) - dq_u - avg(4,0)* W_amp(1,0)) / (2.0*avg(4,0));
  W_amp(3, 0) = dq_rho-(W_amp(0,0)+W_amp(1,0));
  //std::cout<<"wamp cpmlt"<<std::endl;
}
// compute the dissipative flux that is substracted from the central averaged flux at a gauss point
void FDS::RoeFlux::set_diss_flux()
{
  //std::cout<<"diss flux cal"<<std::endl;
  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT); 
  double r = avg(0, 0) / (2 * avg(4, 0));
  dissFlux.init(4, 1);
  // this is the dissipative flux that is substracted from the central average
  double a0 = harten_corrector()(0,0)* W_amp(0,0); //lamda1(0,0) x  wamp1(0,0)
  double a1 = harten_corrector()(1,0)* W_amp(1,0); //lamda2(0,0 ) x wamp2(1,0)
  double a2 = harten_corrector()(2,0)* W_amp(2,0); //lamda3(1,0) x wamp3(2,0)
  double a3 = harten_corrector()(3,0)* W_amp(3,0); //lamda4(3,0) x wamp4(3,0)
  dissFlux(0, 0) = (a0)*Roe_matx(0,0) + a2*Roe_matx(0,2) + a3*Roe_matx(0,3);
  dissFlux(1, 0) = a0*Roe_matx(1,0) + a1*Roe_matx(1,1) + a2*Roe_matx(1,2) + a3*Roe_matx(1,3);
  dissFlux(2, 0) = a0*Roe_matx(2,0) + a1*Roe_matx(2,1) + a2*Roe_matx(2,2) + a3*Roe_matx(2,3);
  dissFlux(3, 0) = a0*Roe_matx(3,0) + a1*Roe_matx(3,1) + a2*Roe_matx(3,2) + a3*Roe_matx(3,3);
  // std::cout<<"dis flx compt"<<std::endl;
}
// get interface flux by calling this function in the rhsboun subs
void FDS::RoeFlux::compute_flux()
{
  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT); 
  //std::cout<<"compute flux called"<<std::endl;
  set_lamda();
  set_vec();
  set_W_amp();
  set_diss_flux();
	intface_flux.init(4,1);
  intface_flux(0, 0) = (avg_flux(0, 0) - dissFlux(0, 0))*0.5;
  intface_flux(1, 0) = (avg_flux(1, 0) - dissFlux(1, 0))*0.5;
  intface_flux(2, 0) = (avg_flux(2, 0) - dissFlux(2, 0))*0.5;
  intface_flux(3, 0) = (avg_flux(3, 0) - dissFlux(3, 0))*0.5;
  //std::cout<<"compute flux complete"<<std::endl;
}

//standard entropy fix for Roes flux
matrix2d FDS::RoeFlux::harten_corrector()
{
  matrix2d H(4,1);
  for(int i=0;i<4;i++) 
  {
    if(std::fabs(lamda(i,0))<deltaVn)
    {
      H(i,0)= (std::pow(lamda(i,0),2) + std::pow(deltaVn,2))/(2*deltaVn);
    }
    else
    {
      H(i,0) = std::fabs(lamda(i,0));
    }
  }

  return H;
}
// end class definition
