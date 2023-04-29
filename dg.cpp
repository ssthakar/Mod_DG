#include "dg.h"
#include "mesh.h"
#include "vmatrix2.h"
#include <cmath>

//constructor for the solution class
soln::soln(grid::mesh &mesh1, std::string s1)
{
	mesh1.RK_stor.init(mesh1.nelem,mesh1.neqns,mesh1.ndegr);
	DG::init_field(mesh1); //intialize the flowfield
	std::vector<std::vector<double>> c = reader::readv(s1); //read in the control file 
	max_iter = c[9][0];
	abstol = c[8][0];

}


//subroutine to intialize the flow field 
void DG::init_field(grid::mesh &mesh1)
{
  mesh1.fvRkstor.init(mesh1.nelem,mesh1.neqns);
  mesh1.fvunkel.init(mesh1.nelem,mesh1.neqns);
  mesh1.fvrhsel.init(mesh1.nelem,mesh1.neqns);
  mesh1.Ltimestep.resize(mesh1.nelem,0);
	mesh1.unkel.init(mesh1.nelem,mesh1.neqns,mesh1.ndegr); //initialize the flow field container 
	mesh1.rhsel.init(mesh1.nelem,mesh1.neqns,mesh1.ndegr); //initialize the RHS vector  container 
	mesh1.res_vec.init(mesh1.neqns,1); //intialize the residual vector
	for(int i=0;i<mesh1.nelem;i++) //loop through all elements  
	{
		for(int m = 0;m<mesh1.neqns;m++)
		{
      mesh1.fvunkel(i,m) = mesh1.U_infty(m,0);
			mesh1.unkel(i,m,0) = mesh1.U_infty(m,0); //pass in free stram values to first coeffcient, assuming no variation of flow field properties in the cell, the gradients will therefore be 0;
		}
	}
}


//function to return the velocity of sound from state at any given point
double DG::vel_sound(matrix2d &Ul)
{
  double eps1 = const_properties::eps;
  //std::cout<<"Ul_from vel"<<std::endl;
  //print2Term(Ul);
  double ci2 = (const_properties::gamma*EOS::perf_gas(Ul)/Ul(0,0)); //return the speed of sound for a given conservative variable vector
  if(ci2<0)
  {
    //std::cout<<"vel"<<ci2<<std::endl;
    //print2Term(Ul);
  }
  //return sqrt(std::max(ci2,eps1));
  return sqrt(ci2);
}

matrix2d DG::U_at_poin(grid::mesh &mesh1, double &gx, double &gy, int i)
{
	matrix2d state;
  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT); 
	state.init(4,1);
  for(int k=0;k<mesh1.neqns;k++)
  {
	  state(k, 0) = mesh1.unkel(i, k, 0) + mesh1.unkel(i, k, 1) * (gx - mesh1.geoel(i, 1))/mesh1.geoel(i,3) + mesh1.unkel(i, k, 2) * (gy - mesh1.geoel(i, 2))/mesh1.geoel(i,4);
  }
  return state;
}

//overloaded function DG::U_at_poin to get state at point when its unclear whether the cell is ahost or ghost, mostly used for local time step calculation 
matrix2d DG::U_at_poin(grid::mesh &mesh1, double &gx, double &gy, int i,int e, double &nx, double &ny) // int i is the value esuel gives, e is the value of the current esuel row being looped through
{
  //std::cout<<"U_at poin overload called"<<std::endl;
	matrix2d state(4,1);
	matrix2d host = DG::U_at_poin(mesh1,gx,gy,e); // state at host cell 
	if(i==-2) //this is a wall cell and therefore will use mirror
	{
    //std::cout<<"overload case -2"<<std::endl;
    state(0,0) = host(0,0); //density
    state(3,0) = host(3,0);//energy
		// compute velocity using mirror images
		state(1,0) = host(1,0) - 2*(host(1,0)*nx + host(2,0)*ny)*nx;
    state(2,0) = host(2,0) - 2*(host(1,0)*nx + host(2,0)*ny)*ny;	
	}
	else if(i == -4) // this is a free stream cell 
	{
    //std::cout<<"overload case -4"<<std::endl;
		state = mesh1.U_infty;
	}
	else
	{
    for(int k=0;k<mesh1.neqns;k++)
    {

      //std::cout<<"normal call"<<std::endl;
		  state(k, 0) = mesh1.unkel(i, k, 0) + mesh1.unkel(i, k, 1) * (gx - mesh1.geoel(i, 1))/mesh1.geoel(i,3) + mesh1.unkel(i, k, 2) * (gy - mesh1.geoel(i, 2))/mesh1.geoel(i,4);
    }
	}	
	return state;
}


// function to get x-direction flux vector from state
matrix2d DG::Fx(matrix2d &U) // used :::in domain integral only
{
  matrix2d fx(4, 1);
  fx(0, 0) = U(1, 0); // rho U
  fx(1, 0) = U(1, 0) * U(1, 0) / U(0, 0) + EOS::perf_gas(U);
  fx(2, 0) = U(1, 0) * U(2, 0) / U(0, 0);
  fx(3, 0) = (U(3, 0) + EOS::perf_gas(U)) * U(1, 0) / U(0, 0);
  return fx;
}

// function to get y-direction flux vector from state
matrix2d DG::Fy(matrix2d &U) // used in domain integral only
{
  matrix2d fy(4, 1);
  fy(0, 0) = U(2, 0); // rho U
  fy(1, 0) = U(1, 0) * U(2, 0) / U(0, 0);
  fy(2, 0) = U(2, 0) * U(2, 0) / U(0, 0) + EOS::perf_gas(U);
  fy(3, 0) = (U(3, 0) + EOS::perf_gas(U)) * U(2, 0) / U(0, 0);
  return fy;
}

//sub to compute and push the contribution of the domain integral to rhs 
void DG::rhsdomn(grid::mesh &mesh1) // pass reference to the mesh object
{
  for (int i = 0; i < mesh1.nelem; i++) // loop over all elems
  {
    for (int j = 0; j < mesh1.ngauss_domn; j++) // loop over all gauss points
    {
      matrix2d U = DG::U_at_poin(mesh1, mesh1.geoel(i, 2*j+5), mesh1.geoel(i, 2*j+6), i);
      matrix2d shit = Fx(U);
      for(int l=0;l<mesh1.neqns;l++)
      {
        mesh1.rhsel(i, l, 1) = mesh1.rhsel(i, l, 1) - mesh1.geoel(i,0)*mesh1.domweight*DG::Fx(U)(l,0)*mesh1.geoel(i,1)/mesh1.geoel(i,3); //DG::Fx and DG::Fy are functions that return the x and y flux vector 
        mesh1.rhsel(i, l, 2) = mesh1.rhsel(i, l, 2) - mesh1.geoel(i,0)*mesh1.domweight*DG::Fy(U)(l,0)*mesh1.geoel(i,2)/mesh1.geoel(i,4);
      }
    }
  }
}
//ensdub

// sub to push the contribution of the boundary integral to rhsel from bface only
void DG::rhsboun_bface(grid::mesh &mesh1)
{
  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT); 
  assert(mesh1.unkel.size() > 1);        // make sure intialization of the flow field is done and the storage container is resized
  for (int i = 0; i < mesh1.nbface; i++) // loop over all boundary faces
  {
    int &le = mesh1.bface(i,2); //host element 
    int &re = mesh1.bface(i,3); //ghost element
    double nx = mesh1.boun_geoface(i,0);
    double ny = mesh1.boun_geoface(i,1);
    if(mesh1.bface(i,4) ==2)
    {
      //std::cout<<"case  wall cell"<<std::endl;
      for(int j=0;j<mesh1.ngauss_boun;j++) //loop over all gauss points of the boundary
      {
        matrix2d Ul = DG::U_at_poin(mesh1,mesh1.boun_geoface(i,2*j+2),mesh1.boun_geoface(i,2*j+3),le-1);
				matrix2d Ur(4,1); //init the right storage, this belongs to a ghost cell
				Ur(0,0) = Ul(0,0); //density
        Ur(3,0) = Ul(3,0);//energy
        // compute velocity using mirror images
				Ur(1,0) = Ul(1,0) - 2*(Ul(1,0)*nx+ Ul(2,0)*ny)*nx;
        Ur(2,0) = Ul(2,0) - 2*(Ul(1,0)*nx + Ul(2,0)*ny)*ny;
        matrix2d flux = FVS::Van_leer(mesh1,Ul,Ur,nx,ny); //reference to interface flux
        print2Term(flux);
        for(int k=0;k<mesh1.neqns;k++)
        {
          mesh1.rhsel(le-1, k, 0)  = -0.5*flux(0,0)*mesh1.bounweight*mesh1.boun_geoface(i,6) + mesh1.rhsel(le-1,k, 0);
          mesh1.rhsel(le-1, k, 1)  = -0.5*flux(0,0)*mesh1.bounweight*mesh1.boun_geoface(i,6)*(mesh1.int_geoface(i,2*j+2)- mesh1.geoel(le-1,1))/mesh1.geoel(le-1,3) + mesh1.rhsel(le-1,k,1);
          mesh1.rhsel(le-1, k, 2)  = -0.5*flux(0,0)*mesh1.bounweight*mesh1.boun_geoface(i,6)*(mesh1.int_geoface(i,2*j+3)- mesh1.geoel(le-1,2))/mesh1.geoel(le-1,4) + mesh1.rhsel(le-1,k,2);
        }
      }
    }
    //far field flag is 4 
    else
    {
      for(int j=0;j<mesh1.ngauss_boun;j++)
      {
        matrix2d Ul = DG::U_at_poin(mesh1,mesh1.boun_geoface(i,2*j+2),mesh1.boun_geoface(i,2*j+3),le-1);
        matrix2d Ur(4,1);
				Ur = mesh1.U_infty;
        matrix2d flux = FVS::Van_leer(mesh1,Ul,Ur,nx,ny); //reference to interface flux
        for(int k =0;k<mesh1.neqns;k++)
        {
          mesh1.rhsel(le-1, k, 0)  = -0.5*flux(0,0)*mesh1.bounweight*mesh1.boun_geoface(i,6) + mesh1.rhsel(le-1,k, 0);
          mesh1.rhsel(le-1, k, 1)  = -0.5*flux(0,0)*mesh1.bounweight*mesh1.boun_geoface(i,6)*(mesh1.int_geoface(i,2*j+2)- mesh1.geoel(le-1,1))/mesh1.geoel(le-1,3) + mesh1.rhsel(le-1,k,1);
          mesh1.rhsel(le-1, k, 2)  = -0.5*flux(0,0)*mesh1.bounweight*mesh1.boun_geoface(i,6)*(mesh1.int_geoface(i,2*j+3)- mesh1.geoel(le-1,2))/mesh1.geoel(le-1,4) + mesh1.rhsel(le-1,k,2);
        }
      }
    }
  }
}

// endsub


//  sub to compute the contribution of the boundary integral at interfaces to RHSel from internal faces and compute local timestep for the next iteration
void DG::rhsboun_iface(grid::mesh &mesh1)
{
  //std::cout<<"iface called"<<std::endl;
  feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT); 
  for (int i = 0; i < mesh1.nintface; i++) // loop over all the internal faces
  {
    double nx = mesh1.int_geoface(i, 0); // xcomponent of area normal vector, needed for the flux
    double ny = mesh1.int_geoface(i, 1); // y component of area normal vector
    int &le = mesh1.intface(i, 2);        // element to the left of the face when going from p1 to p2
    int &re = mesh1.intface(i, 3);        // element to the right of the face when going from p1 to p2
    for (int j = 0; j < mesh1.ngauss_boun; j++)
    {
      matrix2d Ul = DG::U_at_poin(mesh1, mesh1.int_geoface(i, 2*j + 2), mesh1.int_geoface(i, 2*j + 3), le-1);
      matrix2d Ur = DG::U_at_poin(mesh1, mesh1.int_geoface(i, 2*j + 2), mesh1.int_geoface(i, 2*j + 3), re-1);
      matrix2d flux = FVS::Van_leer(mesh1,Ul,Ur,nx,ny);// refrence to flux vector obtained from Roes reimann flux solver
      for(int k=0;k<mesh1.neqns;k++)
      {
        mesh1.rhsel(le-1, k, 0)  = 0.5*flux(0,0)*mesh1.bounweight*mesh1.int_geoface(i,6) + mesh1.rhsel(le-1,k, 0);
        mesh1.rhsel(le-1, k, 1)  = 0.5*flux(0,0)*mesh1.bounweight*mesh1.int_geoface(i,6)*(mesh1.int_geoface(i,2*j+2)- mesh1.geoel(le-1,1))/mesh1.geoel(le-1,3) + mesh1.rhsel(le-1,k,1);
        mesh1.rhsel(le-1, k, 2)  = 0.5*flux(0,0)*mesh1.bounweight*mesh1.int_geoface(i,6)*(mesh1.int_geoface(i,2*j+3)- mesh1.geoel(le-1,2))/mesh1.geoel(le-1,4) + mesh1.rhsel(le-1,k,2);
      
			//right hand side pushes flux contribution substracted due to the opposing direction of the outward normal
        mesh1.rhsel(re-1, k, 0)  =  0.5*flux(0,0)*mesh1.bounweight*mesh1.int_geoface(i,6) + mesh1.rhsel(re-1,k, 0);
        mesh1.rhsel(re-1, k, 1)  =  0.5*flux(0,0)*mesh1.bounweight*mesh1.int_geoface(i,6)*(mesh1.int_geoface(i,2*j+2)- mesh1.geoel(re-1,1))/mesh1.geoel(re-1,3) + mesh1.rhsel(re-1,k,1);
        mesh1.rhsel(re-1, k, 2)  =  0.5*flux(0,0)*mesh1.bounweight*mesh1.int_geoface(i,6)*(mesh1.int_geoface(i,2*j+3)- mesh1.geoel(re-1,2))/mesh1.geoel(re-1,4) + mesh1.rhsel(re-1,k,2);
      }
    }
  }
}
// endsub


matrix2d FVS::Van_leer(grid::mesh &mesh1, matrix2d &Ul, matrix2d &Ur, double &nx, double &ny)
{
  //std::cout<<"flux"<<std::endl;
  double gam = const_properties::gamma;
  matrix2d flux(4,1);
  matrix2d f_plus(4,1);
  matrix2d f_minus(4,1);
  double ci2 = gam*EOS::perf_gas(Ul)/Ul(0,0);
  if(ci2<0)
  {
    //std::cout<<"floatingpoint undef host"<<ci2<<std::endl;
  }
  double cj2 = gam*EOS::perf_gas(Ul)/Ul(0,0);
  if(cj2<0)
  {
    //std::cout<<"floatingpoint undef right or ghost "<<cj2<<std::endl;
  }
  double ci = DG::vel_sound(Ul);
  double cj = DG::vel_sound(Ur);
  double invRhoL = 1.0/Ul(0,0);
  double invRhoR = 1.0/Ur(0,0);
  double Vi_n = (Ul(1,0)*nx + Ul(2,0)*ny)*invRhoL;
  double Vj_n = (Ur(1,0)*nx + Ur(2,0)*ny)*invRhoR;
  double Mi_n = Vi_n/ci;
  double Mj_n = Vj_n/cj;
  //std::cout<<Mi_n<<" "<<Mj_n<<std::endl;
  //compute positive fluxes for left cell i 
  if(Mi_n>1)
  {
   f_plus(0,0) = Ul(0,0)*Vi_n;
   f_plus(1,0) = Ul(0,0)*Vi_n*Ul(1,0)*invRhoL + EOS::perf_gas(Ul)*nx;
   f_plus(2,0) = Ul(0,0)*Vi_n*Ul(2,0)*invRhoL + EOS::perf_gas(Ul)*ny;
   f_plus(3,0) = Vi_n*(Ul(3,0)+EOS::perf_gas(Ul));
  }
  else if(Mi_n<-1)
  {
    //f_plus.reset(); //make all elements zero
  }
  else
  {
    double u = Ul(1,0)*invRhoL;
    double v = Ul(2,0)*invRhoL;
    f_plus(0,0) = 0.25*Ul(0,0)*ci*std::pow((Mi_n+1),2);
    f_plus(1,0) = f_plus(0,0)*(Ul(1,0)*invRhoL+ nx*(-Vi_n+2*ci)/const_properties::gamma);
    f_plus(2,0) = f_plus(0,0)*(Ul(2,0)*invRhoL+ ny*(-Vi_n+2*ci)/const_properties::gamma);
    f_plus(3,0) = (std::pow(((gam-1)*Vi_n+2*ci),2)/(2*gam*gam - 1) + 0.5*(u*u + v*v - Vi_n*Vi_n))*f_plus(0,0);
  }
  // compute negative fluxes for right cell j
  if(Mj_n>1)
  {
    //f_minus.reset();
  }
  else if(Mj_n<-1)
  {

    f_minus(0,0) = Ur(0,0)*Vj_n;
    f_minus(1,0) = Ur(0,0)*Vj_n*Ur(1,0)*invRhoR + EOS::perf_gas(Ur)*nx;
    f_minus(2,0) = Ur(0,0)*Vj_n*Ur(2,0)*invRhoR + EOS::perf_gas(Ur)*ny;
    f_minus(3,0) = Vj_n*(Ur(3,0)+EOS::perf_gas(Ur));
    
  }
  else
  {
    double u = Ur(1,0)*invRhoR;
    double v = Ur(2,0)*invRhoR;
    f_minus(0,0) = -0.25*Ur(0,0)*cj*std::pow((Mj_n-1),2);
    f_minus(1,0) = f_minus(0,0)*(Ur(1,0)*invRhoR+ nx*(-Vj_n-2*cj)/const_properties::gamma);
    f_minus(2,0) = f_minus(0,0)*(Ur(2,0)*invRhoR+ ny*(-Vj_n-2*cj)/const_properties::gamma);
    f_minus(3,0) = (std::pow(((gam-1)*Vj_n-2*cj),2)/(2*gam*gam - 1) + 0.5*(u*u + v*v - Vj_n*Vj_n))*f_minus(0,0);
  }
  for(int m=0;m<mesh1.neqns;m++)
  {
    flux(m,0)= f_plus(m,0) + f_minus(m,0);
  }
 ///print2Term(flux);
  return flux;
}

matrix2d DG::fv_state(grid::mesh &mesh1, int i)
{
  //std::cout<<"normal overload func"<<std::endl;
  matrix2d U(4,1);
  for(int k=0;k<mesh1.neqns;k++)
  {  
    U(k,0) = mesh1.fvunkel(i,k);
  }
  return U;
}

// sub to get boundary face integral
void DG::fv_rhsboun_bface(grid::mesh &mesh1)
{
  //loop through all booundary faces 
  for(int i=0;i<mesh1.nbface;i++)
  {
    int &le = mesh1.bface(i,2); //get the host element
    matrix2d Ul = DG::fv_state(mesh1,le-1);//get the host conservative varsa
    double &nx = mesh1.boun_geoface(i,0);
    double &ny = mesh1.boun_geoface(i,1);
    if(mesh1.bface(i,4)==2) //this is a wall face
    {
      matrix2d Ur(4,1);
      Ur(3,0) = Ul(3,0); //eneergy
      Ur(0,0) = Ul(0,0); //density 
      Ur(1,0) = Ul(1,0)*nx -2*(Ul(1,0)*nx + Ul(2,0)*ny)*nx;
      Ur(2,0) = Ul(2,0)*ny -2*(Ul(1,0)*nx + Ul(2,0)*ny)*ny;
      matrix2d flux = FVS::Van_leer(mesh1,Ul, Ur, nx, ny);                            
      for(int m=0;m<mesh1.neqns;m++)
      {
        mesh1.fvrhsel(le-1,m) = mesh1.fvrhsel(le-1,m)-mesh1.boun_geoface(i,6)*flux(m,0);
      }
    }
    else if(mesh1.bface(i,4) ==4) //this is a far field cell
    {
      matrix2d Ur = mesh1.U_infty; //put the far field values for the ghost cell 
  
      matrix2d flux = FVS::Van_leer(mesh1, Ul, Ur, nx, ny);
      for(int m=0;m<mesh1.neqns;m++)
      {
        mesh1.fvrhsel(le-1,m) = mesh1.fvrhsel(le-1,m)-mesh1.boun_geoface(i,6)*flux(m,0);
      }
    }
  }
}


void DG::fv_rhsboun_iface(grid::mesh &mesh1)
{ 
  for(int i=0;i<mesh1.nintface;i++)
  {
    int &le = mesh1.intface(i,2);
    int &re = mesh1.intface(i,3);
    double &nx = mesh1.int_geoface(i,0);
    double &ny = mesh1.int_geoface(i,1);
    matrix2d Ul = DG::fv_state(mesh1, le-1);
    matrix2d Ur = DG::fv_state(mesh1, re-1);
    matrix2d flux = FVS::Van_leer(mesh1,Ul,Ur,nx,ny);
    for(int m=0;m<mesh1.neqns;m++)
    {
      mesh1.fvrhsel(le-1,m) = mesh1.fvrhsel(le-1,m) - mesh1.int_geoface(i,6)*flux(m,0);
      mesh1.fvrhsel(re-1,m) = mesh1.fvrhsel(re-1,m) + mesh1.int_geoface(i,6)*flux(m,0);
    }
  }
}

// is the value of the elem being looped, e is the value esuel-1 gives
matrix2d DG::fv_state(grid::mesh &mesh1, int i, int e, double &nx, double &ny)
{
  //std::cout<<"overload called for esuel"<<e<<std::endl;
  matrix2d U(4,1);
  matrix2d host = DG::fv_state(mesh1,i);
  if(e==-2)
  {
    U(0,0) = host(0,0);
    U(3,0) = host(3,0);
    U(1,0) = host(1,0)*nx - 2*(host(1,0)*nx + host(2,0)*ny)*nx;
    U(2,0) = host(2,0)*ny - 2*(host(1,0)*nx + host(2,0)*ny)*ny;
    //std::cout<<"wall cell overload complete"<<std::endl;
  }
  else if(e==-4)
  {
    U = mesh1.U_infty;
  }
  else
  {
    U = DG::fv_state(mesh1,e);
  }
  //std::cout<<"overload complete"<<std::endl;
  return U;   
}















