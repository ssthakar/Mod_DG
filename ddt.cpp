#include "dg.h"
#include "mesh.h"
#include "vmatrix2.h"
#include <algorithm>

// function to calculate local time step for cell i
double ddt::local_ts(grid::mesh &mesh1, int &i)
{
  double denom = 0.0;
  //std::cout<<"local timestep function for element called: "<<i<<std::endl;
	double &area = mesh1.geoel(i,0); //area of the triangular element
	const int &ip1 = mesh1.inpoel(i, 0); // get the points that make up the cell 
	const int &ip2 = mesh1.inpoel(i, 1);
	const int &ip3 = mesh1.inpoel(i,2);
	double &p1x = mesh1.coords(ip1 - 1, 0);
	double &p2x = mesh1.coords(ip2 - 1, 0);
	double &p3x = mesh1.coords(ip3 - 1, 0);
	double &p1y = mesh1.coords(ip1 - 1, 1);
	double &p2y = mesh1.coords(ip2 - 1, 1);
	double &p3y = mesh1.coords(ip3 - 1, 1);
	double mag1 = sqrt(std::pow((p2y-p1y),2) + std::pow((p2x-p1x),2));// push the components of Area weighted normal vectors to the geoface matrix
	double mag2 = sqrt(std::pow((p3y-p2y),2) + std::pow((p3x-p2x),2));
	double mag3 = sqrt(std::pow((p3y-p1y),2) + std::pow((p3x-p1x),2));
	double nx1 = (p2y - p1y)/mag1; //unit normal vector components x and y 
	double ny1 = -1*(p2x-p1x)/mag1;
	double nx2 = (p3y - p2y)/mag2; //unit normal vector components x and y 
	double ny2 = -1*(p3x-p2x)/mag2;
	double nx3 = (p1y - p3y)/mag3; //unit normal vector components x and y 
	double ny3 = -1*(p1x-p3x)/mag3;
	double len1 = len(p1x,p2x,p1y,p2y);
	double len2 = len(p2x,p3x,p2y,p3y);
	double len3 = len(p1x,p3x,p1y,p3y);
  //std::cout<<nx1*nx1 + ny1*ny1<<std::endl;
  for(int j=0;j<mesh1.ntype;j++)  //loop over all the edges of the cell 
	{
		switch(j)
		{
			case 0:
			{
				matrix2d Ui = DG::fv_state(mesh1,i);
        matrix2d Uj = DG::fv_state(mesh1,i,mesh1.esuel(i,j) - 1,nx1,ny1);
			  //std::cout<<"break point 68"<<EOS::perf_gas(Uj)<<std::endl;
        //std::cout<<"left cell"<<std::endl;
        double ci = DG::vel_sound(Ui);
				
        //std::cout<<"left cell"<<std::endl;
        double cj = DG::vel_sound(Uj);
        
        denom = denom + len1*0.5*(ci+cj + std::fabs(Ui(1,0)/Ui(0,0)*nx1 + Uj(1,0)/Uj(0,0)*nx1 + Ui(2,0)/Ui(0,0)*ny1 + Uj(2,0)/Uj(0,0)*ny1));
				break;
			}
			case 1:
			{
  			matrix2d Ui = DG::fv_state(mesh1,i);
        //std::cout<<"left cell"<<std::endl;
        matrix2d Uj = DG::fv_state(mesh1,i,mesh1.esuel(i,j)-1,nx2,ny2);      //std::cout<<" case 1"<<std::endl;
        //std::cout<<"right cell"<<std::endl;
				double ci = DG::vel_sound(Ui);
		 	  //std::cout<<"break point 69"<<std::endl;
				double cj = DG::vel_sound(Uj);
				denom = denom + len2*0.5*(ci+cj + std::fabs(Ui(1,0)/Ui(0,0)*nx2 + Uj(1,0)/Uj(0,0)*nx2 + Ui(2,0)/Ui(0,0)*ny2 + Uj(2,0)/Uj(0,0)*ny2));
				break;
			}
			case 2:
			{
 				matrix2d Ui = DG::fv_state(mesh1,i);
        //std::cout<<"left cell"<<std::endl;
        matrix2d Uj = DG::fv_state(mesh1,i,mesh1.esuel(i,j) - 1,nx3,ny3);       //std::cout<<" case 2"<<std::endl;
        //std::cout<<"right cell"<<std::endl;
				double ci = DG::vel_sound(Ui);
				double cj = DG::vel_sound(Uj);
				denom = denom + 0.5*len3*(ci+cj + std::fabs(Ui(1,0)/Ui(0,0)*nx3 + Uj(1,0)/Uj(0,0)*nx3 + Ui(2,0)/Ui(0,0)*ny3 + Uj(2,0)/Uj(0,0)*ny3));
				break;
			}
		}
	}
  //std::cout<<"timestep: "<<i<<" "<<const_properties::CFL*area/denom<<std::endl;
  feenableexcept(FE_INVALID | FE_OVERFLOW);
  //std::cout<<"local timestep function complete"<<std::endl;
	return area/denom;
}

//subroutine for TVD RK3 takes in reference to current mesh object
void ddt::RK3::RK3_outer(grid::mesh &mesh1,soln &soln1)
{
  
	int iter_count = 0;
	while(iter_count<soln1.max_iter)
	{
		mesh1.RK_stor = mesh1.unkel; //store the current solution for further use
    ddt::RK3::RK_s1(mesh1);
    std::cout<<"timestep global is: "<<const_properties::CFL*mesh1.dt<<std::endl;
    ddt::RK3::RK_s2(mesh1);
		ddt::RK3::RK_s3(mesh1);
		DG::residual(mesh1); //compute the residual;
    std::cout<<"Residual: \n "<<mesh1.res_vec(0,0)<<"\n"<<mesh1.res_vec(1,0)<<"\n"<<mesh1.res_vec(2,0)<<"\n"<<mesh1.res_vec(3,0)<<std::endl;
		//calculate the next time step
	  //DG::delta_T(mesh1); //update local time steps for every element 
    //check for convergence and break if converged
		//if(DG::isSolnConverged(mesh1,soln1) == 1)
			//break;
    std::cout<<"Iteration: "<<iter_count<<std::endl;
		iter_count ++;
    
	}
  
}

void ddt::RK3::RK_s1(grid::mesh &mesh1) //this also serves forward euler
{
  //mesh1.fvRkstor = mesh1.fvunkel;
  DG::delta_T(mesh1);
  //std::cout<<"stage 1 started"<<std::endl;
	//DG::delta_T(mesh1);
  DG::fv_rhsboun_iface(mesh1);
	DG::fv_rhsboun_bface(mesh1);
  for(int i = 0;i<mesh1.nelem;i++) //loop through all elements
	{
		double &M1 = mesh1.geoel(i,11);
		double &M2 = mesh1.geoel(i,12);
		double &M3 = mesh1.geoel(i,13);
		double dT = mesh1.dt;
    //double &dT = mesh1.geoel(i,14);
		for(int m=0;m<mesh1.neqns;m++) //loop through all conservative variables
		{
			mesh1.fvunkel(i,m) = mesh1.fvunkel(i,m)  + dT*mesh1.fvrhsel(i,m)/mesh1.geoel(i,0);
			//mesh1.unkel(i,m,1) = mesh1.unkel(i,m,1) + dT*(M3*mesh1.rhsel(i,m,1) - M2*mesh1.rhsel(i,m,2))/(M1*M3 - M2*M2);
			//mesh1.unkel(i,m,2) = mesh1.unkel(i,m,2) + dT*(M1*mesh1.rhsel(i,m,2) - M2*mesh1.rhsel(i,m,1))/(M1*M3 - M2*M2);
		}
		//std::cout<<mesh1.fvunkel(i,1)<<std::endl; 
	}
  //std::cout<<"stage 1 complete: "<<std::endl;
}


void ddt::RK3::RK_s2(grid::mesh &mesh1)
{
  //std::cout<<"stage 2 started"<<std::endl;
	//update RHS to get the contribution of the latest solution state 
  DG::delta_T(mesh1);
	DG::fv_rhsboun_bface(mesh1);
	DG::fv_rhsboun_iface(mesh1);
	for(int i = 0;i<mesh1.nelem;i++) //loop through all elements
	{
		double &M1 = mesh1.geoel(i,11);
		double &M2 = mesh1.geoel(i,12);
		double &M3 = mesh1.geoel(i,13);
		//double &dT = mesh1.geoel(i,14);
		double dT = mesh1.dt;
    for(int m=0;m<mesh1.neqns;m++) //loop through all conservative variables
		{
			mesh1.fvunkel(i,m) = 3*(mesh1.fvRkstor(i,m))/4 + 0.25*(mesh1.fvunkel(i,m) + dT*mesh1.fvrhsel(i,m)/mesh1.geoel(i,0));
			//mesh1.unkel(i,m,1) = 3*(mesh1.RK_stor(i,m,0))/4 + 0.25*(mesh1.unkel(i,m,1) + dT*(M3*mesh1.rhsel(i,m,1) - M2*mesh1.rhsel(i,m,2))/(M1*M3 - M2*M2));
			//mesh1.unkel(i,m,2) = 3*(mesh1.RK_stor(i,m,0))/4 + 0.25*(mesh1.unkel(i,m,2) + dT*(M1*mesh1.rhsel(i,m,2) - M2*mesh1.rhsel(i,m,1))/(M1*M3 - M2*M2));
		}
	}
  //std::cout<<"stage 2 complete"<<std::endl;
}

void ddt::RK3::RK_s3(grid::mesh &mesh1)
{
  //std::cout<<"stage 3 started"<<std::endl;
	//update RHS to get the contribution of the latest solution state 
  DG::delta_T(mesh1);
	DG::fv_rhsboun_bface(mesh1);
	DG::fv_rhsboun_iface(mesh1);
	for(int i = 0;i<mesh1.nelem;i++) //loop through all elements
	{
		double &M1 = mesh1.geoel(i,11);
		double &M2 = mesh1.geoel(i,12);
		double &M3 = mesh1.geoel(i,13);
		//double &dT = mesh1.geoel(i,14);
		double dT = mesh1.dt;
    for(int m=0;m<mesh1.neqns;m++) //loop through all conservative variables
		{
			mesh1.fvunkel(i,m) = 1*(mesh1.fvRkstor(i,m))/3 + 2*(mesh1.fvunkel(i,m) + dT*mesh1.fvrhsel(i,m)/mesh1.geoel(i,0))/3;
		//	mesh1.unkel(i,m,1) = 1*(mesh1.RK_stor(i,m,0))/3 + 2*(mesh1.unkel(i,m,1) + dT*(M3*mesh1.rhsel(i,m,1) - M2*mesh1.rhsel(i,m,2))/(M1*M3 - M2*M2))/3;
			//mesh1.unkel(i,m,2) = 1*(mesh1.RK_stor(i,m,0))/3 + 2*(mesh1.unkel(i,m,2) + dT*(M1*mesh1.rhsel(i,m,2) - M2*mesh1.rhsel(i,m,1))/(M1*M3 - M2*M2))/3;
		}
	}

  //std::cout<<"stage 3 complete"<<std::endl;
}

void DG::delta_T(grid::mesh &mesh1)
{
	for(int i = 0;i<mesh1.nelem;i++)
	{
		mesh1.Ltimestep[i] = ddt::local_ts(mesh1,i);
    //std::cout<<"element: "<<i<<std::endl;  
    feenableexcept(FE_INVALID);
  }
  mesh1.dt = *std::min_element(mesh1.Ltimestep.begin(),mesh1.Ltimestep.end())*const_properties::CFL;
}

void DG::residual(grid::mesh &mesh1)
{
	mesh1.res_vec.reset(); //set all values to 0	
	for(int i=0;i<mesh1.nelem;i++)
	{
		for(int m=0;m<mesh1.neqns;m++)
		{
			mesh1.res_vec(m,0) = mesh1.res_vec(m,0) + sqrt(mesh1.geoel(i,0)*std::pow((std::fabs(mesh1.fvunkel(i,m) - mesh1.fvRkstor(i,m))),2));
		}
	}
}

bool DG::isSolnConverged(grid::mesh &mesh1, soln &soln1)
{
	bool result = true;
	for(int i=0;i<mesh1.neqns;i++)
	{
		if(mesh1.res_vec(i,0)>soln1.abstol)
			result = false;
	}
	return result;
}




















