#include "dg.h"
#include "mesh.h"

// function to calculate local time step for cell i
double ddt::local_ts(grid::mesh &mesh1, int &i)
{
	double denom = 0.0; //declare the place holder that will store in the boundary sum
	double &area = mesh1.geoel(i,0); //area of the triangular element
	int &ip1 = mesh1.inpoel(i, 0); // get the points that make up the cell 
	int &ip2 = mesh1.inpoel(i, 1);
	int &ip3 = mesh1.inpoel(i,2);
	double &p1x = mesh1.coords(ip1 - 1, 0);
	double &p2x = mesh1.coords(ip2 - 1, 0);
	double &p3x = mesh1.coords(ip3 - 1, 0);
	double &p1y = mesh1.coords(ip1 - 1, 1);
	double &p2y = mesh1.coords(ip2 - 1, 1);
	double &p3y = mesh1.coords(ip3 - 1, 1);
	double mag1 = sqrt(std::pow((p2y-p1y),2) + std::pow((p2x-p1x),2));// push the components of Area weighted normal vectors to the geoface matrix
	double mag2 = sqrt(std::pow((p2y-p1y),2) + std::pow((p2x-p1x),2));
	double mag3 = sqrt(std::pow((p2y-p1y),2) + std::pow((p2x-p1x),2));
	double nx1 = (p2y - p1y)/mag1; //unit normal vector components x and y 
	double ny1 = -1*(p2x-p1x)/mag1;
	double nx2 = (p3y - p2y)/mag2; //unit normal vector components x and y 
	double ny2 = -1*(p3x-p2x)/mag2;
	double nx3 = (p1y - p3y)/mag3; //unit normal vector components x and y 
	double ny3 = -1*(p1x-p3x)/mag3;
	double len1 = len(p1x,p2x,p1y,p2y);
	double len2 = len(p2x,p3x,p2y,p3y);
	double len3 = len(p1x,p3x,p1y,p3y);

	for(int j=0;j<mesh1.ntype;j++)  //loop over all the edges of the cell 
	{
		switch(j)
		{
			case 0:
			{// side with p1 and p2 as end points
				matrix2d Ui = DG::U_at_poin(mesh1, mesh1.geoel(i, 2*j+5), mesh1.geoel(i, 2*j+6), i);
				matrix2d Uj = DG::U_at_poin(mesh1, mesh1.geoel(i, 2*j+5), mesh1.geoel(i, 2*j+6), mesh1.esuel(i,j));
				double ci = DG::vel_sound(Ui);
				double cj = DG::vel_sound(Uj);
				denom = denom + len1*0.5*(ci+cj + std::fabs(Ui(1,0)/Ui(0,0)*nx1 + Uj(1,0)/Uj(0,0)*nx1 + Ui(2,0)/Ui(0,0)*ny1 + Uj(2,0)/Uj(0,0)*ny1));
				break;
			}
			case 1:
			{
				matrix2d Ui = DG::U_at_poin(mesh1, mesh1.geoel(i, 2*j+5), mesh1.geoel(i, 2*j+6), i);
				matrix2d Uj = DG::U_at_poin(mesh1, mesh1.geoel(i, 2*j+5), mesh1.geoel(i, 2*j+6), mesh1.esuel(i,j));
				double ci = DG::vel_sound(Ui);
				double cj = DG::vel_sound(Uj);
				denom = denom + len2*0.5*(ci+cj + std::fabs(Ui(1,0)/Ui(0,0)*nx2 + Uj(1,0)/Uj(0,0)*nx2 + Ui(2,0)/Ui(0,0)*ny2 + Uj(2,0)/Uj(0,0)*ny2));
				break;
			}
			case 2:
			{
				matrix2d Ui = DG::U_at_poin(mesh1, mesh1.geoel(i, 2*j+5), mesh1.geoel(i, 2*j+6), i);
				matrix2d Uj = DG::U_at_poin(mesh1, mesh1.geoel(i, 2*j+5), mesh1.geoel(i, 2*j+6), mesh1.esuel(i,j));
				double ci = DG::vel_sound(Ui);
				double cj = DG::vel_sound(Uj);
				denom = denom + 0.5*len3*(ci+cj + std::fabs(Ui(1,0)/Ui(0,0)*nx3 + Uj(1,0)/Uj(0,0)*nx3 + Ui(2,0)/Ui(0,0)*ny3 + Uj(2,0)/Uj(0,0)*ny3));
				break;
			}
		}
	}
	return area/denom;
}

//subroutine for TVD RK3 takes in reference to current mesh object
void ddt::RK3::RK3_outer(grid::mesh &mesh1,soln &soln1)
{
	int iter_count = 0;
	while(iter_count<soln1.max_iter)
	{
		mesh1.RK_stor = mesh1.unkel; //store the current solution for further use
		//go through all the stages of the multi-stage method, update solution at final stage 
		ddt::RK3::RK_s1(mesh1);
		ddt::RK3::RK_s2(mesh1);
		ddt::RK3::RK_s3(mesh1);
		iter_count++;
		//check for convergence and break if converged
		if(grid::run::isSolnConverged(mesh1) == 1)
			break;
	}
}




void ddt::RK3::RK_s1(grid::mesh &mesh1)
{
	//update RHS to get the contribution of the latest solution state 
	DG::rhsboun_bface(mesh1);
	DG::rhsboun_iface(mesh1);
	DG::rhsboun_iface(mesh1);
	grid::run::delta_T(mesh1); //generate local time steps for every element 
	//step in time and update solution storage, keep RK storage the same
	for(int i = 0;i<mesh1.nelem;i++) //loop through all elements
	{
		double &M1 = mesh1.geoel(i,11);
		double &M2 = mesh1.geoel(i,12);
		double &M3 = mesh1.geoel(i,13);
		double &dT = mesh1.geoel(i,15);
		for(int m=0;m<mesh1.neqns;m++) //loop through all conservative variables
		{
			mesh1.unkel(i,m,0) = mesh1.unkel(i,m,0) + mesh1.rhsel(i,m,0);
			mesh1.unkel(i,m,1) = mesh1.unkel(i,m,1) + dT*(M3*mesh1.rhsel(i,m,1) - M2*mesh1.rhsel(i,m,2))/(M1*M3 - M2*M2);
			mesh1.unkel(i,m,2) = mesh1.unkel(i,m,2) + dT*(M1*mesh1.rhsel(i,m,2) - M2*mesh1.rhsel(i,m,1))/(M1*M3 - M2*M2);
		}
	}
}


void ddt::RK3::RK_s2(grid::mesh &mesh1)
{
	//update RHS to get the contribution of the latest solution state 
	DG::rhsboun_bface(mesh1);
	DG::rhsboun_iface(mesh1);
	DG::rhsboun_iface(mesh1);
	grid::run::delta_T(mesh1); //update local time steps for every element 
	for(int i = 0;i<mesh1.nelem;i++) //loop through all elements
	{
		double &M1 = mesh1.geoel(i,11);
		double &M2 = mesh1.geoel(i,12);
		double &M3 = mesh1.geoel(i,13);
		double &dT = mesh1.geoel(i,15);
		for(int m=0;m<mesh1.neqns;m++) //loop through all conservative variables
		{
			mesh1.unkel(i,m,0) = 3*(mesh1.RK_stor(i,m,0))/4 + 0.25*(mesh1.unkel(i,m,0) + mesh1.rhsel(i,m,0));
			mesh1.unkel(i,m,1) = 3*(mesh1.RK_stor(i,m,0))/4 + 0.25*(mesh1.unkel(i,m,1) + dT*(M3*mesh1.rhsel(i,m,1) - M2*mesh1.rhsel(i,m,2))/(M1*M3 - M2*M2));
			mesh1.unkel(i,m,2) = 3*(mesh1.RK_stor(i,m,0))/4 + 0.25*(mesh1.unkel(i,m,2) + dT*(M1*mesh1.rhsel(i,m,2) - M2*mesh1.rhsel(i,m,1))/(M1*M3 - M2*M2));
		}
	}
}

void ddt::RK3::RK_s3(grid::mesh &mesh1)
{
	//update RHS to get the contribution of the latest solution state 
	DG::rhsboun_bface(mesh1);
	DG::rhsboun_iface(mesh1);
	DG::rhsboun_iface(mesh1);
	grid::run::delta_T(mesh1); //update local time steps for every element 
	for(int i = 0;i<mesh1.nelem;i++) //loop through all elements
	{
		double &M1 = mesh1.geoel(i,11);
		double &M2 = mesh1.geoel(i,12);
		double &M3 = mesh1.geoel(i,13);
		double &dT = mesh1.geoel(i,15);
		for(int m=0;m<mesh1.neqns;m++) //loop through all conservative variables
		{
			mesh1.unkel(i,m,0) = 1*(mesh1.RK_stor(i,m,0))/3 + 2*(mesh1.unkel(i,m,0) + mesh1.rhsel(i,m,0))/3;
			mesh1.unkel(i,m,1) = 1*(mesh1.RK_stor(i,m,0))/3 + 2*(mesh1.unkel(i,m,1) + dT*(M3*mesh1.rhsel(i,m,1) - M2*mesh1.rhsel(i,m,2))/(M1*M3 - M2*M2))/3;
			mesh1.unkel(i,m,2) = 1*(mesh1.RK_stor(i,m,0))/3 + 2*(mesh1.unkel(i,m,2) + dT*(M1*mesh1.rhsel(i,m,2) - M2*mesh1.rhsel(i,m,1))/(M1*M3 - M2*M2))/3;
		}
	}

}

void grid::run::delta_T(grid::mesh &mesh1)
{
	for(int i = 0;i<mesh1.nelem;i++)
	{
		mesh1.geoel(i,14) = ddt::local_ts(mesh1,i);
	}
}


