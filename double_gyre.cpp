//This code generates double gyre temporal velocity fields

#include <iostream>
#include "math.h"
#include <vector>
#include <string.h>
#include <fstream>
#include <sstream>
using namespace std;

//velocity function declaration
double u_velocity(double x1, double y1, double t); 
double v_velocity(double x1, double y1, double t); 

//Main function
int main()
{

	std::cout << "\033[1;32m This code generates vtk files for a double gyre flow \033[0m\n"; 
	double x_i=0.0, x_f=2.0, y_i=0.0, y_f=1.0; 
	vector<double> x, y, u, v; 
	double x_temp, y_temp, u_temp, v_temp; 
	double dx=0.01, dy=0.01; 
	int nx, ny;
	double t_i=0.0, t_f=10.0, dt=0.1, t_temp; 
	int i, j, k, m, nt; 
	ofstream vtkfile; 
	
	//Enter all the simulation values
	/**cout<<"enter the value of initial time"<<endl; 
	cin>>t_i;  
	cout<<"enter the value of final time"<<endl;
	cin>>t_f;
	cout<<"enter the time step size"<<endl;
	cin>>dt;
	
	nt = (t_f-t_i)/dt; 

	cout<<"the value of the number of time steps is"<<'\t'<<nt<<endl;

	cout<<"enter the value of dx and dy"<<endl;
	cin>>dx; 
	dy = dx; 

	cout<<"enter the ranges of x"<<endl;
	cin>>x_i>>x_f;

	cout<<"enter the ranges of y"<<endl;
	cin>>y_i>>y_f; 
	**/


	nt = (t_f-t_i)/dt;  			//total number of time steps
	nx = (x_f-x_i)/dx; 			//total number of cells along x-direction
	ny = (y_f-y_i)/dy; 			//total number of cells along y-direction
	

	//Generate the grid for the solution
	for(i=0; i<=nx; i++)
	{
		x_temp = x_i+dx*i; 		
		x.push_back(x_temp); 
	}

	for(j=0; j<=ny; j++)
	{
		y_temp = y_i+dy*j; 
		y.push_back(y_temp); 
	}


	//Outer time loop (Do everything inside this loop -> calculate velocities, and write vtk files)
	for(k=0; k<=nt; k++)
	{
		t_temp = k*dt; 

		//Calculate velocities
		for(j=0; j<=ny; j++)
		{
			for(i=0; i<=nx; i++)
			{
				u_temp = u_velocity(x[i],y[j],t_temp);
				u.push_back(u_temp); 
			}
		}
		
		for(j=0; j<=ny; j++)
		{
			for(i=0; i<=nx; i++)
			{
				v_temp = v_velocity(x[i],y[j],t_temp); 
				v.push_back(v_temp); 
			}
		}


		//WRITE DATA TO VTK FILES
		//Read strings to make variable names
		ostringstream fileNameStream;					//Create a string stream
   		fileNameStream <<"velocity_"<<k<<".vtk";			//Give the name of the string stream
 		std::string fileName = fileNameStream.str();			//Assign the string stream to some other string variable
  		vtkfile.open(fileName.c_str());					//open the output file stream with the string name

		vtkfile<<"# vtk DataFile Version 3.0\n"; 
		vtkfile<<"Header\n"; 
		vtkfile<<"ASCII\n"; 
		vtkfile<<"DATASET RECTILINEAR_GRID\n";	
		vtkfile<<"FIELD FieldData 1\n"; 
		vtkfile<<"TIME 1 1 double\n";
		vtkfile<<k*dt<<"\n";
		vtkfile<<"DIMENSIONS "<<nx+1<<" "<<ny+1<<" "<<"1\n"; 	
		vtkfile<<"X_COORDINATES "<<nx+1<< " double\n"; 
		for(i=0; i<=nx; i++)
		{
			vtkfile<<x[i]<<" ";
		}
		vtkfile<<"\n"; 
		vtkfile<<"Y_COORDINATES "<<ny+1<< " double\n"; 
		for(j=0; j<=ny; j++)
		{
			vtkfile<<y[j]<<" ";
		}
		vtkfile<<"\n"; 
		vtkfile<<"Z_COORDINATES 1 double\n"; 
		vtkfile<<0.0<<"\n";
		//write velocity values
		vtkfile<<"POINT_DATA "<<(nx+1)*(ny+1)<<"\n";
		vtkfile<<"VECTORS velocity double\n"; 
		m=0; 
		for(j=0; j<=ny; j++)
		{
			for(i=0; i<=nx; i++)
			{
				vtkfile<<u[m]<<"\t"<<v[m]<<"\t"<<"0.0 \n"; 
				m=m+1;
			}
		}
	

		vtkfile.close(); 						//close the file stream


		//Clear the vectors so that it will be free for the next outer loop iteration
		u.clear(); 
		v.clear(); 

	}				
	
	//clear the vectors of x and y
	x.clear(); 
	y.clear(); 
			
	return 0; 
}



//function definition
double u_velocity(double x1, double y1, double t)
{
	double u_vel; 
	
	//u_vel=-3.141592*1*sin(3.141592*x1)*cos(3.141592*y1); 			//steady state double_gyre velocity

	u_vel = -3.141592*0.25*sin(3.141592*((0.25*sin(2*3.141592*t)*x1*x1)+(1-(0.5*sin(2*3.141592*t)))*x1))*cos(3.141592*y1);

	return u_vel; 
}

double v_velocity(double x1, double y1, double t)
{
	double v_vel; 

	v_vel= 3.141592*0.25*cos(3.141592*((0.25*sin(2*3.141592*t)*x1*x1)+(1-(0.5*sin(2*3.141592*t)))*x1))*sin(3.141592*y1)*(2*0.25*sin(2*3.141592*t)*x1+(1-(0.5*sin(2*3.141592*t)))); 

	return v_vel; 
}



