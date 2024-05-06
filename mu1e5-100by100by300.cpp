#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/numeric/odeint.hpp> // odeint function definitions
#include <math.h>
#include <cmath>
#include <omp.h>
#include <boost/math/special_functions/cos_pi.hpp>
#include <boost/math/special_functions/sin_pi.hpp>
#include <stdexcept>
#include <boost/throw_exception.hpp>
#include <boost/numeric/odeint/stepper/stepper_categories.hpp>
#include <boost/numeric/odeint/stepper/controlled_step_result.hpp>
#include <boost/numeric/odeint/integrate/max_step_checker.hpp>
#include <boost/numeric/odeint/integrate/detail/integrate_const.hpp>
#include <boost/numeric/odeint/util/bind.hpp>
#include <boost/numeric/odeint/util/unwrap_reference.hpp>
#include <boost/numeric/odeint/util/copy.hpp>
#include <boost/numeric/odeint/util/detail/less_with_sign.hpp>
#include <boost/numeric/odeint/integrate/null_observer.hpp>
#include <iostream>
#include <omp.h>
using namespace std;
using namespace boost::numeric::odeint;
typedef std::vector< double> state_type;
ofstream fp;
#define costhslices 100
#define xslices 400
#define yslices 400
int width=xslices/20;
double theta=1.0e-6;
double omega=1.0e-1;
double lam=0.0e0;
double Lx=20.0e0; // in km
double Ly=20.0e0; // in km
double c=3.0e5; // speed of light in km per sec
double mup=1.0e2;
int tsize = 8*costhslices*xslices*yslices;
double tharray[costhslices];
double dtharray[costhslices];
double d[8*xslices*yslices];
double dc[8*xslices*yslices];
double ds[8*xslices*yslices];
double th1=M_PI/3.0e0;
double th2=M_PI*2.0e0/3.0e0;
//#define energyslices 1
void my_observer( const state_type &x, const double t );
double normal(double x, double sig,double mean)
{
    double ans;
    //ans=1.0e0/(sqrt(2.0e0*M_PI*sig*sig));
    ans=1.0e0*exp(-(x-mean)*(x-mean)/(2.0e0*sig*sig));
    return(ans);
}
void initarrays(state_type &x)
{
    for(int k=0;k<costhslices;k++)
    {
        tharray[k]=(0.5e0+(double)k)*2.0e0*M_PI/(double)costhslices;
        dtharray[k]=2.0e0*M_PI/(double)costhslices;
    }
    for(int i=0;i<xslices;i++)
    {
	    for(int j=0;j<yslices;j++)
	    {
	        for(int k=0;k<costhslices;k++)
	        {
                
		    //double xmean=1.0e0/2.0e0;
		    //double ymean=1.0e0/2.0e0;
		    if(tharray[k] > th1 && tharray[k] < th2)
		    //if(k==20)
		    {
			//if(i==xslices/2 && j==yslices/2)
			//{
                    		x[i*yslices*costhslices*8+j*costhslices*8+k*8+0]=normal(i,width,(double)(xslices/2))*normal(j,width,(double)(yslices/2))*1.0e0/(th2-th1);
			//}
		    }
                    x[i*yslices*costhslices*8+j*costhslices*8+k*8+1]=0.0e0;
                    x[i*yslices*costhslices*8+j*costhslices*8+k*8+2]=0.0e0;
                    x[i*yslices*costhslices*8+j*costhslices*8+k*8+3]=0.0e0;
		    if(tharray[k] > 0.0e0 && tharray[k] < 2.0e0*M_PI)
		    //if(k==20)
                    { 
			//if(i==xslices/2 && j==yslices/2)
                        //{
                        	x[i*yslices*costhslices*8+j*costhslices*8+k*8+4]=normal(i,width,(double)(xslices/2))*normal(j,width,(double)(yslices/2))*0.5e0/(th2-th1);
			//}
		    }
                    x[i*yslices*costhslices*8+j*costhslices*8+k*8+5]=0.0e0;
                    x[i*yslices*costhslices*8+j*costhslices*8+k*8+6]=0.0e0;
                    x[i*yslices*costhslices*8+j*costhslices*8+k*8+7]=0.0e0;
                
	        }
	    }
    }
    //my_observer(x,1.0);
}
void transport(state_type &x, const double dt)
{
    double delip, deljp;
    int i1, i2;
    int j1, j2;
    int temp;
    int size=tsize;
    state_type x1(size);
    //cout<<"dt = "<<dt<<endl;
    for(int k=0;k<costhslices;k++)
    {
        delip=-(double)xslices/Lx*c*cos(tharray[k])*dt;
        deljp=-(double)yslices/Ly*c*sin(tharray[k])*dt;
        for(int i=0;i<xslices;i++)
    	{
            for(int j=0;j<yslices;j++)
            {
                if((i+delip)<0.0e0)
                {
                    i1=(int)(i+delip)-1+xslices;
                }
                else
                {
                    i1=(int)(i+delip);
                }
                if(i1>xslices-1)
                {
                    i1=i1-xslices;
                }
                i2=i1+1;
                if(i2>xslices-1)
                {
                    i2=i2-xslices;
                }
                if(delip<0.0e0)
                {
                    temp=i2;
                    i2=i1;
                    i1=temp;
                }
                if((j+deljp)<0.0e0)
                {
                    j1=(int)(j+deljp)-1+yslices;
                }
                else
                {
                    j1=(int)(j+deljp);
                }
                if(j1>xslices-1)
                {
                    j1=j1-yslices;
                }
                j2=j1+1;
                if(j2>yslices-1)
                {
                    j2=j2-yslices;
                }
                if(deljp<0.0e0)
                {
                    temp=j2;
                    j2=j1;
                    j1=temp;
                }
		
                //-----------------indices set
                int shiftx=i-i1;
                int shifty=j-j1;
                if((i-i1)<0)
                {
                    shiftx=shiftx+xslices;
                }
                if((j-j1)<0)
                {
                    shifty=shifty+yslices;
                }
                for(int l=0;l<8;l++)
                {
                x1[i*yslices*costhslices*8+j*costhslices*8+k*8+l]=x1[i*yslices*costhslices*8+j*costhslices*8+k*8+l]+x[i1*yslices*costhslices*8+j1*costhslices*8+k*8+l]*(1.0e0-abs(delip+shiftx))*(1.0e0-abs(deljp+shifty));
                x1[i*yslices*costhslices*8+j*costhslices*8+k*8+l]=x1[i*yslices*costhslices*8+j*costhslices*8+k*8+l]+x[i1*yslices*costhslices*8+j2*costhslices*8+k*8+l]*(1.0e0-abs(delip+shiftx))*(abs(deljp+shifty));
                x1[i*yslices*costhslices*8+j*costhslices*8+k*8+l]=x1[i*yslices*costhslices*8+j*costhslices*8+k*8+l]+x[i2*yslices*costhslices*8+j1*costhslices*8+k*8+l]*(abs(delip+shiftx))*(1.0e0-abs(deljp+shifty));
                x1[i*yslices*costhslices*8+j*costhslices*8+k*8+l]=x1[i*yslices*costhslices*8+j*costhslices*8+k*8+l]+x[i2*yslices*costhslices*8+j2*costhslices*8+k*8+l]*(abs(delip+shiftx))*(abs(deljp+shifty));
                }
		//cout<<"k="<<k<<" "<<"j="<<j<<" "<<"j1="<<j1<<" "<<"j2="<<j2<<" "<<"deljp="<<deljp<<endl;
            }
	    //cout<<"k="<<k<<" "<<"i="<<i<<" "<<"i1="<<i1<<" "<<"i2="<<i2<<" "<<"delip="<<delip<<endl;
        }
    }
    x=x1;
}
void my_system( const state_type &x , state_type &dxdt , const double t )
{
    //std::cout<<"t="<<t<<endl;
    for(int i=0;i<tsize;i++)
    {
	    dxdt[i]=0.0e0;
    }
    for(int i=0;i<xslices;i++)
    {
        for(int j=0;j<yslices;j++)
        {
            for(int l=0;l<8;l++)
            {
                d[i*yslices*8+j*8+l]=0.0e0;
                ds[i*yslices*8+j*8+l]=0.0e0;
                dc[i*yslices*8+j*8+l]=0.0e0;
            }
        }
    }
int th_id;
int nthreads;
#pragma omp parallel private(th_id)
{
    nthreads=omp_get_num_threads();
    th_id = omp_get_thread_num();
    //for(int i=0;i<xslices;i++)
    for(int ip=0;(ip+th_id)<xslices;ip=ip+nthreads)
    {
	int i=ip+th_id;
        for(int j=0;j<yslices;j++)
        {
            for(int k=0;k<costhslices;k++)
            {
		if(x[i*yslices*costhslices*8+j*costhslices*8+k*8+0]+x[i*yslices*costhslices*8+j*costhslices*8+k*8+1]+x[i*yslices*costhslices*8+j*costhslices*8+k*8+4]+x[i*yslices*costhslices*8+j*costhslices*8+k*8+5] > 1.0e-12)
                {
                    for(int l=0;l<8;l++)
                    {
                    	d[i*yslices*8+j*8+l]=d[i*yslices*8+j*8+l]+mup*x[i*yslices*costhslices*8+j*costhslices*8+k*8+l]*dtharray[k];
                    	ds[i*yslices*8+j*8+l]=ds[i*yslices*8+j*8+l]+mup*sin(tharray[k])*x[i*yslices*costhslices*8+j*costhslices*8+k*8+l]*dtharray[k];
                    	dc[i*yslices*8+j*8+l]=dc[i*yslices*8+j*8+l]+mup*cos(tharray[k])*x[i*yslices*costhslices*8+j*costhslices*8+k*8+l]*dtharray[k];
		    }
                }
            }
        }
    }
} //pragma loop ends here
//int th_id;
//int nthreads;
#pragma omp parallel private(th_id)
{
    nthreads=omp_get_num_threads();
    th_id = omp_get_thread_num();
    //std::cout<<"nthread = "<<nthreads<<"th_id = "<<th_id<<endl;
    //for(int i=0;i<xslices;i++)
    for(int ip=0;(ip+th_id)<xslices;ip=ip+nthreads)
    {
	int i=ip+th_id;
        for(int j=0;j<yslices;j++)
        {
	    for(int k=0;k<costhslices;k++)
            {
		if(x[i*yslices*costhslices*8+j*costhslices*8+k*8+0]+x[i*yslices*costhslices*8+j*costhslices*8+k*8+1]+x[i*yslices*costhslices*8+j*costhslices*8+k*8+4]+x[i*yslices*costhslices*8+j*costhslices*8+k*8+5] > 1.0e-12)
                {

// Vacuum Hamiltonian
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+0]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+0]+2.0e0*omega*x[i*yslices*costhslices*8+j*costhslices*8+k*8+3]*sin(2.0e0*theta);
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+1]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+1]+-2.0e0*omega*x[i*yslices*costhslices*8+j*costhslices*8+k*8+3]*sin(2.0e0*theta);
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+2]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+2]+2.0e0*omega*x[i*yslices*costhslices*8+j*costhslices*8+k*8+3]*cos(2.0e0*theta);
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+3]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+3]+omega*(-x[i*yslices*costhslices*8+j*costhslices*8+k*8+0]*sin(2.0e0*theta) + x[i*yslices*costhslices*8+j*costhslices*8+k*8+1]*sin(2.0e0*theta) - 2.0e0*x[i*yslices*costhslices*8+j*costhslices*8+k*8+2]*cos(2.0e0*theta));
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+4]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+4]+-2.0e0*omega*x[i*yslices*costhslices*8+j*costhslices*8+k*8+7]*sin(2.0e0*theta);
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+5]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+5]+2.0e0*omega*x[i*yslices*costhslices*8+j*costhslices*8+k*8+7]*sin(2.0e0*theta);
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+6]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+6]+-2.0e0*omega*x[i*yslices*costhslices*8+j*costhslices*8+k*8+7]*cos(2.0e0*theta);
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+7]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+7]+omega*(x[i*yslices*costhslices*8+j*costhslices*8+k*8+4]*sin(2.0e0*theta) - x[i*yslices*costhslices*8+j*costhslices*8+k*8+5]*sin(2.0e0*theta) + 2.0e0*x[i*yslices*costhslices*8+j*costhslices*8+k*8+6]*cos(2.0e0*theta));
// Vacuum Hamiltonian
// Matter Hamiltonian
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+0]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+0]+0;
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+1]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+1]+0;
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+2]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+2]+-lam*x[i*yslices*costhslices*8+j*costhslices*8+k*8+3];
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+3]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+3]+lam*x[i*yslices*costhslices*8+j*costhslices*8+k*8+2];
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+4]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+4]+0;
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+5]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+5]+0;
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+6]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+6]+-lam*x[i*yslices*costhslices*8+j*costhslices*8+k*8+7];
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+7]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+7]+lam*x[i*yslices*costhslices*8+j*costhslices*8+k*8+6];
// Matter Hamiltonian
// Self-Interaction Hamiltonian
//------------------------------------------------------//
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+0]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+0]+2.0e0*d[i*yslices*8+j*8+2]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+3] - 2.0e0*d[i*yslices*8+j*8+3]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+2] - 2.0e0*d[i*yslices*8+j*8+6]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+3] + 2.0e0*d[i*yslices*8+j*8+7]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+2] - 2.0e0*cos(tharray[k])*dc[i*yslices*8+j*8+2]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+3] + 2.0e0*cos(tharray[k])*dc[i*yslices*8+j*8+3]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+2] + 2.0e0*cos(tharray[k])*dc[i*yslices*8+j*8+6]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+3] - 2.0e0*cos(tharray[k])*dc[i*yslices*8+j*8+7]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+2] - 2.0e0*sin(tharray[k])*ds[i*yslices*8+j*8+2]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+3] + 2.0e0*sin(tharray[k])*ds[i*yslices*8+j*8+3]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+2] + 2.0e0*sin(tharray[k])*ds[i*yslices*8+j*8+6]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+3] - 2.0e0*sin(tharray[k])*ds[i*yslices*8+j*8+7]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+2];
//------------------------------------------------------//
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+1]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+1]+-2.0e0*d[i*yslices*8+j*8+2]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+3] + 2.0e0*d[i*yslices*8+j*8+3]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+2] + 2.0e0*d[i*yslices*8+j*8+6]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+3] - 2.0e0*d[i*yslices*8+j*8+7]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+2] + 2.0e0*cos(tharray[k])*dc[i*yslices*8+j*8+2]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+3] - 2.0e0*cos(tharray[k])*dc[i*yslices*8+j*8+3]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+2] - 2.0e0*cos(tharray[k])*dc[i*yslices*8+j*8+6]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+3] + 2.0e0*cos(tharray[k])*dc[i*yslices*8+j*8+7]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+2] + 2.0e0*sin(tharray[k])*ds[i*yslices*8+j*8+2]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+3] - 2.0e0*sin(tharray[k])*ds[i*yslices*8+j*8+3]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+2] - 2.0e0*sin(tharray[k])*ds[i*yslices*8+j*8+6]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+3] + 2.0e0*sin(tharray[k])*ds[i*yslices*8+j*8+7]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+2];
//------------------------------------------------------//
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+2]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+2]+x[i*yslices*costhslices*8+j*costhslices*8+k*8+0]*(d[i*yslices*8+j*8+3] - d[i*yslices*8+j*8+7] - cos(tharray[k])*dc[i*yslices*8+j*8+3] + cos(tharray[k])*dc[i*yslices*8+j*8+7] - sin(tharray[k])*ds[i*yslices*8+j*8+3] + sin(tharray[k])*ds[i*yslices*8+j*8+7]) - x[i*yslices*costhslices*8+j*costhslices*8+k*8+1]*(d[i*yslices*8+j*8+3] - d[i*yslices*8+j*8+7] - cos(tharray[k])*dc[i*yslices*8+j*8+3] + cos(tharray[k])*dc[i*yslices*8+j*8+7] - sin(tharray[k])*ds[i*yslices*8+j*8+3] + sin(tharray[k])*ds[i*yslices*8+j*8+7]) - x[i*yslices*costhslices*8+j*costhslices*8+k*8+3]*(d[i*yslices*8+j*8+0] - d[i*yslices*8+j*8+4] - cos(tharray[k])*dc[i*yslices*8+j*8+0] + cos(tharray[k])*dc[i*yslices*8+j*8+4] - sin(tharray[k])*ds[i*yslices*8+j*8+0] + sin(tharray[k])*ds[i*yslices*8+j*8+4]) + x[i*yslices*costhslices*8+j*costhslices*8+k*8+3]*(d[i*yslices*8+j*8+1] - d[i*yslices*8+j*8+5] - cos(tharray[k])*dc[i*yslices*8+j*8+1] + cos(tharray[k])*dc[i*yslices*8+j*8+5] - sin(tharray[k])*ds[i*yslices*8+j*8+1] + sin(tharray[k])*ds[i*yslices*8+j*8+5]);
//------------------------------------------------------//
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+3]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+3]+-x[i*yslices*costhslices*8+j*costhslices*8+k*8+0]*(d[i*yslices*8+j*8+2] - d[i*yslices*8+j*8+6] - cos(tharray[k])*dc[i*yslices*8+j*8+2] + cos(tharray[k])*dc[i*yslices*8+j*8+6] - sin(tharray[k])*ds[i*yslices*8+j*8+2] + sin(tharray[k])*ds[i*yslices*8+j*8+6]) + x[i*yslices*costhslices*8+j*costhslices*8+k*8+1]*(d[i*yslices*8+j*8+2] - d[i*yslices*8+j*8+6] - cos(tharray[k])*dc[i*yslices*8+j*8+2] + cos(tharray[k])*dc[i*yslices*8+j*8+6] - sin(tharray[k])*ds[i*yslices*8+j*8+2] + sin(tharray[k])*ds[i*yslices*8+j*8+6]) + x[i*yslices*costhslices*8+j*costhslices*8+k*8+2]*(d[i*yslices*8+j*8+0] - d[i*yslices*8+j*8+4] - cos(tharray[k])*dc[i*yslices*8+j*8+0] + cos(tharray[k])*dc[i*yslices*8+j*8+4] - sin(tharray[k])*ds[i*yslices*8+j*8+0] + sin(tharray[k])*ds[i*yslices*8+j*8+4]) - x[i*yslices*costhslices*8+j*costhslices*8+k*8+2]*(d[i*yslices*8+j*8+1] - d[i*yslices*8+j*8+5] - cos(tharray[k])*dc[i*yslices*8+j*8+1] + cos(tharray[k])*dc[i*yslices*8+j*8+5] - sin(tharray[k])*ds[i*yslices*8+j*8+1] + sin(tharray[k])*ds[i*yslices*8+j*8+5]);
//------------------------------------------------------//
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+4]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+4]+2.0e0*d[i*yslices*8+j*8+2]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+7] - 2.0e0*d[i*yslices*8+j*8+3]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+6] - 2.0e0*d[i*yslices*8+j*8+6]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+7] + 2.0e0*d[i*yslices*8+j*8+7]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+6] - 2.0e0*cos(tharray[k])*dc[i*yslices*8+j*8+2]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+7] + 2.0e0*cos(tharray[k])*dc[i*yslices*8+j*8+3]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+6] + 2.0e0*cos(tharray[k])*dc[i*yslices*8+j*8+6]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+7] - 2.0e0*cos(tharray[k])*dc[i*yslices*8+j*8+7]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+6] - 2.0e0*sin(tharray[k])*ds[i*yslices*8+j*8+2]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+7] + 2.0e0*sin(tharray[k])*ds[i*yslices*8+j*8+3]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+6] + 2.0e0*sin(tharray[k])*ds[i*yslices*8+j*8+6]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+7] - 2.0e0*sin(tharray[k])*ds[i*yslices*8+j*8+7]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+6];
//------------------------------------------------------//
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+5]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+5]+-2.0e0*d[i*yslices*8+j*8+2]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+7] + 2.0e0*d[i*yslices*8+j*8+3]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+6] + 2.0e0*d[i*yslices*8+j*8+6]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+7] - 2.0e0*d[i*yslices*8+j*8+7]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+6] + 2.0e0*cos(tharray[k])*dc[i*yslices*8+j*8+2]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+7] - 2.0e0*cos(tharray[k])*dc[i*yslices*8+j*8+3]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+6] - 2.0e0*cos(tharray[k])*dc[i*yslices*8+j*8+6]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+7] + 2.0e0*cos(tharray[k])*dc[i*yslices*8+j*8+7]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+6] + 2.0e0*sin(tharray[k])*ds[i*yslices*8+j*8+2]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+7] - 2.0e0*sin(tharray[k])*ds[i*yslices*8+j*8+3]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+6] - 2.0e0*sin(tharray[k])*ds[i*yslices*8+j*8+6]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+7] + 2.0e0*sin(tharray[k])*ds[i*yslices*8+j*8+7]*x[i*yslices*costhslices*8+j*costhslices*8+k*8+6];
//------------------------------------------------------//
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+6]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+6]+x[i*yslices*costhslices*8+j*costhslices*8+k*8+4]*(d[i*yslices*8+j*8+3] - d[i*yslices*8+j*8+7] - cos(tharray[k])*dc[i*yslices*8+j*8+3] + cos(tharray[k])*dc[i*yslices*8+j*8+7] - sin(tharray[k])*ds[i*yslices*8+j*8+3] + sin(tharray[k])*ds[i*yslices*8+j*8+7]) - x[i*yslices*costhslices*8+j*costhslices*8+k*8+5]*(d[i*yslices*8+j*8+3] - d[i*yslices*8+j*8+7] - cos(tharray[k])*dc[i*yslices*8+j*8+3] + cos(tharray[k])*dc[i*yslices*8+j*8+7] - sin(tharray[k])*ds[i*yslices*8+j*8+3] + sin(tharray[k])*ds[i*yslices*8+j*8+7]) - x[i*yslices*costhslices*8+j*costhslices*8+k*8+7]*(d[i*yslices*8+j*8+0] - d[i*yslices*8+j*8+4] - cos(tharray[k])*dc[i*yslices*8+j*8+0] + cos(tharray[k])*dc[i*yslices*8+j*8+4] - sin(tharray[k])*ds[i*yslices*8+j*8+0] + sin(tharray[k])*ds[i*yslices*8+j*8+4]) + x[i*yslices*costhslices*8+j*costhslices*8+k*8+7]*(d[i*yslices*8+j*8+1] - d[i*yslices*8+j*8+5] - cos(tharray[k])*dc[i*yslices*8+j*8+1] + cos(tharray[k])*dc[i*yslices*8+j*8+5] - sin(tharray[k])*ds[i*yslices*8+j*8+1] + sin(tharray[k])*ds[i*yslices*8+j*8+5]);
//------------------------------------------------------//
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+7]=dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+7]+-x[i*yslices*costhslices*8+j*costhslices*8+k*8+4]*(d[i*yslices*8+j*8+2] - d[i*yslices*8+j*8+6] - cos(tharray[k])*dc[i*yslices*8+j*8+2] + cos(tharray[k])*dc[i*yslices*8+j*8+6] - sin(tharray[k])*ds[i*yslices*8+j*8+2] + sin(tharray[k])*ds[i*yslices*8+j*8+6]) + x[i*yslices*costhslices*8+j*costhslices*8+k*8+5]*(d[i*yslices*8+j*8+2] - d[i*yslices*8+j*8+6] - cos(tharray[k])*dc[i*yslices*8+j*8+2] + cos(tharray[k])*dc[i*yslices*8+j*8+6] - sin(tharray[k])*ds[i*yslices*8+j*8+2] + sin(tharray[k])*ds[i*yslices*8+j*8+6]) + x[i*yslices*costhslices*8+j*costhslices*8+k*8+6]*(d[i*yslices*8+j*8+0] - d[i*yslices*8+j*8+4] - cos(tharray[k])*dc[i*yslices*8+j*8+0] + cos(tharray[k])*dc[i*yslices*8+j*8+4] - sin(tharray[k])*ds[i*yslices*8+j*8+0] + sin(tharray[k])*ds[i*yslices*8+j*8+4]) - x[i*yslices*costhslices*8+j*costhslices*8+k*8+6]*(d[i*yslices*8+j*8+1] - d[i*yslices*8+j*8+5] - cos(tharray[k])*dc[i*yslices*8+j*8+1] + cos(tharray[k])*dc[i*yslices*8+j*8+5] - sin(tharray[k])*ds[i*yslices*8+j*8+1] + sin(tharray[k])*ds[i*yslices*8+j*8+5]);
// Self-Interaction Hamiltonian
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+0]=c*dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+0];
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+1]=c*dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+1];
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+2]=c*dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+2];
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+3]=c*dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+3];
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+4]=c*dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+4];
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+5]=c*dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+5];
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+6]=c*dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+6];
dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+7]=c*dxdt[i*yslices*costhslices*8+j*costhslices*8+k*8+7];
		}
            }
        }
    }
}// closing parenthesis for pragma loop
}
void my_observer( const state_type &x, const double t )
{
    std::cout<<t<<" ";
    for(int i=0;i<xslices;i++)
    {
	for( int j=0;j<yslices;j++)
	{
	    for(int k=0;k<costhslices;k=k+10)
	    {
		for(int l=0;l<8;l++)
		{
		    //std::cout<<x[i*yslices*costhslices*8+j*costhslices*8+k*8+l]<<" ";
		}
	    }
	}
    }
    std::cout<<std::endl;
    //std::string fname = __FILE__;
    //size_t pos = fname.find("cpp");
    //boost::replace_all(fname, "cpp", "sum");
    //std::cout<< fname <<endl;
    //ofstream fp;
    //fp.open(fname,ios::app | ios::out);
    for(int i=0;i<xslices;i++)
    {
        for(int j=0;j<yslices;j++)
        {
	    for(int l=0;l<8;l++)
            {
		d[i*yslices*8+j*8+l]=0.0e0;
	    }	
	}
    }
    fp<<t<<" "; 
    for(int i=0;i<xslices;i++)
    {
        for(int j=0;j<yslices;j++)
        {
            for(int k=0;k<costhslices;k++)
            {
		for(int l=0;l<8;l++)
                {
                    d[i*yslices*8+j*8+l]=d[i*yslices*8+j*8+l]+x[i*yslices*costhslices*8+j*costhslices*8+k*8+l];
                }
            }
	    for(int l=0;l<8;l++)
	    {
		fp<<d[i*yslices*8+j*8+l]<<" ";
	    }
        }
    }
    fp<<endl;   
    //fp.close();
}
namespace boost {
namespace numeric {
namespace odeint {
namespace detail {

template< class Stepper , class System , class State , class Time , class Observer >
size_t my_integrate_adaptive(
        Stepper stepper , System system , State &start_state ,
        Time &start_time , Time end_time , Time &dt ,
        Observer observer , controlled_stepper_tag
)
{
    typename odeint::unwrap_reference< Observer >::type &obs = observer;
    typename odeint::unwrap_reference< Stepper >::type &st = stepper;

    failed_step_checker fail_checker;  // to throw a runtime_error if step size adjustment fails
    size_t count = 0;
    double dtold=0.0e0;
    double dtused=0.0e0;
    //cout<<"Start dt = "<<dt<<endl;
    while( less_with_sign( start_time , end_time , dt ) )
    {
        obs( start_state , start_time );
        if( less_with_sign( end_time , static_cast<Time>(start_time + dt) , dt ) )
        {
            dt = end_time - start_time;
        }

        controlled_step_result res;
        do
        {
	    dtold=dt;
            res = st.try_step( system , start_state , start_time , dt );
            fail_checker();  // check number of failed steps
	    //cout<<"Count, dtold and dt = "<<count<<" "<<dtold<<" "<<dt<<endl;
        }
        while( res == fail );
        fail_checker.reset();  // if we reach here, the step was successful -> reset fail checker

        ++count;
	dtused=dtused+dtold;
	//cout<<"Count and dtused = "<<count<<" "<<dtused<<endl;
    }
    obs( start_state , start_time );
    transport(start_state,dtused);
    return count;
}


} // namespace detail
} // namespace odeint
} // namespace numeric
} // namespace boost


namespace boost {
    namespace numeric {
        namespace odeint {


            /*
 *              * the two overloads are needed in order to solve the forwarding problem
 *                           */
            template< class Stepper , class System , class State , class Time , class Observer >
            size_t my_integrate_adaptive(
                                      Stepper stepper , System system , State &start_state ,
                                      Time start_time , Time end_time , Time dt ,
                                      Observer observer )
            {
                typedef typename odeint::unwrap_reference< Stepper >::type::stepper_category stepper_category;
                return detail::my_integrate_adaptive(
                                                  stepper , system , start_state ,
                                                  start_time , end_time , dt ,
                                                  observer , stepper_category() );

                /*
 *                  * Suggestion for a new extendable version:
 *                                   *
 *                                                    * integrator_adaptive< Stepper , System, State , Time , Observer , typename Stepper::stepper_category > integrator;
 *                                                                     * return integrator.run( stepper , system , start_state , start_time , end_time , dt , observer );
 *                                                                                      */
            }
        } // namespace odeint
    } // namespace numeric
} // namespace boost






int main(int argc, char* argv[])
{
    //-----------
    //std::cout<< __FILE__ <<endl;
    std::string fname = __FILE__;
    size_t pos = fname.find("cpp");
    boost::replace_all(fname, "cpp", "sum");
    //std::cout<< fname <<endl;
    //ofstream fp;
    fp.open(fname);
    //-------------
    std::cout.setf ( std::ios::scientific, std::ios::floatfield );
    std::cout.precision(15);
    int size=tsize;
    //printf("%d\n",size);
    state_type x0(size); // Initial condition, vector of 2 elements (position and velocity)
    //state_type dx(size);
    //state_type dy(size);
    //x0[0] = 0.0;
    //x0[1] = 1.0;
    double err_abs = 1.0e-16;
    double err_rel = 1.0e-12;
    double a_x = 1.0;
    double a_dxdt = 1.0;
    initarrays(x0);
    // Integration parameters
    double t0 = 0.0e0; // time in seconds
    int tsteps = 200;
    double t1 = 2.0e-5; // time in seconds
    double dt = (t1-t0)/((double)tsteps);
    //typedef runge_kutta_cash_karp54< state_type > solver;
    //typedef runge_kutta_dopri5< state_type > solver;
    typedef runge_kutta_fehlberg78< state_type > solver;
    typedef controlled_runge_kutta< solver > controller;
    my_observer(x0,t0);
    for(int ts=0;ts<tsteps;ts++)
    {
        //printf("xxxx\n");
        my_integrate_adaptive( make_controlled( err_abs , err_rel , solver() ), my_system, x0 , t0+ts*dt , t0+(1+ts)*dt , dt*1.0e-1,null_observer());
        my_observer(x0,t0+(1+ts)*dt);
    }
    fp.close();
}
