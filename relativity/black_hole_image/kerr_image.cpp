#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
//#include <vector>
#include <fstream>
#include <sstream>
#include <chrono>
#include <omp.h>

//this is a function to create linearly separated points. 
void linspace(double start, double end, double numpoints, double* points) {
	//double * points  = new double [numpoints];
	double dx = (end - start) / (numpoints - 1);
	for (int i = 0; i < numpoints; i++) {
		points[i] = (start + i * dx);
	}
}

//these are metric functions for Kerr
double sigma(double r, double th, double a) {
	return pow(r, 2) + pow(a * cos(th), 2);
}

double delta(double M, double r, double a) {
	return pow(r, 2) + pow(a, 2) - 2 * M * r;
}

double g_tt(double M, double r, double th, double a) {
	return -1 * (1 - (2 * M * r) / (sigma(r, th, a)));
}

double g_tphi(double M, double r, double th, double a) {
	return -(2 * M * r * a * pow(sin(th), 2.0)) / (sigma(r, th, a));
}

double g_rr(double M, double r, double th, double a) {
	return sigma(r, th, a) / delta(M, r, a);
}

double g_thth(double r, double th, double a) {
	return sigma(r, th, a);
}

double g_phiphi(double M, double r, double th, double a) {
	return (pow(a, 2) + pow(r, 2) + (2 * M * r * pow(a * sin(th), 2.0)) / (sigma(r, th, a))) * pow(sin(th), 2);
}

//Slopes for t, phi and the effective_potential equation

double t_slope(double E, double L, double M, double r, double th, double a) {
	return (E * g_phiphi(M, r, th, a) + L * g_tphi(M, r, th, a)) / (pow(g_tphi(M, r, th, a), 2) - g_tt(M, r, th, a) * g_phiphi(M, r, th, a));
}

double phi_slope(double E, double L, double M, double r, double th, double a) {
	return -(E * g_tphi(M, r, th, a) + L * g_tt(M, r, th, a)) / (pow(g_tphi(M, r, th, a), 2) - g_tt(M, r, th, a) * g_phiphi(M, r, th, a));
}

double V_eff(double E, double L, double M, double r, double th, double a) {
	return (pow(g_rr(M, r, th, a), -1)) * ((g_phiphi(M, r, th, a) * pow(E, 2) + g_tt(M, r, th, a) * pow(L, 2) + 2 * E * L * g_tphi(M, r, th, a)) / (-pow(g_tphi(M, r, th, a), 2) + g_tt(M, r, th, a) * g_phiphi(M, r, th, a)));
}

// all the derivatives of the metric required are also coded as functions

//first is for the radial equation of motion
double g_rr_th(double M, double r, double th, double a) {
	return -(sin(2 * th) * pow(a, 2)) / (delta(M, r, a));
}

double g_rr_r(double M, double r, double th, double a) {
	return (2 / pow(delta(M, r, a), 2)) * (r * delta(M, r, a) - (r - M) * sigma(r, th, a));

}

double g_tt_r(double M, double r, double th, double a) {
	return ((2 * M) / (pow(sigma(r, th, a), 2))) * (sigma(r, th, a) - 2 * pow(r, 2));
}

double g_tphi_r(double M, double r, double th, double a) {
	return (2 * M * a * pow(sin(th), 2) / pow(sigma(r, th, a), 2)) * (2 * pow(r, 2) - sigma(r, th, a));
}

double g_phiphi_r(double M, double r, double th, double a) {
	return (pow(sin(th), 2)) * (2 * r + ((2 * M * pow(a * sin(th), 2)) / pow(sigma(r, th, a), 2)) * (sigma(r, th, a) - 2 * pow(r, 2)));
}

double g_thth_r(double r) {
	return 2 * r;
}

//now for the theta equation of motion

double g_tt_th(double M, double r, double th, double a) {
	return (2 * M * r * sin(2 * th) * pow(a, 2)) / pow(sigma(r, th, a), 2);
}

double g_tphi_th(double M, double r, double th, double a) {
	return -((2 * M * r * a * sin(2 * th)) / pow(sigma(r, th, a), 2)) * (sigma(r, th, a) + pow(a * sin(th), 2));
}

double g_thth_th(double th, double a) {
	return -sin(2 * th) * pow(a, 2);
}

double g_phiphi_th(double M, double r, double th, double a) {
	return (sin(2 * th) / pow(sin(th), 2)) * g_phiphi(M, r, th, a) + ((2 * M * r * sin(2 * th) * pow(a * sin(th), 2)) / pow(sigma(r, th, a), 2)) * (sigma(r, th, a) + pow(a * sin(th), 2));
}

double v_r_slope(double M, double r, double th, double a, double v_t, double v_r, double v_th, double v_phi) {
	return  (-1 / (2 * g_rr(M, r, th, a))) * (2 * v_r * v_th * g_rr_th(M, r, th, a) + pow(v_r, 2) * g_rr_r(M, r, th, a) - pow(v_t, 2) * g_tt_r(M, r, th, a) - 2 * v_t * v_phi * g_tphi_r(M, r, th, a) - pow(v_th, 2) * g_thth_r(r) - pow(v_phi, 2) * g_phiphi_r(M, r, th, a));
}

double v_th_slope(double M, double r, double th, double a, double v_t, double v_r, double v_th, double v_phi) {
	return (-1 / (2 * g_thth(r, th, a))) * (2 * v_r * v_th * g_thth_r(r) - pow(v_r, 2) * g_rr_th(M, r, th, a) - pow(v_t, 2) * g_tt_th(M, r, th, a) - 2 * v_t * v_phi * g_tphi_th(M, r, th, a) + pow(v_th, 2) * g_thth_th(th, a) - pow(v_phi, 2) * g_phiphi_th(M, r, th, a));
}

class Kerr_potential_test_photon {

public:
	double mass; //here mass is the mass of the black hole
	double angmom; //here angmom is the angmom per unit mass(a) of the black hole
	//the test particle has unit mass
	double E; //E is the conserved energy set later on
	double L; //L is the conserved angular momentum set later on
	double Q;
	double r_event_min; //inner event horizon
	double r_event_max; //outer event horizon
	double r_plus; //outermost spherical equatorial orbit
	double r_minus; //innermost spherical equatorial orbit
	
	double energy_checker;
	double carter_constant;
	int size;
	bool inspiral;
	bool disk;
	bool escape;

	double g;


	//Store trajectory from RK Cash-Karp integration

	double t;
	double r;
	double th;
	double phi;
	double v_t;
	double v_r;
	double v_th;
	double v_phi;

	double x;
	double y;
	double z;
	double zold;

	Kerr_potential_test_photon(double M, double a)
	{
		mass = M;
		angmom = a;
		E = 0.0;
		L = 0.0;
		Q = 0.0;
		r_event_min = mass - sqrt(pow(mass, 2) - pow(angmom, 2));
		r_event_max = mass + sqrt(pow(mass, 2) - pow(angmom, 2));
		r_plus = 2 * mass * (1 + cos((2.0 / 3.0) * acos(angmom / mass)));
		r_minus = 2 * mass * (1 + cos((2.0 / 3.0) * acos(-angmom / mass)));
	
		energy_checker = 0.0;
		carter_constant  = 0.0;
		size = 0;
		inspiral = false;
		disk = false;
		escape = false;

		g = 0.;

		t = 0.0;
		r = 0.0;
		th = 0.0;
		phi = 0.0;
		v_t = 0.0;
		v_r = 0.0;
		v_th = 0.0;
		v_phi = 0.0;

		x = 0.0;
		y = 0.0;
		z = 0.0;
		zold = 0.0;

	}

	void set_initial_data(double t0, double r0, double th0, double phi0, double En, double Lz, double v_r_0, double v_th_0) {
		
		//re-initialize the stats for each new photon
		E = En;
		L = Lz;
		Q =  pow(sigma(r0,th0,angmom)*v_th_0,2) + pow(cos(th0),2)*(-pow(E*angmom,2) + pow(L/sin(th0),2));
		
		energy_checker =  fabs(pow(v_r_0, 2.0) + pow(v_th_0, 2.0) * g_thth(r0,th0,angmom) / g_rr(mass,r0,th0,angmom) + V_eff(E,L,mass,r0,th0,angmom));
		carter_constant = carter_constant + Q;

		size = 1;
		inspiral = false;
		disk = false;
		escape = false;
		
		t = t0;
		r = r0;
		th = th0;
		phi = phi0;

		v_t = t_slope(E, L, mass, r0, th0, angmom);
		v_phi = phi_slope(E, L, mass, r0,th0,angmom);
		v_r = v_r_0;
		v_th = v_th_0;

		x = (sqrt(pow(r0,2) + pow(angmom,2)) * sin(th0) * cos(phi0));
		y = (sqrt(pow(r0,2) + pow(angmom,2)) * sin(th0) * sin(phi0));
		z = (r0 * cos(th0));
		zold = z;

	}



	void RKCK(double stop) {
		double tau = 0.0;
		double dtau = - pow(10.,-3.);

		double rtol = pow(10.,-11.);
		double atol = pow(10.,-11.);
		double err;

		double dtau_min = -pow(10.,-6.);
		double dtau_max = -pow(10.,-2.);
		
		double r_d = 0.0;
		double th_d = 0.0;
		double v_t_d = 0.0;
		double v_r_d = 0.0;
		double v_th_d = 0.0;
		double v_phi_d = 0.0;

		double v_r_k1;
		double r_k1;
		double v_th_k1;
		double th_k1;
		double t_k1;
		double phi_k1;

		double v_r_k2;
		double r_k2;
		double v_th_k2;
		double th_k2;
		double t_k2;
		double phi_k2;

		double v_r_k3;
		double r_k3;
		double v_th_k3;
		double th_k3;
		double t_k3;
		double phi_k3;

		double v_r_k4;
		double r_k4;
		double v_th_k4;
		double th_k4;
		double t_k4;
		double phi_k4;

		double v_r_k5;
		double r_k5;
		double v_th_k5;
		double th_k5;
		double t_k5;
		double phi_k5;

		double v_r_k6;
		double r_k6;
		double v_th_k6;
		double th_k6;
		double t_k6;
		double phi_k6;

		//lower-order solution
		double t_5;
		double r_5;
		double th_5;
		double phi_5;
		double v_t_5;
		double v_r_5;
		double v_th_5;
		double v_phi_5;

		
		
		while (true) {
			/*
			if (r<=r_event_max || (th>=-0.15 && th<=0.15) || (th>=M_PI-0.15 && th<=M_PI+0.15)){
				dtau = -pow(10,-3);
			}
			else if (r_event_max < r && r <= 10.0*r_event_max){
				dtau = - pow(10,-2);
			}
			else if (10.0*r_event_max < r && r <= 50.0*r_event_max){
				dtau = -5.0*pow(10,-2);
			}
			else{
				dtau = -pow(10,-1);
			}*/
			v_r_k1 = (dtau)*v_r_slope(mass,r, th, angmom, v_t, v_r, v_th, v_phi);
			r_k1 = (dtau)*v_r;
			v_th_k1 = (dtau)*v_th_slope(mass,r, th, angmom, v_t, v_r, v_th, v_phi);
			th_k1 = (dtau)*v_th;
			t_k1 = (dtau)*v_t;
			phi_k1 = (dtau)*v_phi;

			r_d = r + r_k1 / 5.0;


			if (r_d <= r_event_max) {
				inspiral = true;
				break;
			}

			if (r_d>=stop){
				escape = true;
				break;
			}

			th_d = th + th_k1 / 5.0;
			
			v_r_d = v_r + v_r_k1 / 5.0;
			v_th_d = v_th + v_th_k1 / 5.0;
			v_t = t_slope(E,L,mass,r_d, th_d, angmom);
			v_phi = phi_slope(E,L,mass,r_d, th_d,angmom);


			v_r_k2 = (dtau)*v_r_slope(mass,r_d, th_d, angmom, v_t, v_r_d, v_th_d, v_phi);
			r_k2 = (dtau) * (v_r_d);
			v_th_k2 = (dtau)*v_th_slope(mass,r_d, th_d, angmom, v_t, v_r_d, v_th_d, v_phi);
			th_k2 = (dtau) * (v_th_d);
			t_k2 = (dtau)*v_t;
			phi_k2 = (dtau)*v_phi;

			r_d = r + r_k1 * (3.0 / 40.0) + r_k2 * (9.0 / 40.0);


			if (r_d <= r_event_max) {
				inspiral = true;
				break;
			}

			if (r_d>=stop){
				escape = true;
				break;
			}

			th_d = th + th_k1 * (3.0 / 40.0) + th_k2 * (9.0 / 40.0);
			
			v_r_d = v_r + v_r_k1 * (3.0 / 40.0) + v_r_k2 * (9.0 / 40.0);
			v_th_d = v_th + v_th_k1 * (3.0 / 40.0) + v_th_k2 * (9.0 / 40.0);
			v_t = t_slope(E,L,mass,r_d, th_d,angmom);
			v_phi = phi_slope(E,L,mass,r_d, th_d, angmom);

			v_r_k3 = (dtau)*v_r_slope(mass,r_d, th_d,angmom, v_t, v_r_d, v_th_d, v_phi);
			r_k3 = (dtau) * (v_r_d);
			v_th_k3 = (dtau)*v_th_slope(mass,r_d, th_d, angmom, v_t, v_r_d, v_th_d, v_phi);
			th_k3 = (dtau) * (v_th_d);
			t_k3 = (dtau)*v_t;
			phi_k3 = (dtau)*v_phi;

			r_d = r + r_k1 * (3.0 / 10.0) + r_k2 * (-9.0 / 10.0) + r_k3 * (6.0 / 5.0);


			if (r_d <= r_event_max) {
				inspiral = true;
				break;
			}

			if (r_d>=stop){
				escape = true;
				break;
			}

			th_d = th + th_k1 * (3.0 / 10.0) + th_k2 * (-9.0 / 10.0) + th_k3 * (6.0 / 5.0);
			v_r_d = v_r + v_r_k1 * (3.0 / 10.0) + v_r_k2 * (-9.0 / 10.0) + v_r_k3 * (6.0 / 5.0);
			v_th_d = v_th + v_th_k1 * (3.0 / 10.0) + v_th_k2 * (-9.0 / 10.0) + v_th_k3 * (6.0 / 5.0);
			v_t = t_slope(E,L,mass,r_d, th_d, angmom);
			v_phi = phi_slope(E,L,mass,r_d, th_d, angmom);

			v_r_k4 = (dtau)*v_r_slope(mass,r_d, th_d,angmom, v_t, v_r_d, v_th_d, v_phi);
			r_k4 = (dtau) * (v_r_d);
			v_th_k4 = (dtau)*v_th_slope(mass,r_d, th_d, angmom,v_t, v_r_d, v_th_d, v_phi);
			th_k4 = (dtau) * (v_th_d);
			t_k4 = (dtau)*v_t;
			phi_k4 = (dtau)*v_phi;

			r_d = r+ r_k1 * (-11.0 / 54.0) + r_k2 * (5.0 / 2.0) + r_k3 * (-70.0 / 27.0) + r_k4 * (35.0 / 27.0);


			if (r_d <= r_event_max) {
				inspiral = true;
				break;
			}

			if (r_d>=stop){
				escape = true;
				break;
			}

			th_d = th + th_k1 * (-11.0 / 54.0) + th_k2 * (5.0 / 2.0) + th_k3 * (-70.0 / 27.0) + th_k4 * (35.0 / 27.0);
			v_r_d = v_r + v_r_k1 * (-11.0 / 54.0) + v_r_k2 * (5.0 / 2.0) + v_r_k3 * (-70.0 / 27.0) + v_r_k4 * (35.0 / 27.0);
			v_th_d = v_th + v_th_k1 * (-11.0 / 54.0) + v_th_k2 * (5.0 / 2.0) + v_th_k3 * (-70.0 / 27.0) + v_th_k4 * (35.0 / 27.0);
			v_t = t_slope(E,L,mass,r_d, th_d, angmom);
			v_phi = phi_slope(E,L,mass,r_d, th_d, angmom);

			v_r_k5 = (dtau)*v_r_slope(mass, r_d, th_d,angmom, v_t, v_r_d, v_th_d, v_phi);
			r_k5 = (dtau) * (v_r_d);
			v_th_k5 = (dtau)*v_th_slope(mass,r_d, th_d,angmom, v_t, v_r_d, v_th_d, v_phi);
			th_k5 = (dtau) * (v_th_d);
			t_k5 = (dtau)*v_t;
			phi_k5 = (dtau)*v_phi;

			r_d = r + r_k1 * (1631.0 / 55296.0) + r_k2 * (175.0 / 512.0) + r_k3 * (575.0 / 13824.0) + r_k4 * (44275.0 / 110592.0) + r_k5 * (253.0 / 4096.0);


			if (r_d <= r_event_max) {
				inspiral = true;
				break;
			}

			if (r_d>=stop){
				escape = true;
				break;
			}

			th_d = th + th_k1 * (1631.0 / 55296.0) + th_k2 * (175.0 / 512.0) + th_k3 * (575.0 / 13824.0) + th_k4 * (44275.0 / 110592.0) + th_k5 * (253.0 / 4096.0);
			v_r_d = v_r + v_r_k1 * (1631.0 / 55296.0) + v_r_k2 * (175.0 / 512.0) + v_r_k3 * (575.0 / 13824.0) + v_r_k4 * (44275.0 / 110592.0) + v_r_k5 * (253.0 / 4096.0);
			v_th_d = v_th + v_th_k1 * (1631.0 / 55296.0) + v_th_k2 * (175.0 / 512.0) + v_th_k3 * (575.0 / 13824.0) + v_th_k4 * (44275.0 / 110592.0) + v_th_k5 * (253.0 / 4096.0);
			v_t = t_slope(E,L,mass,r_d, th_d, angmom);
			v_phi = phi_slope(E,L,mass,r_d, th_d, angmom);

			v_r_k6 = (dtau)*v_r_slope(mass,r_d, th_d, angmom, v_t, v_r_d, v_th_d, v_phi);
			r_k6 = (dtau) * (v_r_d);
			v_th_k6 = (dtau)*v_th_slope(mass,r_d, th_d, angmom,v_t, v_r_d, v_th_d, v_phi);
			th_k6 = (dtau) * (v_th_d);
			t_k6 = (dtau)*v_t;
			phi_k6 = (dtau)*v_phi;

			
			if (r + (37.0 / 378.0) * r_k1 + (250.0 / 621.0) * r_k3 + (125.0 / 594.0) * r_k4 + (512.0 / 1771.0) * r_k6 <= r_event_max) {
				inspiral = true;
				break;
			}

			if (r + (37.0 / 378.0) * r_k1 + (250.0 / 621.0) * r_k3 + (125.0 / 594.0) * r_k4 + (512.0 / 1771.0) * r_k6 >= stop) {
				escape = true;
				break;
			}


			//Order 6 values/estimate of true value
			
			t = t + (37.0 / 378.0) * t_k1 + (250.0 / 621.0) * t_k3 + (125.0 / 594.0) * t_k4 + (512.0 / 1771.0) * t_k6;
			r = r + (37.0 / 378.0) * r_k1 + (250.0 / 621.0) * r_k3 + (125.0 / 594.0) * r_k4 + (512.0 / 1771.0) * r_k6;
			th = th + (37.0 / 378.0) * th_k1 + (250.0 / 621.0) * th_k3 + (125.0 / 594.0) * th_k4 + (512.0 / 1771.0) * th_k6;
			phi = phi + (37.0 / 378.0) * phi_k1 + (250.0 / 621.0) * phi_k3 + (125.0 / 594.0) * phi_k4 + (512.0 / 1771.0) * phi_k6;
			v_t = t_slope(E,L,mass,r, th, angmom);
			v_phi = phi_slope(E,L,mass,r, th,angmom);
			v_r = v_r + (37.0 / 378.0) * v_r_k1 + (250.0 / 621.0) * v_r_k3 + (125.0 / 594.0) * v_r_k4 + (512.0 / 1771.0) * v_r_k6;
			v_th = v_th + (37.0 / 378.0) * v_th_k1 + (250.0 / 621.0) * v_th_k3 + (125.0 / 594.0) * v_th_k4 + (512.0 / 1771.0) * v_th_k6;
			
			//converting to cartesian coordinates
			x = (sqrt(pow(r,2) + pow(angmom,2)) * sin(th) * cos(phi));
			y = (sqrt(pow(r,2) + pow(angmom,2)) * sin(th) * sin(phi));
			z = (r * cos(th));

			tau = tau + dtau;
			size = size + 1;

			//energy error
			energy_checker = energy_checker + fabs(pow(v_r, 2.0) + pow(v_th, 2.0) * g_thth(r, th, angmom) / g_rr(mass, r, th, angmom) + V_eff(E, L, mass, r, th, angmom));
			carter_constant = carter_constant + pow(sigma(r,th,angmom)*v_th,2) + pow(cos(th),2)*(-pow(E*angmom,2) + pow(L/sin(th),2));
			
			//change time step for next step
			t_5 = t + (2825.0 / 27648.0) * t_k1 + (18575.0 / 48384.0) * t_k3 + (13525.0 / 55296.0) * t_k4 + (277.0/14336.0)* t_k5 + (1.0 / 4.0) * t_k6;
			r_5 = r + (2825.0 / 27648.0) * r_k1 + (18575.0 / 48384.0) * r_k3 + (13525.0 / 55296.0) * r_k4 + (277.0/14336.0)* r_k5 + (1.0 / 4.0) * r_k6;
			th_5 = th + (2825.0 / 27648.0) * th_k1 + (18575.0 / 48384.0) * th_k3 + (13525.0 / 55296.0) * th_k4 + (277.0/14336.0)* th_k5 + (1.0 / 4.0) * th_k6;
			phi_5 = phi + (2825.0 / 27648.0) * phi_k1 + (18575.0 / 48384.0) * phi_k3 + (13525.0 / 55296.0) * phi_k4 + (277.0/14336.0)* phi_k5 + (1.0 / 4.0) * phi_k6;
			v_t_5 = t_slope(E,L,mass,r_5, th_5, angmom);
			v_phi_5 = phi_slope(E,L,mass,r_5, th_5,angmom);
			v_r_5 = v_r + (2825.0 / 27648.0) * v_r_k1 + (18575.0 / 48384.0) * v_r_k3 + (13525.0 / 55296.0) * v_r_k4 + (277.0/14336.0)* v_r_k5 + (1.0 / 4.0) * v_r_k6;
			v_th_5 = v_th + (2825.0 / 27648.0) * v_th_k1 + (18575.0 / 48384.0) * v_th_k3 + (13525.0 / 55296.0) * v_th_k4 + (277.0/14336.0)* v_th_k5 + (1.0 / 4.0) * v_th_k6;
			
			err = sqrt((1./8.)*( pow(abs(t - t_5)/(atol + rtol*abs(t - t_5)), 2.) + pow(abs(r - r_5)/(atol + rtol*abs(r - r_5)), 2.) + pow(abs(th - th_5)/(atol + rtol*abs(th - th_5)), 2.) + pow(abs(phi - phi_5)/(atol + rtol*abs(phi - phi_5)), 2.) + pow(abs(v_t - v_t_5)/(atol + rtol*abs(v_t - v_t_5)), 2.) + pow(abs(v_r - v_r_5)/(atol + rtol*abs(v_r - v_r_5)), 2.) + pow(abs(v_th - v_th_5)/(atol + rtol*abs(v_th - v_th_5)), 2.) + pow(abs(v_phi - v_phi_5)/(atol + rtol*abs(v_phi - v_phi_5)), 2.) ));
			
			dtau = dtau*pow( 1./err ,1./5.);

			if (dtau>dtau_min){
				dtau = dtau_min;
			}
			else if (dtau<dtau_max){
				dtau = dtau_max;
			}


			//check if photon hits disk
			if( (6.0*mass<=r && r<=15.0*mass) && ( zold*z<0)){
				disk = true;
				break;
			}
			else{
				zold = z;
			}
		}
	}

	
		
	

	int colour_coder() {
        int colour;
		

		if (inspiral) {
			colour = 1;
		}

		else if (disk) {
			colour = 0;
			g = E/(E/sqrt(1 - 3*mass/r) - L*sqrt(mass/pow(r,3.))/sqrt(1 - 3*mass/r));
		}

		else if (escape){
			colour = 1;
		}
        return colour;

	}

};
int main() {


	double M = 1.0;
	double J = 0.0;
	double a = J/M;
	//G = c = 1, otherwise my life is difficult
	int angle_range = 1000;
	//this is the OSIRIS parametrization
	double* alpha = new double[angle_range];
	double* beta = new double[angle_range];
	
	double lim_alpha = (30.)*(M_PI/180.);
	double lim_beta =  (52.5)*(M_PI/180.);
	linspace(-lim_alpha, lim_alpha, angle_range, alpha);
	linspace(-lim_beta, lim_beta, angle_range, beta);

	std::ofstream myFile("schwarz_image_data_20M_2.txt");
	
	double phi = 0.0;
	double th = (87.5)*(M_PI/180.);
	double r = 20.0*M;
	double mod_P = 1.0;
	
	#pragma omp parallel
	{	
		std::stringstream buffer;
		#pragma omp for
		for (int i = 0; i < angle_range; i++) {
			for (int j = 0; j < angle_range; j++) {
				
				Kerr_potential_test_photon photon(M, a);
		
				double x;
				double y;
				double A_t;
				double P_r_obs;
				double P_th_obs;
				double P_phi_obs;

				double E;
				double L;
				double v_r_0;
				double v_th_0;

				x = - sin(beta[i])*r;
				y = sin(alpha[j])*r;
				P_r_obs = mod_P*cos(alpha[j])*cos(beta[i]);
				P_th_obs = mod_P*sin(alpha[j]);
				P_phi_obs = mod_P*cos(alpha[j])*sin(beta[i]);
				A_t = sqrt(g_phiphi(M,r,th,a)/(pow(g_tphi(M,r,th,a),2) - g_tt(M,r,th,a)*g_phiphi(M,r,th,a)));

				E = (mod_P/A_t) - P_phi_obs*g_tphi(M,r,th,a)/sqrt(g_phiphi(M,r,th,a));
				L = P_phi_obs*sqrt(g_phiphi(M,r,th,a));
				v_r_0 = P_r_obs/sqrt(g_rr(M,r,th,a));
				v_th_0 = P_th_obs/sqrt(g_thth(r,th,a));
				
				
				//set initial data, check if its in allowed region
				photon.set_initial_data(0.0, r, th, phi, E, L, v_r_0, v_th_0);
				
				//do the RKCK method of and error analysis
				photon.RKCK(20*M);

				buffer << x << "," << y << "," << photon.colour_coder() << ","<<photon.r<< ","<<photon.th<<","<<photon.size<<","<<photon.energy_checker/photon.size << ","<<fabs(photon.Q - photon.carter_constant/photon.size)<<","<<photon.g <<"\n";
			}
		}
		#pragma omp critical
		{
			myFile << buffer.str();
			buffer.str("");   // clear
		}
			
	}
	
	myFile.close();

	delete [] alpha;
	delete [] beta;

}