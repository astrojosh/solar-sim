/*
 ============================================================================
 Name        : SolarSystem.c
 Author      : Josh Morrow
 Version     :
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Define constants
#define G 6.67430e-11
#define AU 1.495978707e11
#define DAY 86400
#define YEAR 365.25*DAY

//Function to calculate acceleration between two bodies given the mass of the other body and the x,y positions of each body
double accn(double m_j, double p_j, double p_i, double q_j, double q_i) {

	//Distance between the two bodies
	double r = sqrt(pow((p_j - p_i), 2) + pow((q_j - q_i), 2));

	//Acceleration on first body due to other body
	return G * m_j * (p_j - p_i) / pow(r, 3);

}


int main(void) {

	setvbuf(stdout,NULL,_IONBF,0);
	setvbuf(stderr,NULL,_IONBF,0);

	//Loop variables
	int i, j;
	double t, dt = 0.01*DAY, endTime = 170*YEAR;
	double flag[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

	//Number of bodies
	int n = 10;

	//Runga-Kutta variables
	double ax[12], avx[12], bx[12], bvx[12], cx[12], cvx[12], dx[12], dvx[12];
	double ay[12], avy[12], by[12], bvy[12], cy[12], cvy[12], dy[12], dvy[12];


	//Arrays of constants
	//Legend: {Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto, Asteroid, Moon}

	//Mass of the body
	double m[12] = {1.98892e30, 3.3011e23, 4.8675e24, 5.9742e24, 6.4171e23, 1.8986e27, 5.6834e26, 8.6810e25, 1.02413e26, 1.303e22, 9.5e20, 7.342e22};

	//Initial distance between Sun and body in the x direction
	double x[12] = {0, 0, 0.723332*AU, 0, -1.523679*AU, 0, 9.5826*AU, 0, -30.11*AU, 0, 2.50525*AU, 0};

	//Initial distance between Sun and body in the y direction
	double y[12] = {0, 0.387098*AU, 0, -1.0*AU, 0, 5.2044*AU, 0, -19.2184*AU, 0, 39.482*AU, 0, -1.00257*AU};

	//Initial speed of the body in the x direction
	double vx[12] = {0, sqrt(G*m[0]/y[1]), 0, -sqrt(-G*m[0]/y[3]), 0, sqrt(G*m[0]/y[5]), 0, -sqrt(-G*m[0]/y[7]), 0, sqrt(G*m[0]/y[9]), 0, -sqrt(-G*m[0]/y[11])};

	//Initial speed of the body in the y direction
	double vy[12] = {0, 0, -sqrt(G*m[0]/x[2]), 0, sqrt(-G*m[0]/x[4]), 0, -sqrt(G*m[0]/x[6]), 0, sqrt(-G*m[0]/x[8]), 0, -sqrt(G*m[0]/x[10]), 0};


	//Output file
	FILE *finout, *finout2;
	finout = fopen("SolarSystem.txt","w");
	finout2 = fopen("Asteroid.txt","w");

	//Run code, advance time step, repeat
	for(t=0; t<=endTime; t+=dt) {

		//Calculate the a values for each body at this time step
		//for(i=0; i<n; i++) {
		for(i=1; i<n; i++) {

			avx[i] = 0, avy[i] = 0;

			ax[i] = vx[i];
			ay[i] = vy[i];

			for(j=0; j<n; j++) { if(i != j) {

				avx[i] += accn(m[j], x[j], x[i], y[j], y[i]);
				avy[i] += accn(m[j], y[j], y[i], x[j], x[i]);

			} }

		}

		//MERGE LOOPS?

		//Calculate the b values for each body at this time step
		for(i=1; i<n; i++) {

			bvx[i] = 0, bvy[i] = 0;

			bx[i] = (vx[i] + ((dt/2.0) * avx[i]));
			by[i] = (vy[i] + ((dt/2.0) * avy[i]));

			for(j=0; j<n; j++) { if(i != j) {

				bvx[i] += accn(m[j], x[j], (x[i] + ((dt/2.0) * ax[i])), y[j], (y[i] + ((dt/2.0) * ay[i])));
				bvy[i] += accn(m[j], y[j], (y[i] + ((dt/2.0) * ay[i])), x[j], (x[i] + ((dt/2.0) * ax[i])));

			} }

		}

		//Calculate the c values for each body at this time step
		for(i=1; i<n; i++) {

			cvx[i] = 0, cvy[i] = 0;

			cx[i] = (vx[i] + ((dt/2.0) * bvx[i]));
			cy[i] = (vy[i] + ((dt/2.0) * bvy[i]));

			for(j=0; j<n; j++) { if(i != j) {

				cvx[i] += accn(m[j], x[j], (x[i] + ((dt/2.0) * bx[i])), y[j], (y[i] + ((dt/2.0) * by[i])));
				cvy[i] += accn(m[j], y[j], (y[i] + ((dt/2.0) * by[i])), x[j], (x[i] + ((dt/2.0) * bx[i])));

			} }

		}

		//Calculate the d values for each body at this time step
		for(i=1; i<n; i++) {

			dvx[i] = 0, dvy[i] = 0;

			dx[i] = (vx[i] + (dt * cvx[i]));
			dy[i] = (vy[i] + (dt * cvy[i]));

			for(j=0; j<n; j++) { if(i != j) {

				dvx[i] += accn(m[j], x[j], (x[i] + (dt * cx[i])), y[j], (y[i] + (dt * cy[i])));
				dvy[i] += accn(m[j], y[j], (y[i] + (dt * cy[i])), x[j], (x[i] + (dt * cx[i])));

			} }

		}

		//Calculate the speed and position of each body at this time step
		for(i=1; i<n; i++) {

			//Current x position of the body
			double x_old = x[i];
			double y_old = y[i];

			//New positions and speeds of the body
			x[i] = x[i] + (dt / 6.0) * (ax[i] + 2*bx[i] + 2*cx[i] + dx[i]);
			y[i] = y[i] + (dt / 6.0) * (ay[i] + 2*by[i] + 2*cy[i] + dy[i]);
			vx[i] = vx[i] + (dt / 6.0) * (avx[i] + 2*bvx[i] + 2*cvx[i] + dvx[i]);
			vy[i] = vy[i] + (dt / 6.0) * (avy[i] + 2*bvy[i] + 2*cvy[i] + dvy[i]);

			//Check if the planet has completed a full rotation for the first time
			//if((x_old < 0) && (x[i] > 0) && (flag[i] == 0)) {
			if((x_old * x[i] < 0) || (y_old * y[i] < 0)) {

				flag[i] += 1;
				double year;

				if(i % 2 == 1) {
					year = t - dt * (fabs(x[i]))/(fabs(x[i]) + fabs(x_old));
				}
				else {
					year = t - dt * (fabs(y[i]))/(fabs(y[i]) + fabs(y_old));
				}
				//printf("Body %d Year = %g days\n", i, year/DAY);
				if(flag[i] == 4){
					printf("%d\t%g\n", i, year/DAY);
					fprintf(finout2, "%d\t%g\n", i, year/DAY);
				}

			}

			//Only print values to file every 20 days
			if(fmod(t, 20.0*DAY) == 0) {
				//Print x, y, z and planet (according to legend above)
				fprintf(finout, "%g\t%g\t%d\t%d\n", x[i], y[i], 0, i);

				//Print out kinetic energy of the asteroid
				//LOOP ASTEROID RADII, HISTOGRAM BIN
				if(i == 10) {
					double ke = 0.5 * m[i] * (pow(vx[i], 2) + pow(vy[i], 2));
					//fprintf(finout2, "%g\n", ke);
				}
			}

		}

	}

	//fprintf(finout, "0\t0\t0\n");
	fclose(finout);

	return 0;

}
