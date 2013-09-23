#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace std;

int pbc_round(double input)
{
	int i  = input;

	if (abs(input - i) >= 0.5)
	{
		if (input > 0) {i += 1;}
		if (input < 0) {i -= 1;}
	}
return i;
}

double min_distance(double array[], int length)
{
double min = array[0];
for (int i = 1; i < length; i ++)
{
	if (array[i] < min && array[i] != 0)
	{
		min = array[i];
	}
}
return min;
}


int main()
{
// file streams
ifstream inputfile;
ofstream hbonds_outputfile;

string infile;
int timesteps, nooa, noha;
double lattice;

cout << "Please enter the filename of your file: ";
cin >> infile;
cout << "Please enter the number of oxygen atoms: ";
cin >> nooa;
cout << "Please enter the number of hydrogen atoms: ";
cin >> noha;
cout << "Please enter the number of timesteps: ";
cin >> timesteps;
cout << "Please enter the lattice constant for your periodic cube: ";
cin >> lattice;
cout << "Program running...please wait a moment.\n\n";
// end of main menu*/

//read the inputfile
inputfile.open(infile.c_str());

int noa = nooa + noha;

double x[noa*timesteps], y[noa*timesteps], z[noa*timesteps];
int counter = 0;

while (!inputfile.eof())
{
	inputfile >> x[counter] >> y[counter] >> z[counter];	
	counter ++;

}


double oxyz[nooa][3], hxyz[noha][3];
int ohindices[nooa][4];

for (int i = 0; i < nooa; i ++)
{
	for (int j = 0; j < 4; j ++)
	{
		ohindices[i][j] = -1;
	}
}

vector <double> everything;
hbonds_outputfile.open("hbonds_contour.dat");

int xcoord, ycoord, zcoord, xyzbin[13][13][13] = {};
double xyzhbin[13][13][13] = {}, hcount[nooa] ;

for (int i = 0; i < timesteps; i ++)
{
	for (int j = 0; j < nooa; j ++)
	{
		oxyz[j][0] = lattice*x[j + i*noa];
		oxyz[j][1] = lattice*y[j + i*noa];
		oxyz[j][2] = lattice*z[j + i*noa];
		hcount[j] = 0;
	}
	for (int j = 0; j < noha; j ++)
	{
		hxyz[j][0] = lattice*x[nooa + j + i*noa];
		hxyz[j][1] = lattice*y[nooa + j + i*noa];
		hxyz[j][2] = lattice*z[nooa + j + i*noa];
	}
	for (int j = 0; j < nooa; j++)
	{
		int count = 0;

		for (int k = 0; k < noha; k ++)
		{
			double dx = oxyz[j][0] - hxyz[k][0];
			double dy = oxyz[j][1] - hxyz[k][1];
			double dz = oxyz[j][2] - hxyz[k][2];
	
			dx -= lattice*pbc_round(dx/lattice);
			dy -= lattice*pbc_round(dy/lattice);
			dz -= lattice*pbc_round(dz/lattice);			
	
			double ohdist = sqrt ( dx*dx + dy*dy + dz*dz );
	
			if (ohdist < 1.15)
			{
				ohindices[j][count] = k;
				count ++;
			}
		}
		hcount[j] += count;
	}

	for (int j = 0; j < nooa; j ++)
	{
		int count = 0;

		for (int k = 0; k < nooa; k ++)
		{
			double odx = oxyz[j][0] - oxyz[k][0];
			double ody = oxyz[j][1] - oxyz[k][1];
			double odz = oxyz[j][2] - oxyz[k][2];

			odx -= lattice*pbc_round(odx/lattice);
                        ody -= lattice*pbc_round(ody/lattice);
                        odz -= lattice*pbc_round(odz/lattice);


			double oodist = sqrt (odx*odx + ody*ody + odz*odz );

			if (oodist > 0.0 && oodist < 3.6)
			{
				for (int n = 0; n < 4; n ++)
				{
					if ( n == 0)
					{
					if (ohindices[k][n] != -1)
					{
					double hdx = oxyz[j][0] - hxyz[ohindices[k][n]][0];
					double hdy = oxyz[j][1] - hxyz[ohindices[k][n]][1];
					double hdz = oxyz[j][2] - hxyz[ohindices[k][n]][2];
				
					hdx -= lattice*pbc_round(hdx/lattice);
                        		hdy -= lattice*pbc_round(hdy/lattice);
                        		hdz -= lattice*pbc_round(hdz/lattice);

					double hdist = sqrt( hdx*hdx + hdy*hdy + hdz*hdz );
					double dot = odx*hdx + ody*hdy + odz*hdz;
					double angle = acos (dot / (oodist*hdist)) * 57.2957795;
					if (angle < 30.0 && hdist < 2.4)
					{
						count ++;
					}

					}
					}
					if (n > 0)
                                        {
                                        if (ohindices[k][n] != -1 && ohindices[k][n] != ohindices[k][n-1])
                                        {
                                        double hdx = oxyz[j][0] - hxyz[ohindices[k][n]][0];
                                        double hdy = oxyz[j][1] - hxyz[ohindices[k][n]][1];
                                        double hdz = oxyz[j][2] - hxyz[ohindices[k][n]][2];

                                        hdx -= lattice*pbc_round(hdx/lattice);
                                        hdy -= lattice*pbc_round(hdy/lattice);
                                        hdz -= lattice*pbc_round(hdz/lattice);

                                        double hdist = sqrt( hdx*hdx + hdy*hdy + hdz*hdz );
                                        double dot = odx*hdx + ody*hdy + odz*hdz;
                                        double angle = acos (dot / (oodist*hdist)) * 57.2957795;
                                        if (angle < 30.0 && hdist < 2.4)
                                        {
                                                count ++;
                                        }

                                        }
                                        }

				}	
			}
		}
		hcount[j] += count;	
		xcoord = round(oxyz[j][0]);
		ycoord = round(oxyz[j][1]);
		zcoord = round(oxyz[j][2]);
		xyzbin[xcoord][ycoord][zcoord] ++;
		xyzhbin[xcoord][ycoord][zcoord] += hcount[j];
		
	}	


}
int linespace = 1;
for (int i = 0; i < 13; i ++)
{
	linespace = 0;
	for (int j = 0; j < 13; j ++)
	{	
		for (int k = 0; k < 13; k ++)
		{
			if (xyzbin[k][j][i] != 0)
			{
				hbonds_outputfile << xyzhbin[k][j][i]/xyzbin[k][j][i] << "  ";
			}
			else
			{
			hbonds_outputfile << xyzhbin[k][j][i] << "  ";
			}
			linespace ++;
			if (linespace % 5 == 0)
			{
				hbonds_outputfile << endl;
			}
		}
	}
}

inputfile.close();
hbonds_outputfile.close();

return 0;
} 

