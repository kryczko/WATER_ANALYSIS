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
double xlat, ylat, zlat;

cout << "XYZ file:\n==>  ";
cin >> infile;
cout << "Number of atoms (e.g. 2 & 4 for 2 oxygen and 4 hydrogens):\n==> ";
cin >> nooa >> noha;
cout << "Unit cell dimensions (e.g. 10.0 10.0 10.0):\n==> ";
cin >> xlat >> ylat >> zlat;
cout << "\nPROGRAM RUNING, PLEASE WAIT.\n\n";
// end of main menu*/

//read the inputfile
inputfile.open(infile.c_str());

int noa = nooa + noha;

vector <double> ox, oy, oz, hx, hy, hz;
double xdat, ydat, zdat;
string data;

while (!inputfile.eof())
{
	inputfile >> data;
	if (data == "O")
	{
		inputfile >> xdat >> ydat >> zdat;
		ox.push_back(xdat);
		oy.push_back(ydat);
		oz.push_back(zdat);
	}
	if (data == "H")
	{
		inputfile >> xdat >> ydat >> zdat;
		ox.push_back(xdat);
		oy.push_back(ydat);
		oz.push_back(zdat);
	}

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
hbonds_outputfile.open("hbonds_histogram.dat");

int xcoord, ycoord, xybin[13] = {};
double xyhbin[13] = {}, hcount[nooa] ;

for (int i = 0; i < timesteps; i ++)
{
	for (int j = 0; j < nooa; j ++)
	{
		oxyz[j][0] = ox[j + i*nooa];
		oxyz[j][1] = oy[j + i*nooa];
		oxyz[j][2] = oz[j + i*nooa];
		hcount[j] = 0;
	}
	for (int j = 0; j < noha; j ++)
	{
		hxyz[j][0] = hx[j + i*noha];
		hxyz[j][1] = hy[j + i*noha];
		hxyz[j][2] = hz[j + i*noha];
	}
	for (int j = 0; j < nooa; j++)
	{
		int count = 0;

		for (int k = 0; k < noha; k ++)
		{
			double dx = oxyz[j][0] - hxyz[k][0];
			double dy = oxyz[j][1] - hxyz[k][1];
			double dz = oxyz[j][2] - hxyz[k][2];
	
			dx -= xlat*pbc_round(dx/xlat);
			dy -= ylat*pbc_round(dy/ylat);
			dz -= zlat*pbc_round(dz/zlat);			
	
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

			odx -= xlat*pbc_round(odx/xlat);
                        ody -= ylat*pbc_round(ody/ylat);
                        odz -= zlat*pbc_round(odz/zlat);


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
				
					hdx -= xlat*pbc_round(hdx/zlat);
                        		hdy -= ylat*pbc_round(hdy/ylat);
                        		hdz -= zlat*pbc_round(hdz/zlat);

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

                                        hdx -= xlat*pbc_round(hdx/xlat);
                                        hdy -= ylat*pbc_round(hdy/ylat);
                                        hdz -= zlat*pbc_round(hdz/zlat);

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
		xybin[xcoord] ++;
		xyhbin[xcoord] += hcount[j];
		
	}	


}

for (int i = 0; i < 13; i ++)
{

		if (xybin[i] != 0)
		{
		hbonds_outputfile << i << "\t" << xyhbin[i]/xybin[i] << endl;
                hbonds_outputfile << i + 1 << "\t" << xyhbin[i]/xybin[i] << endl;

		}
		else
		{
		hbonds_outputfile << i << "\t" << xyhbin[i] << endl;
                hbonds_outputfile << i + 1 << "\t" << xyhbin[i] << endl;
		
		}
	
}

inputfile.close();
hbonds_outputfile.close();

return 0;
} 

