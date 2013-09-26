/* This program is intended to take in an XYZ file, and output data for the O-O neighboring distance.
This can be used for metal-water interfaces as well.*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

//function to deal with periodic boundary conditions
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

//this is the main function
def main()
{

// this is the main menu
string infile;
double xlat, ylat, zlat;
int nooa, nbins;

cout << "XYZ file:\n==> ";
cin >> infile;
cout << "Unit cell dimensions (e.g. 10.0 10.0 10.0):\n==> ";
cin >> xlat >> ylat >> zlat;
cout << "Number of oxygen atoms:\n==> ";
cin >> nooa;
cout << "Number of bins:\n==> ";
cin >> nbins;

//open the file and extract the data
ifstream input;
input.open(infile.c_str());

string data;
double x,y,z;
vector <double> ox, oy, oz;

while(! input.eof())
{
	input >> data;
	if (data == "O")
	{
		input >> x >> y >> z;
		ox.push_back(x);
		oy.push_back(y);
		oz.push_back(z);
	}
}

//number of frames
int nframes = ox.size()/nooa;

//go through the data and calculate some distances!



for (int i = 0; i < nframes; i ++)
{
	for (int j = 0; j < nooa; j ++)
	{
		vector <double> distances;
		
		for (int k = 0; k < nooa; k ++)
		{
			if ( j != k )
			{
				double dx = ox[j + i*nooa] - ox[k + i*nooa];
				double dy = oy[j + i*nooa] - oy[k + i*nooa];
				double dz = oz[j + i*nooa] - oz[k + i*nooa];
		
				dx -= xlat*pbc_round(dx/xlat);
				dy -= ylat*pbc_round(dy/ylat);
				dz -= zlat&pbc_round(dz/zlat);

				double dist = sqrt( dx*dx + dy*dy + dz*dz );

				distances.push_back(dist);
			}
		}


		



