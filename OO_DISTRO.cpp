/* This program is intended to take in an XYZ file, and output data for the O-O neighboring distance.
This can be used for metal-water interfaces as well.*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;


//function to find the minimum of the vector
double minval(vector <double> array, int array_length)
{
	double min = array[0];
	
	for (int i = 1; i < array_length; i ++)
	{
		if (array[i] < min)
		{
			min = array[i];
		}
	}
return min;
}

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
int main()
{

// this is the main menu
string infile;
double max, xlat, ylat, zlat;
int nooa, nbins;

cout << "XYZ file:\n==> ";
cin >> infile;
cout << "Unit cell dimensions (e.g. 10.0 10.0 10.0):\n==> ";
cin >> xlat >> ylat >> zlat;
cout << "Number of oxygen atoms:\n==> ";
cin >> nooa;
cout << "Number of bins:\n==> ";
cin >> nbins;
cout << "Max distance (Angstroms):\n==> ";
cin >> max;

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

cout << "\n Read input file. \n\n";

//number of frames
int nframes = ox.size()/nooa;

//go through the data and calculate some distances!

vector <double> final_distances;
ofstream output;
output.open("OO_DISTRO.dat");
double bin[nbins];
for (int i = 0; i < nbins; i ++)
{
	bin[i] = 0;
}
double incr = max / nbins;
int numsum = 0;
double distsum = 0.0;
for (int i = 0; i < nframes; i ++)
{
	vector <double> second_distances;

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
				dz -= zlat*pbc_round(dz/zlat);

				double dist = sqrt( dx*dx + dy*dy + dz*dz );

				distances.push_back(dist);
			
			}
		}

		second_distances.push_back(minval(distances, distances.size()));
	
	}
	
	
	
	for (int j = 0; j <  second_distances.size(); j ++)
	{
		for (int k = 0; k < second_distances.size(); k ++)
		{
			if (j != k && second_distances[j] == second_distances[k])
			{
				second_distances.erase(second_distances.begin()+k);
			}
		}
	}

	numsum += second_distances.size();

	for (int j = 0; j < second_distances.size(); j ++)
	{
		distsum += second_distances[j];
		int bin_num = second_distances[j]/incr;
		bin[bin_num] ++; 
	}
		
}
double average = distsum / numsum;

for (int i = 0; i < nbins; i ++)
{
	output << i*incr << "\t" << 100*bin[i]/(numsum) << "\t" << average << "\t" << 100*bin[i]/numsum << endl;
	output << (i+1)*incr << "\t" << 100*bin[i]/(numsum) << "\t" << average << "\t" << 100*bin[i]/numsum << endl;
} 

output.close();
input.close();
return 0;
}
				
	

	



