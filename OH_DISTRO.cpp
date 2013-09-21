/*This program prompts the user for an xyz file, and needs some basic info
like how many atoms, unit cell etc. The output can be plotted using software like gnuplot*/

#include <iostream> // input & output
#include <fstream> // input & output files
#include <cmath> // basic math functions
#include <string> // for string capabilities
#include <vector> // use of dynamic arrays

using namespace std; //include all above libraries

// this is the function that deals with periodic boundary conditions
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

int main() // this is the main function
{
// this is the main menu
string infile;
int nooa, noha, nbins;
double xlat, ylat, zlat, max;

cout << "XYZ file:\n==> ";
cin >> infile;
cout << "How many atoms (e.g. 2 4 for 2 oxygen and 4 hydrogen):\n==> ";
cin >> nooa >> noha;
cout << "Unit cell dimensions (e.g. 10.0 10.0 10.0 for x,y,z):\n==> ";
cin >> xlat >> ylat >> zlat;
cout << "Number of bins:\n==> ";
cin >> nbins;
cout << "Maximum distance:\n==> ";
cin >> max;
cout << "\nCODE RUNNING, PLEASE WAIT.\n";

// open the file declared above
ifstream input;
input.open(infile.c_str());

// go through the input file 
string stuff;
vector <double> ox, oy, oz, hx, hy, hz;

while(! input.eof())
{
	input >> stuff;
	
	if(stuff == "O")
	{
		input >> ox >> oy >> oz;
	}
	if(stuff == "H")
	{
		input >> hx >> hy >> hz;
	}
}
cout << "\nREAD INPUTFILE...\n";

// declare the number of frames
int nframes = ox.size() / nooa;

// now the fun part. first we create a neighbor list, then we use that neighbor list to
// find the distances we need, then we put each distance into each bin and output the data.

int nlist1[nooa][4]; // used for indexing 
double nlist2[nooa][4]; // used to keep the distances with the index

// initialize the neighborlists so that they have -1 if no hydrogen is there
for (int i = 0; i < nooa; i ++)
{
	for (int j = 0; j < 4; j ++)
	{
		nlist1[i][j] = -1;
		nlist2[i][j] = -1;
	}
}

double binsize = max/nbins;
int bin[nbins];
for (int i = 0; i < nbins; i ++)
{
	bin[i] = 0;
}

for (int i = 0; i < nframes; i ++)
{
	for (int j = 0; j < nooa; j ++)
	{
		int count = 0;
		for (int k = 0; k < noha; k ++)
		{
			double dx = ox[j + i*nooa] - hx[k + i*noha];
			double dy = oy[j + i*nooa] - hy[k + i*noha];
			double dz = oz[j + i*nooa] - hz[k + i*noha];	

			dx -= xlat*pbc_round(dx/xlat);
			dy -= ylat*pbc_round(dy/ylat);
			dz -= zlat*pbc_round(dz/zlat);
			
			double dist = sqrt(dx*dx + dy*dy + dz*dz);

			if (dist < 1.2) // less than 1.2 angstroms, estimated distance
			{
				nlist1[j][count] = k;			
				nlist2[j][count] = dist;

				count ++;
			}
		}
	}
	
	for (int j = 0; j < nooa; j ++)
	{
		for (int k = 0; k < 4; k ++)
		{
			if (nlist1[j][k] != -1)
			{
				int bin_number = nlist2[j][k]/binsize;
				bin[bin_number] ++;
			}
		}
	}							
}



ofstream output;
output.open("OHD.dat");	

for (int i = 0; i < nbins; i ++)
{
	output << i << "\t" << bin[i] << endl;
	output << i + 1 << "\t" << bin[i] << endl;
}

input.close();
output.close();
return 0;

}





























