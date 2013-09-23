#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main()
{

int timesteps, number_of_atoms;
double lattice;
string infile;

cout << "Please enter the filename of your file: ";
cin >> infile;
cout << "Please enter the number of atoms: ";
cin >> number_of_atoms;
cout << "Please enter the number of timesteps: ";
cin >> timesteps;
cout << "Please enter the lattice constant for your periodic cube: ";
cin >> lattice;
cout << "Program running...please wait a moment.\n\n";

ifstream input;
ofstream output;

input.open(infile.c_str());
output.open("coords.xyz");


double x[number_of_atoms*timesteps], y[number_of_atoms*timesteps], z[number_of_atoms*timesteps];
int counter = 0;

while (!input.eof())
{
	input >> x[counter] >> y[counter] >> z[counter];	
	counter ++;

}

	int nooa = number_of_atoms/3, noha = 2*number_of_atoms/3;

	double oxyz[nooa*timesteps][3], hxyz[noha*timesteps][3];

	for (int i = 0; i < timesteps; i ++)
	{
		for (int j = 0; j < nooa; j ++)
		{
			oxyz[j + i*nooa][0] = lattice*x[j + i*number_of_atoms];
                        oxyz[j + i*nooa][1] = lattice*y[j + i*number_of_atoms];
                        oxyz[j + i*nooa][2] = lattice*z[j + i*number_of_atoms];
                }
	}
	
	   for (int i = 0; i < timesteps; i ++)
        {
                for (int j = 0; j < noha; j ++)
                {
                        hxyz[j + i*noha][0] = lattice*x[nooa + j + i*number_of_atoms];
                        hxyz[j + i*noha][1] = lattice*y[nooa + j + i*number_of_atoms];
                        hxyz[j + i*noha][2] = lattice*z[nooa + j + i*number_of_atoms];
                }
        }
int count = 1;
	for (int i = 0; i < timesteps*nooa; i ++)
	{
		if (i % nooa == 0)
		{
			output << number_of_atoms << "\n";
			output << count << endl;
			count ++;
		}
		output << "O" << "  " << oxyz[i][0] << "  " << oxyz[i][1] << "  " << oxyz[i][2] << endl;
		output << "H" << "  " << hxyz[2*i][0] << "  " << hxyz[2*i][1] << "  " << hxyz[2*i][2] << endl;
		output << "H" << "  " << hxyz[2*i + 1][0] << "  " << hxyz[2*i + 1][1] << "  " << hxyz[2*i + 1][2] << endl;		
	}
output.close();
input.close();

return 0;
}


















