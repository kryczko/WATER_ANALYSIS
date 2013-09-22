/* This program will plot the average dipole orientation of the water molecules wrt the z-axis.
   It is wrt the z-axis because this is perpendicular to the metal surface. To use this
   program you need an xyz file, the lattice parameters for your unit cell and you need
   to know how many oxygens/number of water molecules you have. This will give the program
   info about the number of frames.
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>
#include <vector>

using namespace std;

//This is the function to deal with periodic boundary conditions
//_______________________________________________________________
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
//_____________________________________________________________


// This is the main function
int main()
{
	//____________________________________________
	// This is the main menu
	string infile;
	double xlat, ylat, zlat;
	int nooa, noha, nbins;

	cout << "XYZ file:\n==> ";
	cin >> infile;
	cout << "Unit cell dimensions (x y z):\n==> ";
	cin >> xlat >> ylat >> zlat;
	cout << "Number of O atoms:\n==> ";
	cin >> nooa;
	cout << "Number of H atoms:\n==> ";
	cin >> noha;
	cout << "Number of bins:\n==> ";
	cin >> nbins;
	//____________________________________________

	//_________________________________________________________
	// Go through the input file and extract the data needed
	
	ifstream input;
	input.open(infile.c_str());
	
	string stuff;
	double x, y, z;
	vector <double> ox, oy, oz, hx, hy, hz;	

	while (! input.eof())
	{
		input >> stuff;
		
		if (stuff == "O")
		{
			input >> x >> y >> z;
			ox.push_back(x);
			oy.push_back(y);
			oz.push_back(z);
		}
		if (stuff == "H")
		{
			input >> x >> y >> z;
			hx.push_back(x);
			hy.push_back(y);
			hz.push_back(z);
		}
	}
	//_________________________________________________________
	
	int nframes = ox.size() / nooa;

	//_________________________________________________________________________
	//Firstly, create the neighbor list; water can have a max of 4 neighbors
	
	// initialize every element to -1 to use as a condition later
	double neighbor[nooa][4][4];
	for (int i = 0; i < nooa; i ++)
	{
		for (int j = 0; j < 4; j ++)
		{
			neighbor[i][j][0] = -1.0;
		}
	}
	//_________________________________________________________________________
	
	//________________________________________________________________________
	//now lets go through the data frame by frame, get the neighbors, compute 
	//the vectors and place them into an array/vector.
	
	vector <double> angles, zcoords;

	for (int i = 0; i < nframes; i ++)
	{
		for (int j = 0; j < nooa; j ++)
		{
			int hcount = 0;

			for (int k = 0; k < noha; k ++)
			{
				double dx = ox[j + i*nooa] - hx[k + i*noha];
				double dy = oy[j + i*nooa] - hy[k + i*noha];
				double dz = oz[j + i*nooa] - hz[k + i*noha];		
	
				dx -= xlat*pbc_round(dx/xlat);
				dy -= ylat*pbc_round(dy/ylat);
				dz -+ zlat*pbc_round(dz/zlat);
				
				double distance = sqrt( dx*dx + dy*dy + dz*dz );
				
				if (distance < 1.2)
				{
					neighbor[j][hcount][0] = (double) k;
					neighbor[j][hcount][1] = dx;
					neighbor[j][hcount][2] = dy;
					neighbor[j][hcount][3] = dz;
					hcount ++;
				}
			}
		}
		
		for (int j = 0; j < nooa; j ++)
		{
			double xsum = 0, ysum = 0, zsum = 0;	
		
			for (int k = 0; k < 4; k ++)
			{		
				if (neighbor[j][k][0] != -1.0)
				{
					xsum += neighbor[j][k][1];
					ysum += neighbor[j][k][2];
					zsum += neighbor[j][k][3];
					
				}
			}
			if (xsum != 0 && ysum !=0 && zsum !=0)
			{			
				double newvec = sqrt( xsum*xsum + ysum*ysum + zsum*zsum );
				// perform the dot product with z hat, the other terms become 0
				double angle = acos( (zsum)/(newvec) ) * 57.2957795;
				angles.push_back(angle);
				zcoords.push_back(oz[j + i*nooa]);
			}
		}
	}
	//_____________________________________________________________________________________
			
	//_________________________________________________________________________
	//output data to file

	ofstream output;
	output.open("dipoleangles.dat");
	double binsize = zlat/nbins, anglesum[nbins];
	int zbins[nbins];
	
	for ( int i = 0; i < nbins; i ++)
	{
		zbins[i] = 0;
		anglesum[i] = 0;
	}	

	for ( int i = 0; i < zcoords.size(); i ++)
	{
		int bin_num = zcoords[i]/binsize;
		zbins[bin_num] ++;
		anglesum[bin_num] += angles[i];
	}
	for (int i = 0; i < nbins; i ++)
	{
		if (zbins[i] != 0)
		{
			output << i*binsize << "\t" << anglesum[i]/zbins[i] << endl;
			output << (i+1)*binsize << "\t" << anglesum[i]/zbins[i] << endl;
		}
		else
		{
			output << i*binsize << "\t" << -1.0 << endl;
                        output << (i+1)*binsize << "\t" << -1.0 << endl;
		}
	}	
	//________________________________________________________________________
	input.close();
	output.close();
	return 0;
}				
	






	
