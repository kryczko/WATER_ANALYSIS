/*This analysis code is meant for a single metal-water interface. YOU NEED AN XYZ FILE.*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

using namespace std;

int main()
{
	string infile, atomname, convert;
	double xlat, ylat, zlat;
	int zbins, xybins, natoms;

	cout << "XYZ file:\n==> ";
	cin >> infile;
	cout << "Lattice constants (x y z) in Angstroms:\n==> ";
	cin >> xlat >> ylat >> zlat;
	cout << "Number of O atoms:\n==> ";
	cin >> natoms;
	cout << "Number of bins (zbins x/ybins):\n==> ";
	cin >> zbins >> xybins;
	cout << "Heavy water? (y or n):\n==> ";
	cin >> convert;

	ifstream input;
	input.open(infile.c_str());
	
	string content;
	vector <double> x, y, z;
	double xval, yval, zval;

	// Go through the input file and extract the oxygen atoms
	//################################################
	while (!input.eof())
	{
		input >> content;
		if (content == "O")
		{
			input >> xval >> yval >> zval;
			x.push_back(xval);
			y.push_back(yval);
			z.push_back(zval);
		} 
	}
	int nframes = x.size() / natoms;
	//##############################################
	
	// Wrap the coordinates
	//#########################################
	for (int i = 0; i < x.size(); i ++)
	{
		if (x[i] > xlat)
		{
			x[i] -= xlat;
		}
		if (x[i] < 0.0)
		{
			x[i] += xlat;
		}
		if (y[i] > ylat)
                {
                        y[i] -= ylat;
                }
                if (y[i] < 0.0)
                {
                        y[i] += ylat;
                }
		if (z[i] > zlat)
                {
                        z[i] -= zlat;
                }
                if (z[i] < 0.0)
                {
                        z[i] += zlat;
                }
	}
	//#########################################
	
	// declare the bins; set the values to 0
	//###############################################################################

	double xzbins[xybins][zbins], yzbins[xybins][zbins]; ;
	for (int i = 0; i < xybins; i ++)
	{
		for (int j = 0; j < zbins; j ++)
		{
			xzbins[i][j] = 0;
			yzbins[i][j] = 0;
		}
	}
	//################################################################################

	// write the bins
	//###################################
	double xbinsize = xlat/xybins, ybinsize = ylat/xybins, zbinsize = zlat/zbins;

	for (int i = 0; i < x.size(); i ++)
	{
		int xbin_num = x[i]/xbinsize;
		int zbin_num = z[i]/zbinsize;
		int ybin_num = y[i]/ybinsize;
		xzbins[xbin_num][zbin_num] ++;
		yzbins[ybin_num][zbin_num] ++;
	}
	//#####################################

	// write out the data to the designated data files so it can be plotted
	//##################################################################
	ofstream xzdens, yzdens;
	xzdens.open("xzdensity.dat");
	yzdens.open("yzdensity.dat");

	double conversion = 18.0e-6/(6.023e23*1.0e-30); // go to from mol/A^3 to g/cc
	if (convert == "y")
	{
		conversion = 20.0e-6/(6.023e23*1.0e-30); // go to from mol/A^3 to g/cc;
	}	

	double xinc = xlat/xybins, yinc = ylat/xybins, zinc = zlat/zbins;
	double xzvolume = xinc*ylat*zinc, yzvolume =xlat*yinc*zinc;
	
	int xzcount(0), yzcount(0);
	double xzsum=0, yzsum=0, xzmean, yzmean;
	for (int i = 0; i < xybins; i ++)
	{
		for (int j = 0; j < zbins; j ++)
		{
			if (xzbins[i][j]*conversion/(xzvolume*nframes) != 0)
			{
				xzsum += xzbins[i][j]*conversion/(xzvolume*nframes);
				xzcount ++;
			}
			if ((yzbins[i][j]*conversion)/(yzvolume*nframes) != 0)
			{
				yzsum += (yzbins[i][j]*conversion)/(yzvolume*nframes);
				yzcount ++;
			}
		}
	}
	xzmean = xzsum / xzcount;
	yzmean = yzsum / yzcount;
	

	for (int i = 0; i < xybins; i ++)
	{
		for (int j = 0; j < zbins; j ++)
		{
			xzdens << i*xbinsize << "\t" << j*zbinsize << "\t" <<  xzbins[i][j]*conversion/(xzvolume*nframes) << "\t" << xzmean << endl;
		}
	}
	for (int i = 0; i < xybins; i ++)
        {
		for (int j = 0; j < zbins; j ++)
		{
        	        yzdens << i*ybinsize << "\t" << j*zbinsize << "\t" << (yzbins[i][j]*conversion)/(yzvolume*nframes) << "\t" << yzmean << endl;
        	}
	}
	//###################################################################
	
	xzdens.close();
	yzdens.close();
	input.close();
	return 0;
}		
			
		
