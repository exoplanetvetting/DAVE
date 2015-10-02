// Written by J. Coughlin (jlcough@nmsu.edu). Last modified 9/20/2011 by J. Coughlin.
//
// Bin a file in that has 2 columns of Xvalue  Yvalue
// Uses error-weighted mean and constant bin steps
// Outputs new binned Xvalue, Yvalue, and error (from stdev) of the binned data point
//
// To compile: g++ -O1 -o bin bin.cpp
// To run: ./bin INFILE BINSIZE OUTFILE
// Example: ./bin mydata.dat 0.05 mybinneddata.dat


#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>

#define N 100000  // Increase to more than 100,000 if you have more than 100,000 data points

using namespace std;

double MEAN(int,double[],double[]),MEANERR(int,double[],double[]),STDEV(double[],int);

int i,j,k;  // Counting integers
int n,sw1;  // n is number of points going into current bin - sw1 is a switch
int Nraw,Nbin; // Number of original and binned data points
double binsize,curbin;
double curdat2[N],curdat3[N];  // Temp arrays for calculations
double inpdata[2][N];  // Input data
double bindat[3][N];   // Output binned data
ifstream datain;
ofstream dataout;
string infile,outfile;

int main (int argc, char* argv[]) {

if(argc>1)
  infile = argv[1];
else
  {
  cout << "Name of input file? ";
  cin >> infile;
  }

if(argc>2)
  binsize = atof(argv[2]);
else
  {
  cout << "Enter binsize in xaxis units: ";
  cin >> binsize;
  }
while(binsize<0)
  {
  cout << endl << "Not a valid bin size. Enter binsize (or 0 to replot original data): ";
  cin >> binsize;
  }
	
if(argc>3)
  outfile = argv[3];
else
  {
  cout << "Name of output file? ";
  cin >> outfile;
  }

if(outfile==infile)
  {
  cout << "Output file cannot be the same as input file. Exiting." << endl;
  exit(0);
  }

datain.open(infile.c_str());
i=0;
datain >> inpdata[0][i];
while(!datain.eof())
  {
  datain >> inpdata[1][i];
  i=i+1;
  datain >> inpdata[0][i];
  }
Nraw = i;
datain.close();

Nbin=0;
for(curbin=inpdata[0][0];curbin<inpdata[0][Nraw-1]+binsize;curbin+=binsize)
  {
  n=0;
  sw1=0;
  for(i=0;i<Nraw;i++)
    if(inpdata[0][i]>=curbin && inpdata[0][i]<curbin+binsize)
      {
      curdat2[n] = inpdata[1][i];
      curdat3[n] = 1.0;
      n++;
      sw1=1;
      }
  if(sw1==1)
    {
    bindat[0][Nbin] = curbin+0.5*binsize;
    bindat[1][Nbin] = MEAN(n,curdat2,curdat3);
    bindat[2][Nbin] = STDEV(curdat2,n)/sqrt(n);
    Nbin++;
    }
  }
  
dataout.open(outfile.c_str());
for(i=0;i<Nbin;i++)
  dataout << setprecision(10) << bindat[0][i] << " " << bindat[1][i] << " " << bindat[2][i] << endl;


}

////////////////////////////////////////////////////////////////////////////////////////////////////

double MEAN(int n, double x[], double xerr[]) {  // Return error-weighted mean of first n terms

int i,j,k;
double invvarsum,xmean;

invvarsum=0.0;
for(i=0;i<n;i++)
  invvarsum += 1.0/(xerr[i]*xerr[i]);

xmean=0;
for(i=0;i<n;i++)
  xmean += x[i]/(xerr[i]*xerr[i]);

return(xmean/invvarsum);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

double MEANERR(int n, double x[], double xerr[]) {  // Return the error of the error-weighted mean of first n terms

int i,j,k;
double invvarsum;

invvarsum=0.0;
for(i=0;i<n;i++)
  invvarsum += 1.0/(xerr[i]*xerr[i]);

return(sqrt(1.0/invvarsum));
}

////////////////////////////////////////////////////////////////////////////////////////////////////

double STDEV(double y[],int a) {
/*-------------------------------------------------------------
  Calculate the standard deviation of an array y, given a terms
  -------------------------------------------------------------*/

double ymean,sum;
int z;

ymean=0;
for(z=0;z<a;z++)
  ymean+=y[z];
ymean/=a;

sum=0;
for(z=0;z<a;z++)
  sum+=pow((y[z]-ymean),2);


return sqrt(sum/a);
}