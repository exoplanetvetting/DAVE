/* Jeff Coughl n
   modshift.cpp

TO-DO: MAKE SURE I AM USING STDEV for thefinal rmsmin calc AND NOT RMS FOR NOISE CALCS. Just check it out if there's a big diff. Checked - nope no big diff. Data is already mean subtracted.
TO-DO: 

   awk '{print("echo -n "substr($1,1,12)"\ >> modshift.txt\n../modshift "$1" >> modshift.txt\nmv ModShift.ps "substr($1,1,12)"-modshift.ps")}' alldatfiles > runmodshiftall
     
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <time.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <map>
#include <utility>
#include <string.h>
#include <math.h>
#include <algorithm>
#include <sstream>


using namespace std;


// Declare constants
const int N=100000;  // Max number of data points
const double ONESIGCONF = 0.682689492137;  // One-sigma confidence interval
const double TWOSIGCONF = 0.954499736104;  // Two-sigma confidence interval
const double THREESIGCONF = 0.997300203936740; // Three-sigma confidence interval
// const double EFAC = 2.71828182846;


// Declare Functions
double STDEV(double[],int),RMS(double[],int),MEDIAN(double[],int),SUM(double[],int),SUMSQ(double[],int),SELECT(double[],int,int),INVERFC(double),MAD(double[],int);

// Declare Global - Only using global because of memory allocation problems. If I declare local i get seg faults.
double data[4][2*N];  // 0 = Time; 1 = Original Data; 2 = Original Model; 3 = Original Residuals
double rms[N],sumres2[N],zparm[N],tmpdob1[N],tmpdob2[N],chi2[N],tmpsum[N]; //,runmean[N],depthfac[N];


// Function for sorting via qsort
int cmp(const void *x, const void *y)
  {
  double xx = *(double*)x, yy = *(double*)y;
  if (xx < yy) return -1;
  if (xx > yy) return  1;
  return 0;
  }  
  

int main (int argc, char* argv[]) {

int i,j,k,l,m,n;
ifstream infile;
ofstream outfile;
string basename,infilename;
double rmsmin=9E9,prilowtime,seclowtime,sechightime,terlowtime,chi2med,chi2outlow,chi2outhigh,chi2inlow,chi2low,chi2terlow;
int ndat,sw1,sw2,ntmp;
double tstart,tend;
int midti,startti,endti;  // Index for middle of transit point, start of transit, and end of transit
double sigpri,sigsec,sigter,sigpos;
double period;
string syscmd;
stringstream convert;
int widthfac;
double width,halfwidth,tenthwidth;
double tdepth,nintrans,depfacsec,depfacter,depfacpos,depsig,sigfa,fred;
double med,std,sigreject,mad;
int tmpint1[N];
string tmpstr1;
double baseflux;
double flat[2*N]; // residual assuming flat baseline
double tmpsum1,tmpsum2;

// clock_t startTime = clock();


if(argc>1)
  infilename = argv[1];
else
  {
  cout << "Name of file with data and model fit? ";
  cin >> infilename;
  }

if(argc>2)
  basename = argv[2];
else
  {
  cout << "Name of basename? ";
  cin >> basename;
  }  

if(argc>3)
  period = atof(argv[3]);
else
  {
  cout << "Period? ";
  cin >> period;
  }


infile.open(infilename.c_str());
if(!infile.good())
  {
  cout << "Input file does not exist or cannot read! Exiting..." << endl;
  exit(0);
  }
sw1=0;
sw2=0;
i=0;
while(!infile.eof())
  {
  for(j=0;j<=2;j++)
    infile >> data[j][i];
  
  if(i==0)
    baseflux = data[2][i];
  
  if(sw1==0 && data[2][i]<baseflux)
    {
    tstart=data[0][i];
    sw1=1;
    }
  if(sw1==1 && data[2][i]==baseflux)
    {
    tend=data[0][i-1];
    sw1=2;
    }
  if(sw2==0 && data[0][i] > 0)
    {
    midti = i;
    sw2=1;
    }
  i++;
  }
infile.close();
ndat = i-1;



// Check to make sure model isn't all zeros, or only positive. Should be transit-like.
j=0;
for(i=0;i<ndat;i++)
  if(data[2][i]<0.0)
    j=1;
if(j==0)
  {
  cout << "Model is all zeros! Exiting..." << endl;
  exit(0);
  }



// Record model depth of primary transit
tdepth = data[2][midti];


if(tstart<0)  // Make width symmetrical. Also prevents glitches on cases if no in-transit data at positive phase.
  tend = fabs(tstart); 
else  // If tstart is positive then data or fit is not right, or no in-transit data points before phase 0.0, so go off the end time
  tstart = -1.0*tend;


// if(tstart==0.0)  // If tstart is 0.0, model is flat, or somethign went wrong, so revert to default values
//   {
//   tstart = -0.1;
//   tend = 0.1;
//   }
//   
// Remove very large outliers
// cout << ndat << endl;
//  cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC << " seconds." << endl;


width = tend-tstart;  // Record Width
tenthwidth = 0.1*width;
for(l=0;l<1;l++)  // Number of outlier passes - ONLY NEED ONE WHEN USING MAD - JUST LEAVING LOOP IN CASE I EVER WANT TO CHANGE MY MIND
  {
  for(i=0;i<ndat;i++)
    data[3][i] = data[1][i] - data[2][i];  // Calculate residuals

//   std = STDEV(data[3],ndat);
  mad = MAD(data[3],ndat);  // MAD is Median Absolute Devitation - better than STD - only need one pass.
  std = 1.4826*mad;  // https://www.mathworks.com/help/stats/mad.html  AND https://en.wikipedia.org/wiki/Median_absolute_deviation
    
  sigreject = sqrt(2)*INVERFC(1.0/ndat);  // Calculate sigma threshold based on number of data points
    
  for(i=0;i<ndat;i++)  // Keep track of which data points to remove. 0 is good, if flagged as 1, remove later
    tmpint1[i]=0;

  for(i=0;i<ndat;i++)
    {
    k=0;
    sw1=0;
    for(j=0;j<ndat;j++)
      if(fabs(data[0][i]-data[0][j]) < tenthwidth)  // look at points within a transit duration
        {
        tmpdob1[k]=data[3][j];
        k++;
        }
    med = MEDIAN(tmpdob1,k);
    
    if(fabs(data[3][i] - med) > sigreject*std)  // Keep data points within sigreject standard deviations of median
      tmpint1[i] = 1;
    }

  m=0;
  for(i=0;i<ndat;i++)
    if(tmpint1[i]==0)  // Keep good data points
      {
      for(j=0;j<4;j++)
        data[j][m] = data[j][i];
      m++;
      }
    else
      if(m<midti)
        midti--;  // If I remove a data point, I have to update midti to keep track of it
  ndat=m;  // Update number of data points to those not thrown out
  }

// cout << ndat << endl;
  
// Double up input data for shifting
for(i=0;i<ndat;i++)
  for(j=0;j<4;j++)
    data[j][i+ndat] = data[j][i];


// cout << data[0][63150] << " " << data[1][63150] << " " << data[2][63150] << endl;

  
// Do the convolution

////////////////////////////////////////////////////////////////////
// Save old in case I need to revert
// To speed things up, assume outside of transit is flat
// for(j=0;j<2*ndat;j++)
//   flat[j] = data[1][j] - baseflux;
// 
// 
// 
// for(i=0;i<ndat;i++)  // Perform ndat pertubations
//   {
//   for(j=0;j<ndat;j++)
//     {
//     if(fabs(data[0][j]) > 0.5*width)  // If data point is outside transit, can use the pre-computed value that assumes the baseline is flat.
//       data[3][j] = flat[j+i];
//     else
//       data[3][j] = data[1][j+i] - data[2][j];  // Shitfing data, holding model steady. Moving data points backwards, or model forwards, same thing
//     }
// 
//   rms[i] = RMS(data[3],ndat);  // RMS of the new residuals
//   if(rms[i] < rmsmin)
//     rmsmin = rms[i];
//   }
//   cout << rms[0] << " " << rms[1] << endl;
//  cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC << " seconds." << endl;
////////////////////////////////////////////////////////////////////



halfwidth = 0.5*width;  // Pre-computing so don't have to do inside loop

// Record the indexs of points that correspond to start and end of trasnit.
sw1=0;
for(j=0;j<ndat;j++)
  { 
  if(sw1==0 && data[0][j] > -halfwidth)
    {
    startti = j-1;  // Cadence of the last point before the transit occurs
    sw1=1;
    }
  if(sw1==1 && fabs(data[0][j])>halfwidth)
    {
    endti = j;  // Cadence of the point just after the transit ends.
    sw1=2;
    }
  }


// To speed things up, assume outside of transit is flat
for(j=0;j<2*ndat;j++)
  flat[j] = pow(data[1][j] - baseflux,2);

// Old code next 8 lines - was trying to speed things up, but never quite figured it out
// for(i=0;i<ndat;i++)
//   {
//   tmpsum[i] = 0;
// //   for(j=0;j<startti;j++)
// //     tmpsum[i] += flat[j+i];
// //   for(j=endti;j<ndat;j++)
// //     tmpsum[i] += flat[j+i];
//   }

// Okay do the actual perumtation. 
for(i=0;i<ndat;i++)  // Perform ndat pertubations
  {
  tmpsum2 = 0;
//   tmpsum2=tmpsum[i];  // Start summation at 0 - OLD attempt at speedup
  
  // Before transit, can look up values for compuatation speed increase
  for(j=0;j<startti;j++)
    tmpsum2 += flat[j+i];
  
  // Compute new values inside transit
  for(j=startti;j<endti;j++)
    tmpsum2 += pow(data[1][j+i] - data[2][j],2);  // Shitfing data, holding model steady. Moving data points backwards, or model forwards, same thing
    
  // After transit, can look up values for computation speed increase
  for(j=endti;j<ndat;j++)
    tmpsum2 += flat[j+i];

  rms[i] = sqrt(tmpsum2/ndat);  //RMS(data[3],ndat);  // RMS of the new residuals
  if(rms[i] < rmsmin)
    rmsmin = rms[i];
  }
  
  
 
 
// Now look at convolved ata to find pri, sec, etc. and do other things
nintrans=0; 
chi2low = 9E9;
tmpstr1 = "outfile1-" + basename + ".dat";
outfile.open(tmpstr1.c_str());
for(i=0;i<ndat;i++)  // Calculate chi squared and find lowest chi squared value
  {
  chi2[i] = ndat*pow(rms[i],2)/pow(rmsmin,2);  // Using rms to calculate chi squared
  if(chi2[i] < chi2low)
    chi2low = chi2[i];
  outfile << fixed << setprecision(10) << data[0][i]/period << " " << setprecision(8) << setw(11) << data[1][i] << " " << data[2][i] << endl; // Write out original data while we're at it
  if(data[0][i+midti] > tstart && data[0][i+midti] < tend)  // If within the actual transit window, compute number of in-transit data points while we're at it
    nintrans++;
  }
outfile.close();


// Search for sec. eclipse and compute in-primary params
chi2outlow = chi2inlow = 9E9;
seclowtime = -2*period;
prilowtime = 0;
widthfac=4;
for(i=0;i<ndat;i++)
  {
  if((data[0][i+midti] < widthfac*tstart || data[0][i+midti] > widthfac*tend) && (data[0][i+midti] < widthfac*tstart+period || data[0][i+midti] > widthfac*tend+period) && (data[0][i+midti] < widthfac*tstart-period || data[0][i+midti] > widthfac*tend-period))  // If outside primary transit
    {
    if(chi2[i]<chi2outlow)
      {
      chi2outlow=chi2[i];
      seclowtime = data[0][i+midti];
      }
    }
  else
    {
    if(chi2[i]<chi2inlow)
      {
      chi2inlow=chi2[i];
      prilowtime = data[0][i+midti];
      }
    }
  }
  
// Find tertiary and biggest positive peak that's outside primary and secondary eclipse
chi2terlow = 9E9;
terlowtime = -2*period;
widthfac=4;
for(i=0;i<ndat;i++)
  {
  if((data[0][i+midti] < prilowtime+widthfac*tstart || data[0][i+midti] > prilowtime+widthfac*tend) && (data[0][i+midti] < prilowtime+widthfac*tstart+period || data[0][i+midti] > prilowtime+widthfac*tend+period) && (data[0][i+midti] < prilowtime+widthfac*tstart-period || data[0][i+midti] > prilowtime+widthfac*tend-period) && (data[0][i+midti] < seclowtime+widthfac*tstart || data[0][i+midti] > seclowtime+widthfac*tend) && (data[0][i+midti] < seclowtime+widthfac*tstart+period || data[0][i+midti] > seclowtime+widthfac*tend)+period && (data[0][i+midti] < seclowtime+widthfac*tstart-period || data[0][i+midti] > seclowtime+widthfac*tend-period)  )  // If outside primray and secondary     
    if(chi2[i]<chi2terlow)
      {
      chi2terlow=chi2[i];
      terlowtime=data[0][i+midti];
      }
  }

// Find biggest positive peak that's outside primary and secondary eclipse
chi2outhigh = -9E9;
sechightime = -2*period;
widthfac=6;
for(i=0;i<ndat;i++)
  {
  if((data[0][i+midti] < prilowtime+widthfac*tstart || data[0][i+midti] > prilowtime+widthfac*tend) && (data[0][i+midti] < prilowtime+widthfac*tstart+period || data[0][i+midti] > prilowtime+widthfac*tend+period) && (data[0][i+midti] < prilowtime+widthfac*tstart-period || data[0][i+midti] > prilowtime+widthfac*tend-period) && (data[0][i+midti] < seclowtime+widthfac*tstart || data[0][i+midti] > seclowtime+widthfac*tend) && (data[0][i+midti] < seclowtime+widthfac*tstart+period || data[0][i+midti] > seclowtime+widthfac*tend)+period && (data[0][i+midti] < seclowtime+widthfac*tstart-period || data[0][i+midti] > seclowtime+widthfac*tend-period)  )  // If outside primray and secondary     
    if(chi2[i]>chi2outhigh)
      {
      chi2outhigh=chi2[i];
      sechightime = data[0][i+midti];
      }
  }

// Compute median excluding primary and secondary eclipse. I have built in contigency in cases of sparse data so this shoudl never fail.
ntmp=0;
widthfac=2;
for(i=0;i<ndat;i++)
  if((data[0][i+midti] < prilowtime+widthfac*tstart || data[0][i+midti] > prilowtime+widthfac*tend) && (data[0][i+midti] < prilowtime+widthfac*tstart+period || data[0][i+midti] > prilowtime+widthfac*tend+period) && (data[0][i+midti] < prilowtime+widthfac*tstart-period || data[0][i+midti] > prilowtime+widthfac*tend-period) && (data[0][i+midti] < seclowtime+widthfac*tstart || data[0][i+midti] > seclowtime+widthfac*tend) && (data[0][i+midti] < seclowtime+widthfac*tstart+period || data[0][i+midti] > seclowtime+widthfac*tend)+period && (data[0][i+midti] < seclowtime+widthfac*tstart-period || data[0][i+midti] > seclowtime+widthfac*tend-period)  )  // If outside primray and secondary     
    {
    tmpdob1[ntmp]=chi2[i];
    ntmp++;
    }
if(ntmp>0)
  chi2med = MEDIAN(tmpdob1,ntmp);
else // First contigency - if there were no points outside twice the width of the pri and secondary eclipse, try again by only excluding data inside primary
  {
  ntmp=0;
  widthfac=2;
  for(i=0;i<ndat;i++)
    if((data[0][i+midti] < widthfac*tstart || data[0][i+midti] > widthfac*tend) && (data[0][i+midti] < widthfac*tstart+period || data[0][i+midti] > widthfac*tend+period) && (data[0][i+midti] < widthfac*tstart-period || data[0][i+midti] > widthfac*tend-period))  // If outside primary transit
      {
      tmpdob1[ntmp]=chi2[i];
      ntmp++;
      }
  if(ntmp>0)
    chi2med = MEDIAN(tmpdob1,ntmp);
  else  // Second contigency - if still can't find any data just outside primary, try it by allowing for just a widthfac of 1 instead of 2, i.e., literally anything outside transit
    {
    ntmp=0;
    widthfac=1;
    for(i=0;i<ndat;i++)
      if((data[0][i+midti] < widthfac*tstart || data[0][i+midti] > widthfac*tend) && (data[0][i+midti] < widthfac*tstart+period || data[0][i+midti] > widthfac*tend+period) && (data[0][i+midti] < widthfac*tstart-period || data[0][i+midti] > widthfac*tend-period))  // If outside primary transit
        {
        tmpdob1[ntmp]=chi2[i];
        ntmp++;
        }
    if(ntmp>0)
      chi2med = MEDIAN(tmpdob1,ntmp);
    else  // Really? So the only data points we have are in-transit? Allright just use those...the only way we could fail now is if there is literally no data, in which case we should have crashed long ago.
      {
      ntmp=0;
      for(i=0;i<ndat;i++)
        {
        tmpdob1[ntmp]=chi2[i];
        ntmp++;
        }
      if(ntmp>0)
        chi2med = MEDIAN(tmpdob1,ntmp);
      else
        {
        cout << "You should never see this message. I think you have no data at all in your light curve." << endl;
        exit(0);
        }
      }
    }
  }

// Calculate rmsmin and Fred, ratio of gaussian noise to systematic noise, excluding primary and secondary
ntmp=0;
widthfac=2;
for(i=0;i<ndat;i++)
  if((data[0][i+midti] < prilowtime+widthfac*tstart || data[0][i+midti] > prilowtime+widthfac*tend) && (data[0][i+midti] < prilowtime+widthfac*tstart+period || data[0][i+midti] > prilowtime+widthfac*tend+period) && (data[0][i+midti] < prilowtime+widthfac*tstart-period || data[0][i+midti] > prilowtime+widthfac*tend-period) && (data[0][i+midti] < seclowtime+widthfac*tstart || data[0][i+midti] > seclowtime+widthfac*tend) && (data[0][i+midti] < seclowtime+widthfac*tstart+period || data[0][i+midti] > seclowtime+widthfac*tend)+period && (data[0][i+midti] < seclowtime+widthfac*tstart-period || data[0][i+midti] > seclowtime+widthfac*tend-period)  )  // If outside primray and secondary     
    {
    tmpdob1[ntmp]=data[1][i+midti]-data[2][i+midti];  // Original data
    tmpdob2[ntmp]=tdepth*(chi2[i]-chi2med)/(chi2inlow-chi2med);  // Converted depth value of the convolved data
    ntmp++;
    }
rmsmin = RMS(tmpdob1,ntmp);  // RMS assuming gaussian noise and excluding pri and sec events
fred = sqrt(nintrans)*STDEV(tmpdob2,ntmp)/rmsmin;  // Have to account for the number of points in-transit we're averaging over.


// Calculate Sigma_FA - the false alarm rate. I've set it so we shoudl only see one false alarm in 10,000 KOIs, assuming gaussian noise.
sigfa = sqrt(2)*INVERFC((tend-tstart)/(period*20000));  // Assume 20,000 KOIs


// If a sec, ter, or pos event was not found, set it to zero effectively
if(chi2outlow == 9E9)
  chi2outlow = chi2med;
if(chi2outhigh == -9E9)
  chi2outhigh = chi2med;
if(chi2terlow == 9E9)
  chi2terlow = chi2med;

// Calcualte the depths of the other events
depfacsec = (chi2outlow-chi2med)/(chi2inlow-chi2med);
depfacter = (chi2terlow-chi2med)/(chi2inlow-chi2med);
depfacpos = (chi2outhigh-chi2med)/(chi2inlow-chi2med);

// Significance of primary, secondary, tertiary, and positive events
depsig = rmsmin/sqrt(nintrans);  // Noise level within transit duration assuming gaussian noise
sigpri = (-tdepth*(chi2inlow-chi2med)/(chi2inlow-chi2med))/depsig;
sigsec = (-tdepth*(chi2outlow-chi2med)/(chi2inlow-chi2med))/depsig;
sigter = (-tdepth*(chi2terlow-chi2med)/(chi2inlow-chi2med))/depsig;
sigpos = (-tdepth*((2*chi2med-chi2outhigh)-chi2med)/(chi2inlow-chi2med))/depsig;

// If significance is really low just call it zero.
if(fabs(sigpri)<1E-2) sigpri=0;
if(fabs(sigsec)<1E-2) sigsec=0;
if(fabs(sigter)<1E-2) sigter=0;
if(fabs(sigpos)<1E-2) sigpos=0;


// Make all phases fit within -0.25 to 0.75
if(prilowtime<-0.25*period)
  prilowtime+=period;

if(seclowtime<-0.25*period)
  seclowtime+=period;

if(terlowtime<-0.25*period)
  terlowtime+=period;

if(sechightime<-0.25*period)
  sechightime+=period;


// Terminal output
cout << fixed << setprecision(10) << basename << " " << sigpri << " " << sigsec << " " << sigter << " " << sigpos << " " << sigfa << " " << fred << " " << prilowtime/period << " " << seclowtime/period << " " << terlowtime/period << " " << sechightime/period << " " << -depfacsec*tdepth << " " << depsig << endl;
//cout << fixed << setprecision(6) << basename <<  " " << depfacsec*tdepth << " " << depsig << endl;

// Uncoment for Plotting
// Output file for plotting
tmpstr1 = "outfile2-" + basename + ".dat";
outfile.open(tmpstr1.c_str());
for(i=0;i<ndat;i++)
  outfile << fixed << setprecision(10) << data[0][midti+i]/period << " " << rms[i] << " " << chi2[i] << " " << tdepth*(chi2[i]-chi2med)/(chi2inlow-chi2med) << endl;
outfile.close();

// Sort for plotting
tmpstr1 = "sort -k 1n outfile2-" + basename + ".dat > outfile3-" + basename + ".dat";
system(tmpstr1.c_str());

// Make a file for binning
syscmd = "awk '{print($1,$2)}' outfile1-" + basename + ".dat > " + basename +"-bininput.dat";
system(syscmd.c_str());

// Bin data for the binned data
convert << "./bin-noerr " << basename << "-bininput.dat " << setprecision(10) << ((tend-tstart)/period)/8 << " " << basename << "-binned2.dat"; syscmd = convert.str(); convert.str(""); convert.clear();
system(syscmd.c_str());

// Make plot file
tmpstr1 = basename + "-ModShift.gnu";
outfile.open(tmpstr1.c_str());
outfile << "set term pdfcairo enhanced dashed size 6in,8in" << endl; //pdfcairo enhanced size 5.5in,4.25in" << endl;
outfile << "set output '" << basename << "-modshift.pdf'" << endl;
outfile << "set multiplot layout 4,1 title \"TCE " << basename << ", P = " << setprecision(6) << period << " Days\" font ',20'" << endl;

// Make the table containing the results
outfile << "set arrow from screen 0.75,0.025 to screen 0.75,0.475 nohead lt 1 lw 5 lc 7" << endl;
outfile << "set arrow from screen 0.965,0.025 to screen 0.965,0.475 nohead lt 1 lw 5 lc 7" << endl;
outfile << "set arrow from screen 0.75,0.025 to screen 0.965,0.025 nohead lt 1 lw 5 lc 7" << endl;
outfile << "set arrow from screen 0.75,0.475 to screen 0.965,0.475 nohead lt 1 lw 5 lc 7" << endl;
outfile << "set arrow from screen 0.85,0.025 to screen 0.85,0.475 nohead lt 1 lw 5 lc 7" << endl;

outfile << "set arrow from screen 0.75,0.425 to screen 0.965,0.425 nohead lt 1 lw 5 lc 7" << endl;
outfile << "set arrow from screen 0.75,0.375 to screen 0.965,0.375 nohead lt 1 lw 5 lc 7" << endl;
outfile << "set arrow from screen 0.75,0.325 to screen 0.965,0.325 nohead lt 1 lw 5 lc 7" << endl;
outfile << "set arrow from screen 0.75,0.275 to screen 0.965,0.275 nohead lt 1 lw 5 lc 7" << endl;
outfile << "set arrow from screen 0.75,0.225 to screen 0.965,0.225 nohead lt 1 lw 5 lc 7" << endl;
outfile << "set arrow from screen 0.75,0.175 to screen 0.965,0.175 nohead lt 1 lw 5 lc 7" << endl;
outfile << "set arrow from screen 0.75,0.125 to screen 0.965,0.125 nohead lt 1 lw 5 lc 7" << endl;
outfile << "set arrow from screen 0.75,0.075 to screen 0.965,0.075 nohead lt 1 lw 5 lc 7" << endl;


outfile << "set label \"{/Symbol s}_{FA}\"          at screen 0.76,0.45 left font ', 14'" << endl;
outfile << "set label \"{/Symbol s}_{Pri}\"         at screen 0.76,0.40 left font ', 14'" << endl;
outfile << "set label \"{/Symbol s}_{Sec}\"         at screen 0.76,0.35 left font ', 14'" << endl;
outfile << "set label \"{/Symbol s}_{Ter}\"         at screen 0.76,0.30 left font ', 14'" << endl;
outfile << "set label \"{/Symbol s}_{Pos}\"         at screen 0.76,0.25 left font ', 14'" << endl;
outfile << "set label \"{/Symbol s}_{Pri-Sec}\"     at screen 0.76,0.20 left font ', 14'" << endl;
outfile << "set label \"{/Symbol s}_{Pri-Pos}\"     at screen 0.76,0.15 left font ', 14'" << endl;
outfile << "set label \"F_{Red}\"                   at screen 0.76,0.10 left font ', 14'" << endl;
outfile << "set label \"{/Symbol s}_{Pri}/F_{red}\" at screen 0.76,0.05 left font ', 14'" << endl;

outfile << "set label \"" << setprecision(1) << sigfa         << "\" at screen 0.9075,0.45 center font ',14' textcolor lt 7" << endl;
outfile << "set label \"" << setprecision(1) << sigpri        << "\" at screen 0.9075,0.40 center font ',14'";
if(sigpri<sigfa)        outfile << " textcolor lt 1" << endl; else outfile << " textcolor lt 7" << endl;
outfile << "set label \"" << setprecision(1) << sigsec        << "\" at screen 0.9075,0.35 center font ',14'";
if(sigsec>sigfa)        outfile << " textcolor lt 1" << endl; else outfile << " textcolor lt 7" << endl;
outfile << "set label \"" << setprecision(1) << sigter        << "\" at screen 0.9075,0.30 center font ',14'";
if(sigter>sigfa)        outfile << " textcolor lt 1" << endl; else outfile << " textcolor lt 7" << endl;
outfile << "set label \"" << setprecision(1) << sigpos        << "\" at screen 0.9075,0.25 center font ',14'";
if(sigpos>sigfa)        outfile << " textcolor lt 1" << endl; else outfile << " textcolor lt 7" << endl;
outfile << "set label \"" << setprecision(1) << sigpri-sigsec << "\" at screen 0.9075,0.20 center font ',14'";
if(sigpri-sigsec<3.0) outfile << " textcolor lt 1" << endl; else outfile << " textcolor lt 7" << endl;
outfile << "set label \"" << setprecision(1) << sigpri-sigpos << "\" at screen 0.9075,0.15 center font ',14'";
if(sigpri-sigpos<3.0) outfile << " textcolor lt 1" << endl; else outfile << " textcolor lt 7" << endl;
outfile << "set label \"" << setprecision(1) << fred          << "\" at screen 0.9075,0.10 center font ',14' textcolor lt 7" << endl;
outfile << "set label \"" << setprecision(1) << sigpri/fred   << "\" at screen 0.9075,0.05 center font ',14'";
if(sigpri/fred<sigfa)   outfile << " textcolor lt 1" << endl; else outfile << " textcolor lt 7" << endl;


// First plot
outfile << "set xlabel 'Phase'" << endl;
outfile << "set xrange [-0.25 to 1.25]" << endl;
outfile << "set x2range [-0.25 to 1.25]" << endl;
outfile << "set xtics 0.25" << endl;
outfile << "set format x '%4.2f'" << endl;
outfile << "set format y '%6.0f'" << endl;
outfile << "set ylabel 'Flux (ppm)'" << endl;
outfile << "set object rect from 0.75,-1E7 to 1.25,1E7 fc rgb '#D3D3D3' lw 0" << endl;
outfile << "set label '' at " << setprecision(10) << prilowtime/period  << ", graph 0.015 front point pt 9  ps 0.8" << endl;
outfile << "set label '' at " << setprecision(10) << seclowtime/period  << ", graph 0.015 front point pt 9  ps 0.8" << endl;
outfile << "set label '' at " << setprecision(10) << sechightime/period << ", graph 0.985 front point pt 11 ps 0.8" << endl;
outfile << "set label '' at " << setprecision(10) << terlowtime/period  << ", graph 0.015 front point pt 8  ps 0.8" << endl;
outfile << "stats '" << basename << "-binned2.dat' u 2 nooutput" << endl;
// outfile << "set yrange []" << endl;
outfile << "set yrange [1.0E6*STATS_min-0.5*(STATS_max-STATS_min) to 1.0E6*STATS_max+0.5*(STATS_max-STATS_min)]" << endl;
outfile << "set autoscale y" << endl;
outfile << "plot 'outfile1-" << basename << ".dat' u 1:($2*1.0E6) pt 7 ps 0.1 lc 1 notitle, '' u ($1+1.0):($2*1.0E6) pt 7 ps 0.1 lc 1 notitle, '' u ($1+2.0):($2*1.0E6) pt 7 ps 0.1 lc 1 notitle, '" << basename << "-binned2.dat' u 1:($2*1.0E6) pt 7 ps 0.1 lc 3 notitle, '' u ($1+1.0):($2*1.0E6) pt 7 ps 0.1 lc 3 notitle, 'outfile1-" << basename << ".dat' u 1:($3*1.0E6) with lines lt 1 lc 7 lw 5 notitle, '' u ($1+1.0):($3*1.0E6) with lines lt 1 lc 7 lw 5 notitle" << endl;

// Second Plot
outfile << "set autoscale y" << endl;
outfile << "set xlabel 'Phase'" << endl;
outfile << "set xtics format '%3.1f'" << endl;
outfile << "unset arrow" << endl;
outfile << "unset xlabel" << endl;
outfile << "unset label" << endl;
outfile << "set xtics ''" << endl;
outfile << "set x2tics 0.25 mirror" << endl;
outfile << "set format x2 '%4.2f'" << endl;
outfile << "set ylabel 'Flux (ppm)'" << endl;
outfile << "set label '' at " << setprecision(10) << prilowtime/period  << ", graph 0.015 front point pt 9  ps 0.8" << endl;
outfile << "set label '' at " << setprecision(10) << seclowtime/period  << ", graph 0.015 front point pt 9  ps 0.8" << endl;
outfile << "set label '' at " << setprecision(10) << sechightime/period << ", graph 0.985 front point pt 11 ps 0.8" << endl;
outfile << "set label '' at " << setprecision(10) << terlowtime/period  << ", graph 0.015 front point pt 8  ps 0.8" << endl;
outfile << "set object rect from 0.75,-1E7 to 1.25,1E7 fc rgb '#D3D3D3' lw 0" << endl;
outfile << "plot 'outfile3-" << basename << ".dat' u 1:(1E6*$4) with lines lt 1 lc 7 notitle, '' u ($1+1.0):(1E6*$4) with lines lt 1 lc 7 notitle, 0 with lines lt 2 lc 1 lw 5 notitle, " << setprecision(10) << sigfa*1E6*depsig << " with lines lt 2 lc 3 lw 5 notitle, " << -1.0*sigfa*1.0E6*depsig << " with lines lt 2 lc 3 lw 5 notitle" << endl;

outfile << "unset arrow" << endl;
outfile << "unset object" << endl;
outfile << "unset label" << endl;

// Pri Zoom
outfile << "set size square 0.41,0.41" << endl;
outfile << "set origin 0.0,0.17" << endl;
outfile << "set label 'Primary' at graph 0.5,0.925 center front" << endl;
outfile << "set xrange [" << setprecision(10) << 3*tstart/period << " to " << 3*tend/period << "]" << endl;
outfile << "set xtics " << setprecision(10) << (3*tend/period-3*tstart/period)/3.0 << endl;
outfile << "set xtics format ''" << endl;
outfile << "unset xlabel" << endl;
outfile << "set x2range [" << setprecision(10) << 3*tstart/period << " to " << 3*tend/period << "]" << endl;
outfile << "set x2tics " << setprecision(10) << (3*tend/period-3*tstart/period)/3.0 << endl;
outfile << "set x2tics format '%5.3f'" << endl;
outfile << "set x2label 'Phase'" << endl;
outfile << "unset y2tics" << endl;
outfile << "unset y2label" << endl;
outfile << "set autoscale y" << endl;
outfile << "set ytics format '%6.0f' mirror" << endl;
outfile << "set ytics auto" << endl;
outfile << "set ylabel 'Flux (ppm)'" << endl;
outfile << "plot '" << basename << "-binned2.dat' u 1:($2*1E6):($3*1E6) with yerrorbars lt 1 pt 7 ps 0.5 lc 3 notitle, '' u ($1+1.0):($2*1E6):($3*1E6) with yerrorbars lt 1 pt 7 ps 0.5 lc 3 notitle, '' u ($1-1.0):($2*1E6):($3*1E6) with yerrorbars lt 1 pt 7 ps 0.5 lc 3 notitle, 'outfile1-" << basename << ".dat' u 1:($3*1E6) with lines lt 1 lc 7 notitle, '' u ($1+1.0):($3*1E6) with lines lt 1 lc 7 notitle, '' u ($1-1.0):($3*1E6) with lines lt 1 lc 7 notitle" << endl;

outfile << "unset label" << endl;

// Sec Zoom
outfile << "set size square 0.41,0.41" << endl;
outfile << "set origin 0.375,0.17" << endl;
outfile << "set label 'Secondary' at graph 0.5,0.925 center front" << endl;
outfile << "set xrange [" << setprecision(10) << seclowtime/period+3*tstart/period << " to " << seclowtime/period+3*tend/period << "]" << endl;
outfile << "set xtics " << setprecision(10) << (3*tend/period-3*tstart/period)/3.0 << endl;
outfile << "set xtics format ''" << endl;
outfile << "unset xlabel" << endl;
outfile << "set x2range [" << setprecision(10) << seclowtime/period+3*tstart/period << " to " << seclowtime/period+3*tend/period << "]" << endl;
outfile << "set x2tics " << setprecision(10) << (3*tend/period-3*tstart/period)/3.0 << endl;
outfile << "set x2tics format '%5.3f' mirror" << endl;
outfile << "set x2label 'Phase'" << endl;
outfile << "unset ytics" << endl;
outfile << "unset ylabel" << endl;
outfile << "set autoscale y" << endl;
outfile << "set y2tics format '% -6.0f' mirror" << endl;
outfile << "set y2tics auto" << endl;
outfile << "unset ylabel" << endl;
// outfile << "set y2label 'Flux (ppm)'" << endl;
outfile << "plot '" << basename << "-binned2.dat' u 1:($2*1E6):($3*1E6) with yerrorbars lt 1 pt 7 ps 0.5 lc 3 notitle, '' u ($1+1.0):($2*1E6):($3*1E6) with yerrorbars lt 1 pt 7 ps 0.5 lc 3 notitle, '' u ($1-1.0):($2*1E6):($3*1E6) with yerrorbars lt 1 pt 7 ps 0.5 lc 3 notitle, 'outfile1-" << basename << ".dat' u ($1+" << setprecision(10) << seclowtime/period << "):(" << depfacsec << "*$3*1E6) with lines lt 1 lc 7 notitle, '' u ($1+1.0+" << seclowtime/period << "):(" << depfacsec << "*$3*1E6) with lines lt 1 lc 7 notitle, '' u ($1-1.0+" << seclowtime/period << "):(" << depfacsec << "*$3*1E6) with lines lt 1 lc 7 notitle" << endl;

outfile << "unset label" << endl;

// Ter Zoom
outfile << "set size square 0.41,0.41" << endl;
outfile << "set origin 0.0,-0.075" << endl;
outfile << "set label 'Tertiary' at graph 0.5,0.925 center front" << endl;
outfile << "set xtics format '%5.3f' mirror" << endl;
outfile << "set xtics " << setprecision(10) << (3*tend/period-3*tstart/period)/3.0 << endl;
outfile << "set xrange [" << setprecision(10) << terlowtime/period+3*tstart/period << " to " << terlowtime/period+3*tend/period << "]" << endl;
outfile << "set xlabel 'Phase'" << endl;
outfile << "unset x2tics" << endl;
outfile << "unset x2label" << endl;
outfile << "unset y2tics" << endl;
outfile << "unset y2label" << endl;
outfile << "set autoscale y" << endl;
outfile << "set ytics format '%6.0f' mirror" << endl;
outfile << "set ytics auto" << endl;
outfile << "set ylabel 'Flux (ppm)'" << endl;
outfile << "plot '" << basename << "-binned2.dat' u 1:($2*1E6):($3*1E6) with yerrorbars lt 1 pt 7 ps 0.5 lc 3 notitle, '' u ($1+1.0):($2*1E6):($3*1E6) with yerrorbars lt 1 pt 7 ps 0.5 lc 3 notitle, '' u ($1-1.0):($2*1E6):($3*1E6) with yerrorbars lt 1 pt 7 ps 0.5 lc 3 notitle, 'outfile1-" << basename << ".dat' u ($1+" << setprecision(10) << terlowtime/period << "):(" << depfacter << "*$3*1E6) with lines lt 1 lc 7 notitle, '' u ($1+1.0+" << terlowtime/period << "):(" << depfacter << "*$3*1E6) with lines lt 1 lc 7 notitle, '' u ($1-1.0+" << terlowtime/period << "):(" << depfacter << "*$3*1E6) with lines lt 1 lc 7 notitle" << endl;

outfile << "unset label" << endl;

// Pos Zoom
outfile << "set size square 0.41,0.41" << endl;
outfile << "set origin 0.375,-0.075" << endl;
outfile << "set label 'Positive' at graph 0.5,0.075 center front" << endl;
outfile << "set xtics format '%5.3f' mirror" << endl;
outfile << "set xtics " << setprecision(10) << (3*tend/period-3*tstart/period)/3.0 << endl;
outfile << "set xrange [" << setprecision(10) << sechightime/period+3*tstart/period << " to " << sechightime/period+3*tend/period << "]" << endl;
outfile << "set xlabel 'Phase'" << endl;
outfile << "unset x2tics" << endl;
outfile << "unset x2label" << endl;
outfile << "unset ytics" << endl;
outfile << "unset ylabel" << endl;
outfile << "set autoscale y" << endl;
outfile << "set y2tics format '% -6.0f' mirror" << endl;
outfile << "set y2tics auto" << endl;
outfile << "unset ylabel" << endl;
// outfile << "set y2label 'Flux (ppm)'" << endl;
outfile << "plot '" << basename << "-binned2.dat' u 1:($2*1E6):($3*1E6) with yerrorbars lt 1 pt 7 ps 0.5 lc 3 notitle, '' u ($1+1.0):($2*1E6):($3*1E6) with yerrorbars lt 1 pt 7 ps 0.5 lc 3 notitle, '' u ($1-1.0):($2*1E6):($3*1E6) with yerrorbars lt 1 pt 7 ps 0.5 lc 3 notitle, 'outfile1-" << basename << ".dat' u ($1+" << setprecision(10) << sechightime/period << "):(" << depfacpos << "*$3*1E6) with lines lt 1 lc 7 notitle, '' u ($1+1.0+" << sechightime/period << "):(" << depfacpos << "*$3*1E6) with lines lt 1 lc 7 notitle, '' u ($1-1.0+" << sechightime/period << "):(" << depfacpos << "*$3*1E6) with lines lt 1 lc 7 notitle" << endl;

// And make the whole thing
outfile << "unset multiplot" << endl;
tmpstr1 = "gnuplot " + basename + "-ModShift.gnu";
system(tmpstr1.c_str());


// Clean up files
syscmd = "cp outfile1-" + basename + ".dat " + basename + ".cln";  // Make clean file for use later
system(syscmd.c_str());
syscmd = "rm outfile?-" + basename + ".dat";
system(syscmd.c_str());
syscmd = "rm " + basename + "-bin*dat";
system(syscmd.c_str());
syscmd = "rm " + basename + "-ModShift.gnu";
system(syscmd.c_str());

// cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC << " seconds." << endl;
// */

return 0;
}

///////////////////////////////////////////////////////////////////////////////
double MAD(double x[], int n)  // Median Absolute Deviation
  {
  double xmed,w[N];
  int z;
  
  xmed = MEDIAN(x,n);
  for(z=0;z<n;z++)
    w[z] = fabs(x[z] - xmed);
  
  return MEDIAN(w,n);
  }

///////////////////////////////////////////////////////////////////////////////
double INVERFC(double p)
  {
  double x,err,t,pp;
  if(p >= 2.0) return -100.0;
  if(p <= 0.0) return 100.0;
  pp=(p<1.0)? p : 2.0-p;
  t = sqrt(-2.*log(pp/2.0));
  x = -0.70711*((2.30753+t*0.27061)/(1.0+t*(0.99229+t*0.04481)) - t);
  for(int j=0;j<2;j++)
    {
      err=erfc(x)-pp;
      x+=err/(M_2_SQRTPI*exp(-pow(x,2))-x*err);
    }
  return(p<1.0? x : -x);
  }

///////////////////////////////////////////////////////////////////////////////

double MEDIAN(double x[], int n)
  {
  double w[N];
  int z;
  
  for(z=0;z<n;z++)
    w[z]=x[z];
  
  qsort(w,n,sizeof(double),cmp);

  if(n%2==0)
    return (w[n/2]+w[n/2-1])/2;  // Find median
  else
    return w[n/2];
  }

///////////////////////////////////////////////////////////////////////////////

double RMS(double y[],int a) {
/*-------------------------------------------------------------
  Calculate the root mean squared of an array y, given a terms
  -------------------------------------------------------------*/

double sum=0;
int z;

for(z=0;z<a;z++)
  sum+=pow(y[z],2);

return sqrt(sum/a);
}


//////////////////////////////////

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

////////////////////////////////////////////////////////////////////////////////////////////////////

double SUM(double x[],int a) {  // Return sum of first a terms

int z;
double xsum;

xsum=0;
for(z=0;z<a;z++)
  xsum += x[z];

return(xsum);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

double SUMSQ(double x[],int a) {  // Return sum of first a terms squared

int z;
double xsum;

xsum=0;
for(z=0;z<a;z++)
  xsum += x[z]*x[z];

return(xsum);
}

///////////////////////////////////////////////////////////////////////////////

double SELECT(double x[],int n, int m)    // Returns the mth smallest value in array x, given n terms in array x
  {
  int z;
  double w[N];
  
  for(z=0;z<n;z++)
    w[z]=x[z];
  
  qsort(w,n,sizeof(double),cmp);
    
  return w[m];
  }



