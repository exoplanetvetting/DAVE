/* Jeff Coughlin
   modshift-v4.cpp

Last modified: 2016-05-12
   
*/


// Include needed libraries
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include <algorithm>
#include <sstream>

using namespace std;


// Declare constants
const int N=100000;  // Max number of data points. Increase if needed.


// Declare Functions
void DO_SHIFT(),PLOT(),BIN_NOERR(string, double, string), DMM(), OE();
double STDEV(double[],int),RMS(double[],int),MEDIAN(double[],int),MAX(double[],int),MIN(double[],int),INVERFC(double),MAD(double[],int),MEAN(double[],int),MODSCALE(double[],double[],int);
string FORMAT(double);

struct datastruct {double time; double phase; double flux; double model; double resid;};  // Time, Phase, Flux, Model, Residual

// Declare global variables (Using global for large array so memory is statically allocated for these large arrays.)

datastruct data[2*N],dataorig[2*N];  // The main data array
double rms[N],chi2[N];
double tmpdobarry1[N],tmpdobarry2[N],tmpdob3[N];  // Generic temporary double arrays
int tmpintarry[N];  // Array of ints for temporary use
double flat[2*N]; // residual assuming flat baseline
datastruct inpdata[N];  // Input data for binning
datastruct bindat[N];   // Output binned data
double convolveddepth[N]={0};  // Array of the convovled depth values

struct resultstruct {double sigpri=0; double sigsec=0; double sigter=0; double sigpos=0; double sigfa1=1; double sigfa2=1; double fred=1; double prilowtime=0; double seclowtime=0; double terlowtime=0; double sechightime=0; double depfacsec=0; double depfacter=0; double depfacpos=0; double depsig=1; double tdepth=0; double odddepth=0; double evndepth=0; double odddeptherr=1; double evndeptherr=1; double sigoe=0; double depthmedmeanrat=0; double shape=0; double duration=0; double baseflux=0; double fullmoddepth=0;};
resultstruct results;

int ndat,ndatorig;
double period,epoch,periodorig,epochorig;
string basename,infilename,objectname;

int i,j,k,l,m,n; // General purpose counting integers
int sw1,sw2;  // Switches
ifstream infile;  // Generic ifstream
ofstream outfile;  // Generic ofstream
string syscmd;  // Generic system command string
int modplotint; // 1 = do plot, else no plot

double tstart,tend;  // Start and end times of model fit
string tmpstr1;  // Generic string for names
double tmpdob1;  // Generic temporary double
int tmpint1;

// Sorting functions
bool double_sort (double lhs, double rhs) {return lhs < rhs;}
bool phase_sort (datastruct lhs, datastruct rhs) {return lhs.phase < rhs.phase;}
bool time_sort (datastruct lhs, datastruct rhs) {return lhs.time < rhs.time;}



// And main prog here
int main (int argc, char* argv[]) {


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
  objectname = argv[3];
else
  {
  cout << "Name of object? ";
  cin >> objectname;
  }
  

if(argc>4)
  period = atof(argv[4]);
else
  {
  cout << "Period? ";
  cin >> period;
  }
  
if(argc>5)
  epoch = atof(argv[5]);
else
  {
  cout << "Epoch? ";
  cin >> epoch;
  }

if(argc>6)
  modplotint = atoi(argv[6]);
else
  {
  cout << "Plot? (1==yes, else no) ";
  cin >> modplotint;
  }



// Open input file and check to make sure file exists
infile.open(infilename.c_str());
if(!infile.good())
  {
  cout << "Input file does not exist or cannot read! Exiting..." << endl;
  exit(1);
  }
  
// Read in file
i=0;
infile >> data[i].time;
while(!infile.eof())
  {
  infile >> data[i].flux;
  infile >> data[i].model;
  i++;
  infile >> data[i].time;
  }
infile.close();
ndat = i;  // Record original number of datapoints

// Sort input data on time
sort(data,data+ndat,time_sort);

// Make sure epoch occurs before earliest data point
while(epoch>data[0].time) 
  epoch-=period;

// Back up original data
ndatorig = ndat;
periodorig = period;
epochorig = epoch;
for(i=0;i<ndat;i++)
  dataorig[i] = data[i];

// Phase Data given period and epoch
for(i=0;i<ndat;i++)
  {
  data[i].phase = period*((data[i].time-epoch)/period - int((data[i].time-epoch)/period));
  if(data[i].phase > 0.75*period)  // Make so ranges from -0.25*period to 0.75*period
    data[i].phase-=period;
  }

  
// Sort data on phase
sort(data,data+ndat,phase_sort);

// Record width in phase space
sw1=0;
for(i=0;i<ndat;i++)
  {
  if(i==0)
    results.baseflux = data[i].model; // Out of transit baseline flux of the model
  
  if(sw1==0 && data[i].model!=results.baseflux)
    {
    tstart=data[i].phase;
    sw1=1;
    }
  if(sw1==1 && data[i].model==results.baseflux)
    {
    tend=data[i-1].phase;
    sw1=2;
    }
  if(sw2==0 && data[i].phase > 0)
    {
    tmpint1 = i;
    sw2=1;
    }
  i++;
  }

// Check for odd cases of no model data before or after midpoint
if(tstart>0) // If tstart is positive then there is no model before center, so go off the end time
  tstart = -1.0*tend;
if(tend<0)  // If tend is negative, then there is no model after center of transit, so go off end time.
  tend = -1.0*tstart;

// Record duration
results.duration = tend-tstart;

// Record nominal model depth
results.fullmoddepth = fabs(data[tmpint1].model - results.baseflux);


// Run Odd-Even Test
OE();

// Run DMM test
DMM();


// Okay now do it normally
ndat = ndatorig;  // Restore original ndat
for(i=0;i<ndat;i++)  // Load original data
  data[i] = dataorig[i];
period = periodorig;  // Load original period
epoch = epochorig;  // Load original epoch

DO_SHIFT();  // Run shift


if(MAX(convolveddepth,ndat) - MIN(convolveddepth,ndat) != 0)  // Add safety check to avoid divide by zero
  results.shape = MAX(convolveddepth,ndat) / (MAX(convolveddepth,ndat) - MIN(convolveddepth,ndat));   // Same as above assuming a median of 0.0


// Terminal output
cout << basename << " " << fixed << setprecision(10) << results.sigpri << " " << results.sigsec << " " << results.sigter << " " << results.sigpos << " " << results.sigoe << " " << results.depthmedmeanrat << " " << results.shape << " " << results.sigfa1 << " " << results.sigfa2 << " " << results.fred << " " << results.prilowtime/period << " " << results.seclowtime/period << " " << results.terlowtime/period << " " << results.sechightime/period << " " << -results.depfacsec*results.tdepth << " " << results.depsig << endl;


PLOT();

// cout << double( clock() - startTime ) / (double)CLOCKS_PER_SEC << " seconds." << endl;

return 0;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////

void DMM()  // Compute individual depths for depth mean to median ratio test
  {
  // Depth mean to median params here
  int ninda,nindb;
  double indivdepths[5000]; // Big just in case - max possible should be ~2800 given 0.5 day TCE and 1400 day baseline
  double depthsmed,depthsmean;
  double dmmwidth;
  



  period = periodorig;  // Load original period
  epoch = epochorig;  // Load original epoch

  ninda=0; // Index to loop through all possible individual events
  nindb=0; // index to keep track of indiv measurements, for transits which actually exist in the data

  dmmwidth = 1.5*results.duration;  // Use data that is within 1.5 transit duration of transit, i.e., transit width either side for baseline. Units are in days.

  while(epoch + ninda*period <= dataorig[ndatorig-1].time)  // Go through all possible individual transit events
    {
    ndat=0;
    for(i=0;i<ndatorig;i++)  // Load from original data...
      if(dataorig[i].time > epoch + ninda*period - dmmwidth && dataorig[i].time < epoch + ninda*period + dmmwidth)  // ...if it falls within dmmwidth of an indidivual transit
        {
        data[ndat] = dataorig[i];
        ndat++;
        }
      
    // Check to see if model is all flat or not, i.e., if there is data in-transit
    j=0;
    for(i=0;i<ndat;i++)
      if(data[i].model!=data[0].model)
        j=1;  // a j value of 1 means the model is not all flat, and we can proceed
    
    // If there is data, and data in-transit, compute the depth of the individual transit
    if(ndat>0 && j==1)
      {
      // Scale model to fit the data
      for(i=0;i<ndat;i++)
        {
        tmpdobarry1[i] = data[i].flux;   // This is just so the data can feed into MODSCALE
        tmpdobarry2[i] = data[i].model;  // This is just so the data can feed into MODSCALE
        }
      tmpdob1 = MODSCALE(tmpdobarry1,tmpdobarry2,ndat);
    
      // Record model depth of primary transit
      indivdepths[nindb] = tmpdob1*results.fullmoddepth;  // -fabs(data[midti].model - results.baseflux);
      
      // Advance counter for number of successsful individual event measurements
      nindb++;
      }
    
    ninda++;  // Try looking at next indidivual event
    }

  depthsmed = MEDIAN(indivdepths,nindb);
  depthsmean = MEAN(indivdepths,nindb);
  results.depthmedmeanrat = depthsmean/depthsmed;
  }
  
  
//////////////////////////////////////////////////////////////////////////////////////////////////////

void OE()  // Compute odd/even statistic
  {
  // Compute for Odd transits only first
  ndat = ndatorig;  // Restore original ndat
  for(i=0;i<ndat;i++)  // Load original data
    data[i] = dataorig[i];
  period = periodorig;  // Load original period
  epoch = epochorig;  // Load original epoch

  // epoch-=0.5*period
  period*=2;   // Compute phase assuming twice period
  for(i=0;i<ndat;i++)  // Phase data
    {
    data[i].phase = ((data[i].time-epoch)/period - int((data[i].time-epoch)/period));
    if(data[i].phase > 0.75)
      data[i].phase-=1.0;  // Make pahse go from -0.25 to 0.75
    }

  sort(data,data+ndat,phase_sort);  // Sort data on phase

  i=0;
  while(data[i].phase < 0.25)  // Find data that now have phase less than 0.5 - these are the odd transits
    i++;
  ndat = i;

  period = periodorig;  // Set period back to normal
  DO_SHIFT();  // Run shift
  results.odddepth = results.tdepth;
  results.odddeptherr = results.depsig;

  syscmd = "mv " + basename + "-outfile1.dat " + basename + "-outfile1-odd.dat";
  system(syscmd.c_str());

  

  // Now compute for even only
  ndat = ndatorig;  // Restore original ndat
  for(i=0;i<ndat;i++)  // Load original data
    data[i] = dataorig[i];
  period = periodorig;  // Load original period
  epoch = epochorig;  // Load original epoch

  epoch-=period;   // Select only even transits. Doing minus so epoch stays earlier than first data time
  period*=2;   // Compute phase assuming twice period
  for(i=0;i<ndat;i++)  // Phase data
    {
    data[i].phase = ((data[i].time-epoch)/period - int((data[i].time-epoch)/period));
    if(data[i].phase > 0.75)
      data[i].phase-=1.0;  // Make phase go from -0.25 to 0.75
    }

  sort(data,data+ndat,phase_sort);  // Sort data on phase

  i=0;
  while(data[i].phase < 0.25)  // Find data that now have phase less than 0.5 - these are the even transits, since we've shifted epoch above
    i++;
  ndat = i;

  period = periodorig;  // Set period back to normal
  DO_SHIFT();  // Run shift
  results.evndepth = results.tdepth;
  results.evndeptherr = results.depsig;
  syscmd = "mv " + basename + "-outfile1.dat " + basename + "-outfile1-evn.dat";
  system(syscmd.c_str());


  
  // Compute odd/even sigma
  if(results.odddepth!=0 && results.evndepth!=0)
    results.sigoe = fabs(results.odddepth - results.evndepth)/sqrt(pow(results.odddeptherr,2) + pow(results.evndeptherr,2));
  }
  
  

/////////////////////////////////////////////////////////////////////////////////////

// Scale model so that it best-fits the data and return the scale factor
// Not the best written search algorithm, but works fast enough I think.

double MODSCALE(double flux[], double model[], int ndat)
  {
  double modscale,rmsscale,modscaletmp,rmsorig,rmsscaletmp,dinc;
//   int sw1; // 0 is try up, 1 is try down
  
  // Starting values
  modscale = modscaletmp = 1.0;
  dinc = 1.0;
  sw1=0;
  
  // Compute current rms
  for(i=0;i<ndat;i++)
    tmpdob3[i] = flux[i] - modscale*model[i];  // Calculate residuals
  rmsscale = RMS(tmpdob3,ndat);  // Compute rms of residuals

  while(dinc>1E-8)  // 1E-8 seems like plenty of good accuracy not that long to converge. Allows for ppm level precision.
    {
    if(sw1==0)
      modscaletmp += dinc;  // Try going positive
    if(sw1==1)
      modscaletmp -= dinc;  // Try going negative
      
    for(i=0;i<ndat;i++)
      tmpdob3[i] = flux[i] - modscaletmp*model[i];  // Calculate residuals
    rmsscaletmp = RMS(tmpdob3,ndat);  // Compute rms of residuals
        
    if(rmsscaletmp < rmsscale)  // If we found a better fit
      {
      modscale = modscaletmp;  // Update modscale with that of best-fit
      rmsscale = rmsscaletmp;  // Update rmsscale with that of best-fit
      if(sw1==0)
        sw1=-1;  // Set switch to -1. Will get incremented to 0 in next iteration. Keeps searching upward
      if(sw1==1)
        sw1=0;   // Set switch to 0. Will get incremented to 1 in next iteration. Keeps searching downward
      }
    else // If better fit was not found, set modescaletmp back to what it was before it was modified
      {
      if(sw1==0)
        modscaletmp -= dinc;
      if(sw1==1)
        modscaletmp += dinc;
      }
      
    sw1++;  // Increment switch
    
    if(sw1==2) // Neither up nor down produced  better fit, search at a smaller increment
      {
      dinc /= 10;  // Use smaller increment
      sw1=0;  // Will start next smaller search going positive
      }
    }
    
  return modscale;  // Return the thing we were fitting for
  }


////////////////////////////////////////////////////////////////////////////////////

void DO_SHIFT()
  {
  double rmsmin=9E9,chi2med,chi2outlow,chi2outhigh,chi2inlow,chi2low,chi2terlow;
  int ntmp;
  int midti,startti,endti;  // Index for middle of transit point, start of transit, and end of transit
  double widthfac;
  double width,halfwidth,detrenddur;
  double nintrans;
  double med,std,sigreject,mad;
  double tmpsum1,tmpsum2;
  
//   clock_t startTime = clock();  
  
  // Phase Data given period and epoch
  for(i=0;i<ndat;i++)
    {
    data[i].phase = period*((data[i].time-epoch)/period - int((data[i].time-epoch)/period));
    if(data[i].phase > 0.75*period)  // Make so ranges from -0.25*period to 0.75*period
      data[i].phase-=period;
    }

  // Sort data on phase
  sort(data,data+ndat,phase_sort);

  // Find middle index
  sw1=0;
  for(i=0;i<ndat;i++)
    {
    if(sw1==0 && data[i].phase > 0)
      {
      midti = i;
      sw1=1;
      }
    i++;
    }

  // Check to make sure model isn't all flat
  j=0;
  for(i=0;i<ndat;i++)
    if(data[i].model!=results.baseflux)
      j++;
    
  if(j==0)
    return;  // Exit function since there are no points in-transit
    
 
  // Scale model to fit data best (Should be the case, but may not due to odd/even stuff.)
  for(i=0;i<ndat;i++)
    {
    tmpdobarry1[i] = data[i].flux;;  // Calculate residuals
    tmpdobarry2[i] = data[i].model;
    }
  tmpdob1 = MODSCALE(tmpdobarry1,tmpdobarry2,ndat);
  
  for(i=0;i<ndat;i++)
    data[i].model *= tmpdob1;
    

  // Record model depth of primary transit
  results.tdepth = -fabs(data[midti].model - results.baseflux);


  // Remove very large outliers
  // detrenddur = 0.1*results.duration;  // I used to detrend at one tenth of a tranist duration
  // if(detrenddur < 0.02043981)
  detrenddur = 0.02043981;  // Set detrending window to one cadence duration. Try for now
   
  for(l=0;l<1;l++)  // Number of outlier passes - ONLY NEED ONE WHEN USING MAD - JUST LEAVING LOOP IN CASE I EVER WANT TO CHANGE MY MIND
    {
    for(i=0;i<ndat;i++)
      data[i].resid = data[i].flux - data[i].model;  // Calculate residuals

  // std = STDEV(data.resid,ndat);  // Saving this line in case I go back to std instead of mad
    for(i=0;i<ndat;i++)
      tmpdob3[i] = data[i].resid;
    mad = MAD(tmpdob3,ndat);  // MAD is Median Absolute Devitation - better than STD - only need one pass.
    std = 1.4826*mad;  // https://www.mathworks.com/help/stats/mad.html  AND https://en.wikipedia.org/wiki/Median_absolute_deviation
      
    sigreject = sqrt(2)*INVERFC(1.0/ndat);  // Calculate sigma threshold based on number of data points
      
    for(i=0;i<ndat;i++)  // Keep track of which data points to remove. 0 is good, if flagged as 1, remove later
      tmpintarry[i]=0;

    for(i=0;i<ndat;i++)
      {
      k=0;
      sw1=0;
      for(j=0;j<ndat;j++)
        if(fabs(data[i].phase-data[j].phase) < detrenddur)  // look at points within the detrend duration
          {
          tmpdobarry1[k]=data[j].resid;
          k++;
          }
      med = MEDIAN(tmpdobarry1,k);
      
      if(fabs(data[i].resid - med) > sigreject*std)  // Keep data points within sigreject standard deviations of median
        tmpintarry[i] = 1;
      }

    m=0;
    for(i=0;i<ndat;i++)
      if(tmpintarry[i]==0)  // Keep good data points
        {
        data[m].time = data[i].time;
        data[m].phase = data[i].phase;
        data[m].flux = data[i].flux;
        data[m].model = data[i].model;
        data[m].resid = data[i].resid;
        m++;
        }
      else
        if(m<midti)
          midti--;  // If I remove a data point, I have to update midti to keep track of it
    ndat=m;  // Update number of data points to those not thrown out
    }

    
  // Double up input data for shifting
  for(i=0;i<ndat;i++)
    {
    data[i+ndat].time = data[i].time;
    data[i+ndat].phase = data[i].phase;
    data[i+ndat].flux = data[i].flux;
    data[i+ndat].model = data[i].model;
    data[i+ndat].resid = data[i].resid;
    }  
  
  
  halfwidth = 0.5*results.duration;  // Pre-computing so don't have to do inside comp. intensive loop

  // Record the indexs of points that correspond to start and end of transit.
  sw1=0;
  for(j=0;j<ndat;j++)
    { 
    if(sw1==0 && data[j].phase > -halfwidth)
      {
      startti = j-1;  // Cadence of the last point before the transit occurs
      sw1=1;
      }
    if(sw1==1 && fabs(data[j].phase) > halfwidth)
      {
      endti = j;  // Cadence of the point just after the transit ends.
      sw1=2;
      }
    }
  
  
  // To speed things up, assume outside of transit is flat
  for(j=0;j<2*ndat;j++)
    flat[j] = pow(data[j].flux - results.baseflux,2);

  // Okay do the actual perumtation. 
  for(i=0;i<ndat;i++)  // Perform ndat pertubations
    {
    tmpsum2 = 0;
  
    // Before transit, can look up values for compuatation speed increase
    for(j=0;j<startti;j++)
      tmpsum2 += flat[j+i];
    
    // Compute new values inside transit
    for(j=startti;j<endti;j++)
      tmpsum2 += pow(data[j+i].flux - data[j].model,2);  // Shitfing data, holding model steady. Moving data points backwards, or model forwards, same thing
      
    // After transit, can look up values for computation speed increase
    for(j=endti;j<ndat;j++)
      tmpsum2 += flat[j+i];
    
    rms[i] = sqrt(tmpsum2/ndat);  // RMS of the new residuals
    if(rms[i] < rmsmin)
      rmsmin = rms[i];
    }
  
  
  // Now look at convolved data to find pri, sec, etc. and do other things
  nintrans=0; 
  chi2low = 9E9;
  tmpstr1 = basename + "-outfile1.dat";
  outfile.open(tmpstr1.c_str());
  for(i=0;i<ndat;i++)  // Calculate chi squared and find lowest chi squared value
    {
    chi2[i] = ndat*pow(rms[i],2)/pow(rmsmin,2);  // Using rms to calculate chi squared
    if(chi2[i] < chi2low)
      chi2low = chi2[i];
    outfile << fixed << setprecision(10) << data[i].phase/period << " " << setprecision(8) << setw(11) << data[i].flux << " " << data[i].model << endl; // Write out original data while we're at it
    if(data[i+midti].phase > tstart && data[i+midti].phase < tend)  // If within the actual transit window, compute number of in-transit data points while we're at it
      nintrans++;
    }
  outfile.close();

  // Search for sec. eclipse and compute in-primary params
  chi2outlow = chi2inlow = 9E9;
  results.seclowtime = -2*period;
  results.prilowtime = 0;
  widthfac=4;
  for(i=0;i<ndat;i++)
    {
    if((data[i+midti].phase < widthfac*tstart || data[i+midti].phase > widthfac*tend) && (data[i+midti].phase < widthfac*tstart+period || data[i+midti].phase > widthfac*tend+period) && (data[i+midti].phase < widthfac*tstart-period || data[i+midti].phase > widthfac*tend-period))  // If outside primary transit
      {
      if(chi2[i]<chi2outlow)
        {
        chi2outlow=chi2[i];
        results.seclowtime = data[i+midti].phase;
        }
      }
    else
      {
      if(chi2[i]<chi2inlow)
        {
        chi2inlow=chi2[i];
        results.prilowtime = data[i+midti].phase;
        }
      }
    }
  // Make sure we find a secondary no matter what
  if(results.seclowtime == -2*period)  // If no secondary was found
    for(i=0;i<ndat;i++)
      if(data[i+midti].phase > 0.4*period && data[i+midti].phase < 0.6*period)
        {
        if(chi2[i]<chi2outlow)
          {
          chi2outlow=chi2[i];
          results.seclowtime = data[i+midti].phase;
          }
        }
  
  // Find tertiary and biggest positive peak that's outside primary and secondary eclipse
  chi2terlow = 9E9;
  results.terlowtime = -2*period;
  widthfac=4;
  for(i=0;i<ndat;i++)
    {
    if((data[i+midti].phase < results.prilowtime+widthfac*tstart || data[i+midti].phase > results.prilowtime+widthfac*tend) && (data[i+midti].phase < results.prilowtime+widthfac*tstart+period || data[i+midti].phase > results.prilowtime+widthfac*tend+period) && (data[i+midti].phase < results.prilowtime+widthfac*tstart-period || data[i+midti].phase > results.prilowtime+widthfac*tend-period) && (data[i+midti].phase < results.seclowtime+widthfac*tstart || data[i+midti].phase > results.seclowtime+widthfac*tend) && (data[i+midti].phase < results.seclowtime+widthfac*tstart+period || data[i+midti].phase > results.seclowtime+widthfac*tend)+period && (data[i+midti].phase < results.seclowtime+widthfac*tstart-period || data[i+midti].phase > results.seclowtime+widthfac*tend-period)  )  // If outside primray and secondary     
      if(chi2[i]<chi2terlow)
        {
        chi2terlow=chi2[i];
        results.terlowtime=data[i+midti].phase;
        }
    }

  // Find biggest positive peak that's outside primary and secondary eclipse
  chi2outhigh = -9E9;
  results.sechightime = -2*period;
  widthfac=6;
  for(i=0;i<ndat;i++)
    {
    if((data[i+midti].phase < results.prilowtime+widthfac*tstart || data[i+midti].phase > results.prilowtime+widthfac*tend) && (data[i+midti].phase < results.prilowtime+widthfac*tstart+period || data[i+midti].phase > results.prilowtime+widthfac*tend+period) && (data[i+midti].phase < results.prilowtime+widthfac*tstart-period || data[i+midti].phase > results.prilowtime+widthfac*tend-period) && (data[i+midti].phase < results.seclowtime+widthfac*tstart || data[i+midti].phase > results.seclowtime+widthfac*tend) && (data[i+midti].phase < results.seclowtime+widthfac*tstart+period || data[i+midti].phase > results.seclowtime+widthfac*tend)+period && (data[i+midti].phase < results.seclowtime+widthfac*tstart-period || data[i+midti].phase > results.seclowtime+widthfac*tend-period)  )  // If outside primray and secondary     
      if(chi2[i]>chi2outhigh)
        {
        chi2outhigh=chi2[i];
        results.sechightime = data[i+midti].phase;
        }
    }

  // Compute median excluding primary and secondary eclipse. I have built in contigency in cases of sparse data so this should never fail.  
  widthfac=2;
  ntmp=0;
  while(ntmp<0.1*ndat)  // Require at least 3 points to exist out of transit & eclipse
    {
    ntmp=0;
    for(i=0;i<ndat;i++)
      if((data[i+midti].phase < results.prilowtime+widthfac*tstart || data[i+midti].phase > results.prilowtime+widthfac*tend) && (data[i+midti].phase < results.prilowtime+widthfac*tstart+period || data[i+midti].phase > results.prilowtime+widthfac*tend+period) && (data[i+midti].phase < results.prilowtime+widthfac*tstart-period || data[i+midti].phase > results.prilowtime+widthfac*tend-period) && (data[i+midti].phase < results.seclowtime+widthfac*tstart || data[i+midti].phase > results.seclowtime+widthfac*tend) && (data[i+midti].phase < results.seclowtime+widthfac*tstart+period || data[i+midti].phase > results.seclowtime+widthfac*tend)+period && (data[i+midti].phase < results.seclowtime+widthfac*tstart-period || data[i+midti].phase > results.seclowtime+widthfac*tend-period)  )  // If outside primray and secondary     
        {
        tmpdobarry1[ntmp]=chi2[i];
        ntmp++;
        }
      widthfac-=0.01;  // Reduce width factor by -0.1 and search again if no points were found out of transit & eclipse
    }
  chi2med = MEDIAN(tmpdobarry1,ntmp);
      
  // Calculate rmsmin and Fred, ratio of gaussian noise to systematic noise, excluding primary and secondary
  widthfac=2;
  ntmp=0;
  while(ntmp<0.1*ndat)  // Require at least 3 points to exist out of transit & eclipse
    {
    ntmp=0;
    for(i=0;i<ndat;i++)
      if((data[i+midti].phase < results.prilowtime+widthfac*tstart || data[i+midti].phase > results.prilowtime+widthfac*tend) && (data[i+midti].phase < results.prilowtime+widthfac*tstart+period || data[i+midti].phase > results.prilowtime+widthfac*tend+period) && (data[i+midti].phase < results.prilowtime+widthfac*tstart-period || data[i+midti].phase > results.prilowtime+widthfac*tend-period) && (data[i+midti].phase < results.seclowtime+widthfac*tstart || data[i+midti].phase > results.seclowtime+widthfac*tend) && (data[i+midti].phase < results.seclowtime+widthfac*tstart+period || data[i+midti].phase > results.seclowtime+widthfac*tend)+period && (data[i+midti].phase < results.seclowtime+widthfac*tstart-period || data[i+midti].phase > results.seclowtime+widthfac*tend-period)  )  // If outside primary and secondary     
        {
        tmpdobarry1[ntmp]=data[i+midti].flux-data[i+midti].model;  // Original data
        tmpdobarry2[ntmp]=results.tdepth*(chi2[i]-chi2med)/(chi2inlow-chi2med);  // Converted depth value of the convolved data
        ntmp++;
        }
    widthfac-=0.01;  // Reduce width factor by -0.1 and search again if no points were found out of transit & eclipse
    }
  rmsmin = RMS(tmpdobarry1,ntmp);  // RMS assuming gaussian noise and excluding pri and sec events  
  results.fred = sqrt(nintrans)*STDEV(tmpdobarry2,ntmp)/rmsmin;  // Have to account for the number of points in-transit we're averaging over.

  // Calculate Sigma_FA - the false alarm rate. I've set it so we shoudl only see one false alarm in 10,000 KOIs, assuming gaussian noise.
  results.sigfa1 = sqrt(2)*INVERFC((tend-tstart)/(period*20000));  // Assume 20,000 KOIs
  results.sigfa2 = sqrt(2)*INVERFC((tend-tstart)/(period));

  // If a sec, ter, or pos event was not found, set it to zero effectively
  if(chi2outlow == 9E9)
    chi2outlow = chi2med;
  if(chi2outhigh == -9E9)
    chi2outhigh = chi2med;
  if(chi2terlow == 9E9)
    chi2terlow = chi2med;

  // Calcualte the depths of the other events
  results.depfacsec = (chi2outlow-chi2med)/(chi2inlow-chi2med);
  results.depfacter = (chi2terlow-chi2med)/(chi2inlow-chi2med);
  results.depfacpos = (chi2outhigh-chi2med)/(chi2inlow-chi2med);
  
  // Significance of primary, secondary, tertiary, and positive events
  results.depsig = rmsmin/sqrt(nintrans);  // Noise level within transit duration assuming gaussian noise
  results.sigpri = (-results.tdepth*(chi2inlow-chi2med)/(chi2inlow-chi2med))/results.depsig;
  results.sigsec = (-results.tdepth*(chi2outlow-chi2med)/(chi2inlow-chi2med))/results.depsig;
  results.sigter = (-results.tdepth*(chi2terlow-chi2med)/(chi2inlow-chi2med))/results.depsig;
  results.sigpos = (-results.tdepth*((2*chi2med-chi2outhigh)-chi2med)/(chi2inlow-chi2med))/results.depsig;

  // If significance is really low just call it zero.
  if(fabs(results.sigpri)<1E-2) results.sigpri=0;
  if(fabs(results.sigsec)<1E-2) results.sigsec=0;
  if(fabs(results.sigter)<1E-2) results.sigter=0;
  if(fabs(results.sigpos)<1E-2) results.sigpos=0;

  // Make all phases fit within -0.25 to 0.75
  if(results.prilowtime<-0.25*period)
    results.prilowtime+=period;

  if(results.seclowtime<-0.25*period)
    results.seclowtime+=period;

  if(results.terlowtime<-0.25*period)
    results.terlowtime+=period;

  if(results.sechightime<-0.25*period)
    results.sechightime+=period; 
    
  // Output file for plotting
  tmpstr1 = basename + "-outfile2.dat";
  outfile.open(tmpstr1.c_str());
  for(i=0;i<ndat;i++)
    {
    convolveddepth[i] = results.tdepth*(chi2[i]-chi2med)/(chi2inlow-chi2med);
    outfile << fixed << setprecision(10) << data[midti+i].phase/period << " " << rms[i] << " " << chi2[i] << " " << convolveddepth[i] << endl;
    }
  outfile.close();
  }

//////////////////////////////////////////////////////////////////////////////////

void PLOT()
  {
  // Sort for plotting
  syscmd = "sort -k 1n " + basename + "-outfile2.dat > " + basename + "-outfile3.dat";
  system(syscmd.c_str());

  // Make files for binning
  syscmd = "awk '{print($1,$2)}' " + basename + "-outfile1.dat > " + basename +"-bininput.dat";
  system(syscmd.c_str());

  syscmd = "awk '{print($1,$2)}' " +  basename + "-outfile1-odd.dat > " + basename +"-bininput-odd.dat";
  system(syscmd.c_str());

  syscmd = "awk '{print($1,$2)}' " + basename + "-outfile1-evn.dat > " + basename +"-bininput-evn.dat";
  system(syscmd.c_str());

  // Bin data for the binned data
  BIN_NOERR(basename+"-bininput.dat",((tend-tstart)/period)/10,basename+"-binned2.dat");
  BIN_NOERR(basename+"-bininput-odd.dat",((tend-tstart)/period)/10,basename+"-binned2-odd.dat");
  BIN_NOERR(basename+"-bininput-evn.dat",((tend-tstart)/period)/10,basename+"-binned2-evn.dat");
  
  // Make plot file
  tmpstr1 = basename + "-ModShift.gnu";
  outfile.open(tmpstr1.c_str());
  outfile << "set term pdfcairo enhanced dashed size 8.5in,11in font ',16'" << endl; //pdfcairo enhanced size 5.5in,4.25in" << endl;
  outfile << "set output '" << basename << "-modshift.pdf'" << endl;
  outfile << "set multiplot layout 4,1 title \"" << objectname << ", P = " << setprecision(6) << period << " Days, E = " << setprecision(6) << epoch << " Days\" font ',20'" << endl;

  // Make the table containing the results
  outfile << "a = 0.045" << endl;
  outfile << "b = 0.0575" << endl;
  outfile << "c = a+7.25*b" << endl;
  outfile << "d = 0.08" << endl;
  outfile << "e = 0.955" << endl;
  
  outfile << "set arrow from screen (a-0.5*b),(e+0.01)       to screen (c+6.5*d),(e+0.01)       nohead lt 1 lw 5 lc 7" << endl;  // Top outside line
  outfile << "set arrow from screen (a-0.5*b),(e-0.025-0.01) to screen (c+6.5*d),(e-0.025-0.01) nohead lt 1 lw 5 lc 7" << endl;  // Bottom outside line
  outfile << "set arrow from screen (a-0.5*b),(e-0.0125) to screen (c+6.5*d),(e-0.0125) nohead lt 1 lw 5 lc 7" << endl;  // Middle line
  outfile << "set arrow from screen (a-0.5*b),(e+0.01) to screen (a-0.5*b),(e-0.025-0.01) nohead lt 1 lw 5 lc 7" << endl;  // Left end
  outfile << "set arrow from screen (c+6.5*d),(e+0.01) to screen (c+6.5*d),(e-0.025-0.01) nohead lt 1 lw 5 lc 7" << endl;  // Right end

  outfile << "set arrow from screen (a+0.5*b),(e-0.025-0.01) to screen (a+0.5*b),(e+0.01) nohead lt 1 lw 5 lc 7" << endl;
  outfile << "set arrow from screen (a+1.5*b),(e-0.025-0.01) to screen (a+1.5*b),(e+0.01) nohead lt 1 lw 5 lc 7" << endl;
  outfile << "set arrow from screen (a+2.5*b),(e-0.025-0.01) to screen (a+2.5*b),(e+0.01) nohead lt 1 lw 5 lc 7" << endl;
  outfile << "set arrow from screen (a+3.5*b),(e-0.025-0.01) to screen (a+3.5*b),(e+0.01) nohead lt 1 lw 5 lc 7" << endl;
  outfile << "set arrow from screen (a+4.5*b),(e-0.025-0.01) to screen (a+4.5*b),(e+0.01) nohead lt 1 lw 5 lc 7" << endl;
  outfile << "set arrow from screen (a+5.5*b),(e-0.025-0.01) to screen (a+5.5*b),(e+0.01) nohead lt 1 lw 5 lc 7" << endl;
  outfile << "set arrow from screen (a+6.5*b),(e-0.025-0.01) to screen (a+6.5*b),(e+0.01) nohead lt 1 lw 5 lc 7" << endl;
  
  outfile << "set arrow from screen (c+0.5*d),(e-0.025-0.01) to screen (c+0.5*d),(e+0.01) nohead lt 1 lw 5 lc 7" << endl;
  outfile << "set arrow from screen (c+1.5*d),(e-0.025-0.01) to screen (c+1.5*d),(e+0.01) nohead lt 1 lw 5 lc 7" << endl;
  outfile << "set arrow from screen (c+2.5*d),(e-0.025-0.01) to screen (c+2.5*d),(e+0.01) nohead lt 1 lw 5 lc 7" << endl;
  outfile << "set arrow from screen (c+3.5*d),(e-0.025-0.01) to screen (c+3.5*d),(e+0.01) nohead lt 1 lw 5 lc 7" << endl;
  outfile << "set arrow from screen (c+4.5*d),(e-0.025-0.01) to screen (c+4.5*d),(e+0.01) nohead lt 1 lw 5 lc 7" << endl;
  outfile << "set arrow from screen (c+5.5*d),(e-0.025-0.01) to screen (c+5.5*d),(e+0.01) nohead lt 1 lw 5 lc 7" << endl;


  // If I want to put Sigma symbol back, it's {/Symbol s}_
  outfile << "set label \"{Pri}\"         at screen (a+0*b),e center font ', 16'" << endl;
  outfile << "set label \"{Sec}\"         at screen (a+1*b),e center font ', 16'" << endl;
  outfile << "set label \"{Ter}\"         at screen (a+2*b),e center font ', 16'" << endl;
  outfile << "set label \"{Pos}\"         at screen (a+3*b),e center font ', 16'" << endl;
  outfile << "set label \"{FA_{1}}\"          at screen (a+4*b),e center font ', 16'" << endl;
  outfile << "set label \"{FA_{2}}\"         at screen (a+5*b),e center font ', 16'" << endl;
  outfile << "set label \"F_{Red}\"                   at screen (a+6*b),e center font ', 16'" << endl;
//   outfile << "set label \"{Pri}/F_{red}\" at screen (c+0*d),e center font ', 16'" << endl;
  outfile << "set label \"{Pri-Ter}\"     at screen (c+0*d),e center font ', 16'" << endl;
  outfile << "set label \"{Pri-Pos}\"     at screen (c+1*d),e center font ', 16'" << endl;
//   outfile << "set label \"{Sec}/F_{red}\" at screen (c+2*d),e center font ', 16'" << endl;
  outfile << "set label \"{Sec-Ter}\"     at screen (c+2*d),e center font ', 16'" << endl;
  outfile << "set label \"{Sec-Pos}\"     at screen (c+3*d),e center font ', 16'" << endl;
  outfile << "set label \"{Odd-Evn}\"    at screen (c+4*d),e center font ', 16'" << endl;
  outfile << "set label \"{DMM}\"    at screen (c+5*d),e center font ', 16'" << endl;
  outfile << "set label \"{Shape}\"    at screen (c+6*d),e center font ', 16'" << endl;
  
  outfile << "set label \"" << FORMAT(results.sigpri)                        << "\" at screen (a+0*b),(e-0.025) center font ',16'";
  if(results.sigpri/results.fred<results.sigfa1)        outfile << " textcolor lt 1" << endl; else outfile << " textcolor lt 7" << endl;
  
  outfile << "set label \"" <<  FORMAT(results.sigsec)                       << "\" at screen (a+1*b),(e-0.025) center font ',16'";
  if(results.sigsec/results.fred>results.sigfa1)        outfile << " textcolor lt 1" << endl; else outfile << " textcolor lt 7" << endl;
  
  outfile << "set label \"" <<  FORMAT(results.sigter)                       << "\" at screen (a+2*b),(e-0.025) center font ',16'";
  if(results.sigter/results.fred>results.sigfa1)        outfile << " textcolor lt 1" << endl; else outfile << " textcolor lt 7" << endl;
  
  outfile << "set label \"" <<  FORMAT(results.sigpos)                       << "\" at screen (a+3*b),(e-0.025) center font ',16'";
  if(results.sigpos/results.fred>results.sigfa1)        outfile << " textcolor lt 1" << endl; else outfile << " textcolor lt 7" << endl;
    
  outfile << "set label \"" <<  FORMAT(results.sigfa1)                       << "\" at screen (a+4*b),(e-0.025) center font ',16' textcolor lt 7" << endl;
  
  outfile << "set label \"" <<  FORMAT(results.sigfa2)                       << "\" at screen (a+5*b),(e-0.025) center font ',16' textcolor lt 7" << endl;
  
  outfile << "set label \"" <<  FORMAT(results.fred)                         << "\" at screen (a+6*b),(e-0.025) center font ',16' textcolor lt 7" << endl;

//   outfile << "set label \"" <<  FORMAT(results.sigpri/results.fred)          << "\" at screen (c+0*d),(e-0.025) center font ',16'";
//   if(results.sigpri/results.fred<results.sigfa1)   outfile << " textcolor lt 1" << endl; else outfile << " textcolor lt 7" << endl;
  
  outfile << "set label \"" <<  FORMAT(results.sigpri-results.sigter)        << "\" at screen (c+0*d),(e-0.025) center font ',16'";
  if(results.sigpri-results.sigter<results.sigfa2) outfile << " textcolor lt 1" << endl; else outfile << " textcolor lt 7" << endl;
  
  outfile << "set label \"" <<  FORMAT(results.sigpri-results.sigpos)        << "\" at screen (c+1*d),(e-0.025) center font ',16'";
  if(results.sigpri-results.sigpos<results.sigfa2) outfile << " textcolor lt 1" << endl; else outfile << " textcolor lt 7" << endl;
  
//   outfile << "set label \"" <<  FORMAT(results.sigsec/results.fred)          << "\" at screen (c+2*d),(e-0.025) center font ',16'";
//   if(results.sigsec/results.fred>results.sigfa1)   outfile << " textcolor lt 1" << endl; else outfile << " textcolor lt 7" << endl;
  
  outfile << "set label \"" <<  FORMAT(results.sigsec-results.sigter)        << "\" at screen (c+2*d),(e-0.025) center font ',16'";
  if(results.sigsec-results.sigter>results.sigfa2) outfile << " textcolor lt 1" << endl; else outfile << " textcolor lt 7" << endl;
  
  outfile << "set label \"" <<  FORMAT(results.sigsec-results.sigpos)        << "\" at screen (c+3*d),(e-0.025) center font ',16'";
  if(results.sigsec-results.sigpos>results.sigfa2) outfile << " textcolor lt 1" << endl; else outfile << " textcolor lt 7" << endl;

  outfile << "set label \"" <<  FORMAT(results.sigoe) << "\" at screen (c+4*d),(e-0.025) center font ',16'";
  if(results.sigoe>results.sigfa1) outfile << " textcolor lt 1" << endl; else outfile << " textcolor lt 7" << endl;

  outfile << "set label \"" <<  FORMAT(results.depthmedmeanrat) << "\" at screen (c+5*d),(e-0.025) center font ',16'";
  if(results.depthmedmeanrat>1.5) outfile << " textcolor lt 1" << endl; else outfile << " textcolor lt 7" << endl;
  
  outfile << "set label \"" <<  FORMAT(results.shape) << "\" at screen (c+6*d),(e-0.025) center font ',16'";
  if(results.shape>0.3) outfile << " textcolor lt 1" << endl; else outfile << " textcolor lt 7" << endl;
  
  // First plot
  outfile << "set origin 0.0,0.67" << endl;
  outfile << "set xlabel 'Phase' offset 0,0.4" << endl;
  outfile << "set xrange [-0.25 to 1.25]" << endl;
  outfile << "set x2range [-0.25 to 1.25]" << endl;
  outfile << "set xtics 0.25" << endl;
  outfile << "set format x '%4.2f'" << endl;
  outfile << "set format y '%6.0f'" << endl;
  outfile << "set ylabel 'Flux (ppm)'" << endl;
  outfile << "set object rect from 0.75,-1E7 to 1.25,1E7 fc rgb '#D3D3D3' lw 0" << endl;
  outfile << "set label '' at " << setprecision(10) << results.prilowtime/period  << ", graph 0.015 front point pt 9  ps 0.8" << endl;
  outfile << "set label '' at " << setprecision(10) << results.seclowtime/period  << ", graph 0.015 front point pt 9  ps 0.8" << endl;
  outfile << "set label '' at " << setprecision(10) << results.sechightime/period << ", graph 0.985 front point pt 11 ps 0.8" << endl;
  outfile << "set label '' at " << setprecision(10) << results.terlowtime/period  << ", graph 0.015 front point pt 8  ps 0.8" << endl;
  outfile << "stats '" << basename << "-binned2.dat' u 2 nooutput" << endl;
  // outfile << "set yrange []" << endl;
  outfile << "set yrange [1.0E6*STATS_min-0.5*(STATS_max-STATS_min) to 1.0E6*STATS_max+0.5*(STATS_max-STATS_min)]" << endl;
  outfile << "set autoscale y" << endl;
  outfile << "plot '" << basename << "-outfile1.dat' u 1:($2*1.0E6) pt 7 ps 0.1 lc 1 notitle, '' u ($1+1.0):($2*1.0E6) pt 7 ps 0.1 lc 1 notitle, '' u ($1+2.0):($2*1.0E6) pt 7 ps 0.1 lc 1 notitle, '" << basename << "-binned2.dat' u 1:($2*1.0E6) pt 7 ps 0.1 lc 3 notitle, '' u ($1+1.0):($2*1.0E6) pt 7 ps 0.1 lc 3 notitle, '" << basename << "-outfile1.dat' u 1:($3*1.0E6) with lines lt 1 lc 7 lw 5 notitle, '' u ($1+1.0):($3*1.0E6) with lines lt 1 lc 7 lw 5 notitle" << endl;

  // Second Plot
  outfile << "set origin 0.0,0.435" << endl;
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
  outfile << "set label '' at " << setprecision(10) << results.prilowtime/period  << ", graph 0.015 front point pt 9  ps 0.8" << endl;
  outfile << "set label '' at " << setprecision(10) << results.seclowtime/period  << ", graph 0.015 front point pt 9  ps 0.8" << endl;
  outfile << "set label '' at " << setprecision(10) << results.sechightime/period << ", graph 0.985 front point pt 11 ps 0.8" << endl;
  outfile << "set label '' at " << setprecision(10) << results.terlowtime/period  << ", graph 0.015 front point pt 8  ps 0.8" << endl;
  outfile << "set object rect from 0.75,-1E7 to 1.25,1E7 fc rgb '#D3D3D3' lw 0" << endl;
  outfile << "plot '" << basename << "-outfile3.dat' u 1:(1E6*$4) with lines lt 1 lc 7 notitle, '' u ($1+1.0):(1E6*$4) with lines lt 1 lc 7 notitle, 0 with lines lt 2 lc 1 lw 5 notitle, " << setprecision(10) << results.sigfa1*1E6*results.depsig << " with lines lt 2 lc 3 lw 5 notitle, " << -1.0*results.sigfa1*1.0E6*results.depsig << " with lines lt 2 lc 3 lw 5 notitle" << endl;

  outfile << "unset arrow" << endl;
  outfile << "unset object" << endl;
  outfile << "unset label" << endl;

  // Pri Zoom
  outfile << "set size square 0.375,0.275" << endl;
  outfile << "set origin 0.0,0.2" << endl;
  outfile << "set label 'Primary' at graph 0.5,0.925 center front" << endl;
  outfile << "set xrange [" << setprecision(10) << results.prilowtime/period+3*tstart/period << " to " << results.prilowtime/period+3*tend/period << "]" << endl;
  outfile << "set xtics " << setprecision(10) << (3*tend/period-3*tstart/period)/3.0 << endl;
  outfile << "set xtics format '%5.3f' mirror" << endl;
  outfile << "set xlabel ' '" << endl;
  outfile << "unset x2tics" << endl;
  outfile << "unset x2label" << endl;
  outfile << "unset y2tics" << endl;
  outfile << "unset y2label" << endl;
  outfile << "set autoscale y" << endl;
  outfile << "set ytics format '%6.0f' mirror" << endl;
  outfile << "set ytics auto" << endl;
  outfile << "set ylabel 'Flux (ppm)'" << endl;
  outfile << "plot '" << basename << "-binned2.dat' u 1:($2*1E6):($3*1E6) with yerrorbars lt 1 pt 7 ps 0.5 lc 3 notitle, '' u ($1+1.0):($2*1E6):($3*1E6) with yerrorbars lt 1 pt 7 ps 0.5 lc 3 notitle, '' u ($1-1.0):($2*1E6):($3*1E6) with yerrorbars lt 1 pt 7 ps 0.5 lc 3 notitle, '" << basename << "-outfile1.dat' u 1:($3*1E6) with lines lt 1 lc 7 notitle, '' u ($1+1.0):($3*1E6) with lines lt 1 lc 7 notitle, '' u ($1-1.0):($3*1E6) with lines lt 1 lc 7 notitle" << endl;

  outfile << "unset label" << endl;

    // Odd Zoom
  outfile << "set size square 0.375,0.275" << endl;
  outfile << "set origin 0.315,0.2" << endl;

  outfile << "set label 'Odd' at graph 0.5,0.925 center front" << endl;
  outfile << "set xrange [" << setprecision(10) << results.prilowtime/period+3*tstart/period << " to " << results.prilowtime/period+3*tend/period << "]" << endl;
  outfile << "set xtics " << setprecision(10) << (3*tend/period-3*tstart/period)/3.0 << endl;
  outfile << "set xtics format '%5.3f' mirror" << endl;
  outfile << "set xlabel ' '" << endl;
  outfile << "unset x2tics" << endl;
  outfile << "unset x2label" << endl;
  outfile << "unset y2tics" << endl;
  outfile << "unset y2label" << endl;
  outfile << "unset autoscale y" << endl;
  outfile << "set ytics format '%6.0f' mirror" << endl;
  outfile << "set ytics auto" << endl;
  outfile << "set ylabel ' '" << endl;
  outfile << "plot '" << basename << "-binned2-odd.dat' u 1:($2*1E6):($3*1E6) with yerrorbars lt 1 pt 7 ps 0.5 lc 3 notitle, '' u ($1+1.0):($2*1E6):($3*1E6) with yerrorbars lt 1 pt 7 ps 0.5 lc 3 notitle, '' u ($1-1.0):($2*1E6):($3*1E6) with yerrorbars lt 1 pt 7 ps 0.5 lc 3 notitle, '" << basename << "-outfile1.dat' u 1:($3*1E6*" << results.odddepth/results.tdepth << ") with lines lt 1 lc 7 notitle, '' u ($1+1.0):($3*1E6*" << results.odddepth/results.tdepth << ") with lines lt 1 lc 7 notitle, '' u ($1-1.0):($3*1E6*" << results.odddepth/results.tdepth << ") with lines lt 1 lc 7 notitle" << endl;
  
  outfile << "unset label" << endl;
  
  // Even Zoom
  outfile << "set size square 0.375,0.275" << endl;
  outfile << "set origin 0.63,0.2" << endl;
  outfile << "set label 'Even' at graph 0.5,0.925 center front" << endl;
  outfile << "set xrange [" << setprecision(10) << results.prilowtime/period+3*tstart/period << " to " << results.prilowtime/period+3*tend/period << "]" << endl;
  outfile << "set xtics " << setprecision(10) << (3*tend/period-3*tstart/period)/3.0 << endl;
  outfile << "set xtics format '%5.3f' mirror" << endl;
  outfile << "set xlabel ' '" << endl;
  outfile << "unset x2tics" << endl;
  outfile << "unset x2label" << endl;
  outfile << "unset y2tics" << endl;
  outfile << "unset y2label" << endl;
//   outfile << "set autoscale y" << endl;
  outfile << "set ytics format '%6.0f' mirror" << endl;
  outfile << "set ytics auto" << endl;
  outfile << "set ylabel ' '" << endl;
  outfile << "plot '" << basename << "-binned2-evn.dat' u 1:($2*1E6):($3*1E6) with yerrorbars lt 1 pt 7 ps 0.5 lc 3 notitle, '' u ($1+1.0):($2*1E6):($3*1E6) with yerrorbars lt 1 pt 7 ps 0.5 lc 3 notitle, '' u ($1-1.0):($2*1E6):($3*1E6) with yerrorbars lt 1 pt 7 ps 0.5 lc 3 notitle, '" << basename << "-outfile1.dat' u 1:($3*1E6*" << results.evndepth/results.tdepth << ") with lines lt 1 lc 7 notitle, '' u ($1+1.0):($3*1E6*" << results.evndepth/results.tdepth << ") with lines lt 1 lc 7 notitle, '' u ($1-1.0):($3*1E6*" << results.evndepth/results.tdepth << ") with lines lt 1 lc 7 notitle" << endl;

  outfile << "unset label" << endl;
  
  // Sec Zoom
  outfile << "set size square 0.375,0.275" << endl;
  outfile << "set origin 0.0, -0.015" << endl;
  outfile << "set label 'Secondary' at graph 0.5,0.925 center front" << endl;
  outfile << "set xrange [" << setprecision(10) << results.seclowtime/period+3*tstart/period << " to " << results.seclowtime/period+3*tend/period << "]" << endl;
  outfile << "set xtics " << setprecision(10) << (3*tend/period-3*tstart/period)/3.0 << endl;
  outfile << "set xtics format '%5.3f' mirror" << endl;
  outfile << "set xlabel 'Phase'" << endl;
  outfile << "unset x2tics" << endl;
  outfile << "unset x2label" << endl;
  outfile << "unset y2tics" << endl;
  outfile << "unset y2label" << endl;
  outfile << "set autoscale y" << endl;
  outfile << "set ytics format '%6.0f' mirror" << endl;
  outfile << "set ytics auto" << endl;
  outfile << "unset ylabel" << endl;
  outfile << "set ylabel 'Flux (ppm)'" << endl;
  outfile << "plot '" << basename << "-binned2.dat' u 1:($2*1E6):($3*1E6) with yerrorbars lt 1 pt 7 ps 0.5 lc 3 notitle, '' u ($1+1.0):($2*1E6):($3*1E6) with yerrorbars lt 1 pt 7 ps 0.5 lc 3 notitle, '' u ($1-1.0):($2*1E6):($3*1E6) with yerrorbars lt 1 pt 7 ps 0.5 lc 3 notitle, '" << basename << "-outfile1.dat' u ($1+" << setprecision(10) << results.seclowtime/period << "):((" << results.baseflux << " + " << results.depfacsec << "*($3-" << results.baseflux << "))*1E6) with lines lt 1 lc 7 notitle, '' u ($1+1.0+" << results.seclowtime/period << "):((" << results.baseflux << " + " << results.depfacsec << "*($3-" << results.baseflux << "))*1E6) with lines lt 1 lc 7 notitle, '' u ($1-1.0+" << results.seclowtime/period << "):((" << results.baseflux << " + " << results.depfacsec << "*($3-" << results.baseflux << "))*1E6) with lines lt 1 lc 7 notitle" << endl;

  outfile << "unset label" << endl;


  // Ter Zoom
  if(results.terlowtime/period>-0.25)
    {
    outfile << "set size square 0.375,0.275" << endl;
    outfile << "set origin 0.315, -0.015" << endl;
    outfile << "set label 'Tertiary' at graph 0.5,0.925 center front" << endl;
    outfile << "set xtics " << setprecision(10) << (3*tend/period-3*tstart/period)/3.0 << endl;
    outfile << "set xtics format '%5.3f' mirror" << endl;
    outfile << "set xrange [" << setprecision(10) << results.terlowtime/period+3*tstart/period << " to " << results.terlowtime/period+3*tend/period << "]" << endl;
    outfile << "set xlabel 'Phase'" << endl;
    outfile << "unset x2tics" << endl;
    outfile << "unset x2label" << endl;
    outfile << "unset y2tics" << endl;
    outfile << "unset y2label" << endl;
    outfile << "set autoscale y" << endl;
    outfile << "set ytics format '%6.0f' mirror" << endl;
    outfile << "set ytics auto" << endl;
    outfile << "set ylabel ' '" << endl;
    outfile << "plot '" << basename << "-binned2.dat' u 1:($2*1E6):($3*1E6) with yerrorbars lt 1 pt 7 ps 0.5 lc 3 notitle, '' u ($1+1.0):($2*1E6):($3*1E6) with yerrorbars lt 1 pt 7 ps 0.5 lc 3 notitle, '' u ($1-1.0):($2*1E6):($3*1E6) with yerrorbars lt 1 pt 7 ps 0.5 lc 3 notitle, '" << basename << "-outfile1.dat' u ($1+" << setprecision(10) << results.terlowtime/period << "):((" << results.baseflux << " + " << results.depfacter << "*($3-" << results.baseflux << "))*1E6) with lines lt 1 lc 7 notitle, '' u ($1+1.0+" << results.terlowtime/period << "):((" << results.baseflux << " + " << results.depfacter << "*($3-" << results.baseflux << "))*1E6) with lines lt 1 lc 7 notitle, '' u ($1-1.0+" << results.terlowtime/period << "):((" << results.baseflux << " + " << results.depfacter << "*($3-" << results.baseflux << "))*1E6) with lines lt 1 lc 7 notitle" << endl;
    }
  
  outfile << "unset label" << endl;

  // Pos Zoom
  if(results.sechightime/period>-0.25)
    {
    outfile << "set size square 0.375,0.375" << endl;
    outfile << "set origin 0.63, -0.065" << endl;
    outfile << "set label 'Positive' at graph 0.5,0.075 center front" << endl;
    outfile << "set xtics " << setprecision(10) << (3*tend/period-3*tstart/period)/3.0 << endl;
    outfile << "set xtics format '%5.3f' mirror" << endl;
    outfile << "set xrange [" << setprecision(10) << results.sechightime/period+3*tstart/period << " to " << results.sechightime/period+3*tend/period << "]" << endl;
    outfile << "set xlabel 'Phase'" << endl;
    outfile << "unset x2tics" << endl;
    outfile << "unset x2label" << endl;
    outfile << "unset y2tics" << endl;
    outfile << "unset y2label" << endl;
    outfile << "set autoscale y" << endl;
    outfile << "set ytics format '%6.0f' mirror" << endl;
    outfile << "set ytics auto" << endl;
    outfile << "set ylabel ' '" << endl;
    outfile << "plot '" << basename << "-binned2.dat' u 1:($2*1E6):($3*1E6) with yerrorbars lt 1 pt 7 ps 0.5 lc 3 notitle, '' u ($1+1.0):($2*1E6):($3*1E6) with yerrorbars lt 1 pt 7 ps 0.5 lc 3 notitle, '' u ($1-1.0):($2*1E6):($3*1E6) with yerrorbars lt 1 pt 7 ps 0.5 lc 3 notitle, '" << basename << "-outfile1.dat' u ($1+" << setprecision(10) << results.sechightime/period << "):((" << results.baseflux << " + " << results.depfacpos << "*($3-" << results.baseflux << "))*1E6) with lines lt 1 lc 7 notitle, '' u ($1+1.0+" << results.sechightime/period << "):((" << results.baseflux << " + " << results.depfacpos << "*($3-" << results.baseflux << "))*1E6) with lines lt 1 lc 7 notitle, '' u ($1-1.0+" << results.sechightime/period << "):((" << results.baseflux << " + " << results.depfacpos << "*($3-" << results.baseflux << "))*1E6) with lines lt 1 lc 7 notitle" << endl;
    }
  
  outfile << "unset label" << endl;
  
  
  // And make the whole thing
  outfile << "unset multiplot" << endl;
  tmpstr1 = "gnuplot " + basename + "-ModShift.gnu";
  if(modplotint==1)
    system(tmpstr1.c_str());


  // Clean up files
  // syscmd = "cp outfile1-" + basename + ".dat " + basename + ".cln";  // Make clean file for use later - COMMENTING OUT - NOT NEEDED AT PRESENT
  // system(syscmd.c_str());
//  syscmd = "rm " + basename + "-outfile*.dat";
//  system(syscmd.c_str());
//  syscmd = "rm " + basename + "-bin*dat";
//  system(syscmd.c_str());
//  syscmd = "rm " + basename + "-ModShift.gnu";
//  system(syscmd.c_str());
  }
  
//////////////////////////////////////////////////////////////////////////////

string FORMAT(double f)  // MAKES IT FIT WITHIN 5 DIGITS - CAN RE-WRITE IN FUTURE TO MAKE MORE GENERIC
  {
  stringstream ss;
  string outstring;
  
  if(f==0)
    return "0";
  
  int d = (int)::ceil(::log10(f < 0 ? -f : f)); /*digits before decimal point*/
  
  if(d<2)
    {
    ss << fixed << setprecision(2) << f;
    outstring = ss.str();
    }
  
  if(d>=2 && d<4)
    {
    ss << fixed << setprecision(1) << f;
    outstring = ss.str();
    }
    
  if(d>=4 && d<6)
    {
    ss << fixed << setprecision(0) << f;
    outstring = ss.str();
    }
    
  if(d>=6)
    {
    ss << setiosflags(ios::uppercase) << scientific << setprecision(1) << f;
    outstring = ss.str();
    size_t f = outstring.find("E+0");
    outstring.replace(f,string("E+0").length(),"E");
    }
    
  return outstring;
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
  
  sort(w,w+n,double_sort);
  
  if(n%2==0)
    return (w[n/2]+w[n/2-1])/2;  // Find median
  else
    return w[n/2];
  }
  
///////////////////////////////////////////////////////////////////////////////

double MAX(double x[], int n)
  {
  double max=-9E99;
  int z;
  
  for(z=0;z<n;z++)
    if(x[z]>max)
      max=x[z];
    
  return max;
  }


///////////////////////////////////////////////////////////////////////////////

double MIN(double x[], int n)
  {
  double min=9E99;
  int z;
  
  for(z=0;z<n;z++)
    if(x[z]<min)
      min=x[z];
    
  return min;
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

double MEAN(double x[], int n) {  // Return mean of first n terms

int z;
double xmean;

xmean=0.0;
for(z=0;z<n;z++)
  xmean += x[z];

return(xmean/n);
}
  

////////////////////////////////////////////////////////////////////////////////////////////////////
  
void BIN_NOERR(string bininfile, double binsize, string binoutfile)
  {
  // Bin a file in that has 2 columns of Xvalue  Yvalue
  // Uses error-weighted mean and constant bin steps
  // Outputs new binned Xvalue, Yvalue, and error (from stdev) of the binned data point
  
  int i;  // Counting integers
  int n;  // n is number of points going into current bin - sw1 is a switch
  int Nraw,Nbin; // Number of original and binned data points
  double curbin;
  ifstream datain;
  ofstream dataout;
  
  
  if(binoutfile==bininfile)
    {
    cout << "Binned output file cannot be the same as input file. Exiting." << endl;
    exit(1);
    }

    
  i=0;
  datain.open(bininfile.c_str());
  datain >> inpdata[i].phase;
  while(!datain.eof())
    {
    datain >> inpdata[i].flux;
    i++;
    datain >> inpdata[i].phase;
    }
  datain.close();
  Nraw = i;

  
  Nbin=0;
  for(curbin=inpdata[0].phase; curbin<inpdata[Nraw-1].phase+binsize; curbin+=binsize)
    {
    n=0;
    sw1=0;
    for(i=0;i<Nraw;i++)
      if(inpdata[i].phase>=curbin && inpdata[i].phase<curbin+binsize)
        {
        tmpdob3[n] = inpdata[i].flux;
        n++;
        sw1=1;
        }

    if(sw1==1)
      {
      bindat[Nbin].phase = curbin+0.5*binsize;
      bindat[Nbin].flux = MEAN(tmpdob3,n);
      bindat[Nbin].resid = STDEV(tmpdob3,n)/sqrt(n);
      Nbin++;
      }
    }
  

  dataout.open(binoutfile.c_str());
  for(i=0;i<Nbin;i++)
    dataout << setprecision(10) << bindat[i].phase << " " << bindat[i].flux << " " << bindat[2].resid << endl;
  dataout.close();
  
  }

