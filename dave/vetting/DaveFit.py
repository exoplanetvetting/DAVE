# TO-DO: Modify modshift-dave to find the widths of primary and secondary. Can use that to estimate esinw.
# Make it so first fit from clipped data just use a JOIN to extend to entire data - baseline is flat.  Nice and quick.
# Odd-Even implement. Can just use first clipped data


import subprocess
import os
import sys
import math
from subprocess import check_output
import time

# Record start time of program
program_start_time = time.time()

# Constants
PI = 3.14159265359


lcname = sys.argv[1]
basename = sys.argv[2]
user_period = sys.argv[3]
user_epoch = sys.argv[4]
user_depth = sys.argv[5]  # Fractional Depth
user_duration = sys.argv[6]
user_teff = sys.argv[7]
user_logg = sys.argv[8]
user_metal = sys.argv[9]


#global period
period = float(user_period)
epoch = float(user_epoch)
depth = float(user_depth)
duration = float(user_duration)
teff = float(user_teff)
logg = float(user_logg)
metal = float(user_metal)
metal = round(2*metal)/2.0
if teff < 3500:
  teff = 3500
if teff > 40000:
  teff = 40000
if logg < 0.0:
    logg = 0.0
if logg > 5.0:
    logg = 5.0
if metal < -5.0:
    metal = -5.0
if metal > 1.0:
    metal = 1.0

# File names
paroutfile     = basename + ".par"
lcoutfile      = basename + ".out"
modeloutfile   = basename + ".fit"
fitfile        = basename + ".in"
clippedlcname  = basename + ".clp"
maglcname      = basename + ".mag"
modshiftinfile = basename + ".mod"
cleanlcname    = basename + ".cln"




def RUN_JKTEBOP(task,nint,timefluxfile):
    with open(fitfile, "wt") as infile:
        infile.write("%i %f Task to do (from 1 to 9)   Integ. ring size (deg)\n" % (task, intring))
        infile.write("%f %f Sum of the radii           Ratio of the radii\n" % (radsum, radrat))
        infile.write("%f %f Orbital inclination (deg)  Mass ratio of system\n" % (incl, massrat))
        infile.write("%f %f Orbital eccentricity       Periastron longitude deg\n" % (ecosw, esinw))
        infile.write("%f %f Gravity darkening (star A) Grav darkening (star B)\n" % (gravda, gravdb))
        infile.write("%f %f Surface brightness ratio   Amount of third light\n" % (surfbr, tlight))
        infile.write("%s %s LD law type for star A     LD law type for star B\n" % (ldlawa, ldlawb))
        
       #infile.write("%f %f LD star A (linear coeff)   LD star B (linear coeff)\n" % (lda1, ldb1))
       #infile.write("%f %f LD star A (nonlin coeff)   LD star B (nonlin coeff)\n" % (lda2, ldb2))
        
        infile.write("%f %f LD star A (coefficient 1) LD star B (coefficient 1)\n" % (lda1, ldb1))
        infile.write("%f %f LD star A (coefficient 2) LD star B (coefficient 2)\n" % (lda2, ldb2))
        infile.write("%f %f LD star A (coefficient 3) LD star B (coefficient 3)\n" % (lda3, ldb3))
        infile.write("%f %f LD star A (coefficient 4) LD star B (coefficient 4)\n" % (lda4, ldb4))

        infile.write("%f %f Reflection effect star A   Reflection effect star B\n" % (refla, reflb))
        infile.write("%f %f Phase shift of primary min Light scale factor (mag)\n" % (pshift, lscale))
        infile.write("%.10f    Orbital period of eclipsing binary system (days)\n" % period)
        infile.write("%.10f    Reference time of primary minimum (HJD)\n" % epoch)
        if task==4:
            infile.write("%.10f Sigma value to reject discrepant observations\n" % 5.0)
        infile.write("%i %i Adjust RADII SUM  or  RADII RATIO     (0 or 1 or 2)\n" % (fit_radsum, fit_radrat))
        infile.write("%i %i Adjust INCLINATION  or  MASSRATIO     (0 or 1 or 2)\n" % (fit_incl, fit_massrat))
        infile.write("%i %i Adjust ECCENTRICITY  or  OMEGA        (0 or 1 or 2)\n" % (fit_ecosw, fit_esinw))
        infile.write("%i %i Adjust GRAVDARK1  or  GRAVDARK2       (0 or 1 or 2)\n" % (fit_gravda, fit_gravdb))
        infile.write("%i %i Adjust SURFACEBRIGHT2  or  THIRDLIGHT (0 or 1 or 2)\n" % (fit_surfbr, fit_tlight))
        
       #infile.write("%i %i Adjust LD-lin1  or  LD-lin2           (0 or 1 or 2)\n" % (fit_lda1, fit_ldb1))
       #infile.write("%i %i Adjust LD-nonlin1  or  LD-nonlin2     (0 or 1 or 2)\n" % (fit_lda2, fit_ldb2))
        
        infile.write("%i %i Adjust LDcoeff-A1 or LDcoeff-B1        (0, 1, 2, 3)\n" % (fit_lda1, fit_ldb1))
        infile.write("%i %i Adjust LDcoeff-A2 or LDcoeff-B2        (0, 1, 2, 3)\n" % (fit_lda2, fit_ldb2))
        infile.write("%i %i Adjust LDcoeff-A3 or LDcoeff-B3        (0, 1, 2, 3)\n" % (fit_lda3, fit_ldb3))
        infile.write("%i %i Adjust LDcoeff-A4 or LDcoeff-B4        (0, 1, 2, 3)\n" % (fit_lda4, fit_ldb4))

        
        infile.write("%i %i Adjust REFLECTION COEFFS 1 and 2      (-1, 0, 1 ,2)\n" % (fit_refla, fit_reflb))
        infile.write("%i %i Adjust PHASESHIFT  or  SCALE FACTOR   (0 or 1 or 2)\n" % (fit_pshift, fit_lscale))
        infile.write("%i %i Adjust PERIOD  or  TZERO (min light)  (0 or 1)\n" % (fit_period, fit_epoch))
        infile.write("%s Name of file containing light curve\n" % timefluxfile)
        infile.write("%s Name of output parameter file\n" % paroutfile)
        infile.write("%s Name of output light curve file\n" % lcoutfile)
        infile.write("%s Name of output model light curve fit file\n" % modeloutfile)
        if nint > 0:
            infile.write("NUMI " + str(nint) + " 1764\n")

    # Delete old fit files, if they exist
    try:
        os.remove(paroutfile)
    except OSError:
        pass

    try:
        os.remove(lcoutfile)
    except OSError:
        pass
      
    try:
        os.remove(modeloutfile)
    except OSError:
        pass
      


    jkt_start_time = time.time()
    #return check_output(["./jktebop",fitfile])
    subprocess.call(["./jktebop",fitfile])
    jkt_end_time = time.time()
    print("--- %s second jktebop runtime ---" % (jkt_end_time - jkt_start_time))


def READ_JKTEBOP():
    global surfbr
    global radsum
    global radrat
    global lda1
    global lda2
    global lda3
    global lda4
    global ldb1
    global ldb2
    global ldb3
    global ldb4
    global incl
    global ecosw
    global esinw
    global gravda
    global gravdb
    global refla
    global reflb
    global massrat
    global tlight
    global pshift
    global lscale
    global period
    global epoch
  
    sw1=0
    for line in open(paroutfile,'r'):
        if line.find("Final") > -1:
            sw1=1
        if(sw1>0):
            if line.find("Surf. bright. ratio") > -1:
                info = line.split()
                surfbr = float(info[4])
                
            if line.find("Sum of frac radii") > -1:
                info = line.split()
                radsum = float(info[5])
                
            if line.find("Ratio of the radii") > -1:
                info = line.split()
                radrat = float(info[5])
                
            if line.find("Frac radius star A") > -1:
                info = line.split()
                radsum = float(info[5])
                
            if line.find("Frac radius star B") > -1:
                info = line.split()
                radrat = float(info[5])
                
            if line.find("Limb darkening A1") > -1:
                info = line.split()
                lda1 = float(info[4])
                
            if line.find("Limb darkening A2") > -1:
                info = line.split()
                lda2 = float(info[4])

            if line.find("Limb darkening A3") > -1:
                info = line.split()
                lda3 = float(info[4])
                
            if line.find("Limb darkening A4") > -1:
                info = line.split()
                lda4 = float(info[4])                
                
            if line.find("Limb darkening B1") > -1:
                info = line.split()
                ldb1 = float(info[4])
                
            if line.find("Limb darkening B2") > -1:
                info = line.split()
                ldb2 = float(info[4])
                
            if line.find("Limb darkening B3") > -1:
                info = line.split()
                ldb3 = float(info[4])
                
            if line.find("Limb darkening B4") > -1:
                info = line.split()
                ldb4 = float(info[4])                
                
            if line.find("Orbit inclination") > -1:
                info = line.split()
                incl = float(info[3])
                
            if line.find("ecc * cos(omega)") > -1:
                info = line.split()
                ecosw = float(info[4])
                
            if line.find("ecc * sin(omega)") > -1:
                info = line.split()
                esinw = float(info[4])
                
            if line.find("Grav darkening A") > -1:
                info = line.split()
                gravda = float(info[4])
                
            if line.find("Grav darkening B") > -1:
                info = line.split()
                gravdb = float(info[4])
                
            if line.find("Reflected light A") > -1:
                info = line.split()
                refla = float(info[4])
                
            if line.find("Reflected light B") > -1:
                info = line.split()
                reflb = float(info[4])
                
            if line.find("Phot mass ratio") > -1:
                info = line.split()
                massrat = float(info[4])
                
            if line.find("Third light (L_3)") > -1:
                info = line.split()
                tlight = float(info[4])
                
            if line.find("Phase correction") > -1:
                info = line.split()
                pshift = float(info[3])
                
            if line.find("Light scale factor") > -1:
                info = line.split()
                lscale = float(info[4])
                
            #if line.find("Integration ring") > -1:
                #info = line.split()
                #intring = float(info[3])             
                
            if line.find("Orbital period (P)") > -1:
                info = line.split()
                period = float(info[4])   
                
            if line.find("Ephemeris timebase") > -1:
                info = line.split()
                epoch = float(info[3])   
              
            # If we've read through all the values, as a fail safe stop parsing file
            if line.find("Limb darkening law for primary star") > -1:
                sw1 = 0
              
def PLOT(plotdata):
    sortedlcoutfile = basename + "-sort.out"
    subprocess.call("sort -k 4n %s > %s " % (plotdata,sortedlcoutfile),shell=True)


    plotfile = "plot.gnu"

    with open(plotfile, "wt") as infile:
        infile.write("set term pdfcairo enhanced dashed\n")
        infile.write("set output '%s.pdf'\n" % basename)
        infile.write("set multiplot layout 2, 1 title 'Kepler " + basename + "' font ',14'\n")
        infile.write("set xrange [-0.25 to 1.25]\n")
        infile.write("set yrange [] reverse\n")
        infile.write("set xlabel 'Phase'\n")
        infile.write("set ylabel 'Flux'\n")
        infile.write("set format y '%7.3f'\n")
        infile.write("plot '%s' u ($4-1.0):(10**($2/2.5)-1) pt 7 ps 0.2 lc 1 notitle, '' u ($4):(10**($2/2.5)-1) pt 7 ps 0.2 lc 1 notitle, '' u ($4+1.0):(10**($2/2.5)-1) pt 7 ps 0.2 lc 1 notitle, '' u ($4-1.0):(10**($5/2.5)-1) with lines lt 1 lw 2 lc 7 notitle, '' u ($4):(10**($5/2.5)-1) with lines lt 1 lw 2 lc 7 notitle,'' u ($4+1.0):(10**($5/2.5)-1) with lines lt 1 lw 2 lc 7 notitle\n" % sortedlcoutfile)
        infile.write("set xrange [%f to %f ]\n" % (-2.5*duration/period,2.5*duration/period))
        
        infile.write("set size 0.5,0.5\n")
        infile.write("set origin 0.0,0.0\n")
        infile.write("plot '%s' u ($4-1.0):(10**($2/2.5)-1) pt 7 ps 0.2 lc 1 notitle, '' u ($4):(10**($2/2.5)-1) pt 7 ps 0.2 lc 1 notitle, '' u ($4+1.0):(10**($2/2.5)-1) pt 7 ps 0.2 lc 1 notitle, '' u ($4-1.0):(10**($5/2.5)-1) with lines lt 1 lw 2 lc 7 notitle, '' u ($4):(10**($5/2.5)-1) with lines lt 1 lw 2 lc 7 notitle,'' u ($4+1.0):(10**($5/2.5)-1) with lines lt 1 lw 2 lc 7 notitle\n" % sortedlcoutfile)
        infile.write("unset multiplot\n")
      
    subprocess.call("gnuplot plot.gnu",shell=True)


def RUN_MODSHIFT():
    for line in open(lcoutfile,'r'):  # There's probably a better way of doing this rather than reading the entire file, but it works
        info = line.split()
    baseflux = info[4]

    # Get phase and flux for entire light curve. Using period/epoch from best fit
    script = "awk 'BEGIN{epoch=" + str("%.10f" % epoch) + "; period=" + str("%.10f" % period) + "} {if(NR==1) while(epoch>$1) epoch-=period; printf(\"%12.10f %11.8f\\n\",($1-epoch)/period - int(($1-epoch)/period),$2)}' " + maglcname + " | sort -k 1,1 > tmp1"
    subprocess.call(script,shell=True)

    # Get phase and model for clipped light curve.
    script = "awk 'BEGIN{epoch=" + str("%.10f" % epoch) + "; period=" + str("%.10f" % period) + "} {if(NR==2) while(epoch>$1) epoch-=period; if(NR>1) printf(\"%12.10f %11.8f\\n\",($1-epoch)/period - int(($1-epoch)/period),$5)}' " + lcoutfile + " | sort -k 1,1 > tmp2"
    subprocess.call(script,shell=True)
    
    # Combine them so that all data outside of clip gets baseline model flux
    script = "join -a 1 -o 1.1 1.2 2.2 -e " + baseflux + " tmp1 tmp2 > tmp3"
    subprocess.call(script,shell=True)
    
    # Now make input for model-shift plot, where the phase is in days, starting at 0.75.
    script = "awk '{if($1<0.75) phase = $1*"+str("%.10f" % period)+"; else phase = $1*"+str("%.10f" % period)+"-"+str("%.10f" % period)+"; printf(\"%13.10f %11.8f %11.8f\\n\",phase,-1*10**($2/2.5)+1,-1*10**($3/2.5)+1)}' tmp3 | sort -k 1n > " + modshiftinfile
    subprocess.call(script,shell=True)

    # Clean up tmp files
    os.remove("tmp1")
    os.remove("tmp2")
    os.remove("tmp3")
    
    # Run modshift, and return the output
    mod_start_time = time.time()
    return check_output(["./modshift",modshiftinfile,basename,str(period)])
    mod_end_time = time.time()
    print("--- %s second modshift runtime ---" % (mod_end_time - mod_start_time))


def GET_LDS(teff,logg,metal):
    return check_output(["./jktld",str(teff),str(logg),str(metal),"2","t","8","Kp"])
    # t is for Sing, 8 is for sing's table, Kp is for kepler bandpass


# Turn flux data file into magnitudes for JKTEBOP input
script = "awk '{print($1,-2.5*log(1.0+$2)/log(10))}' " + lcname + " > " + maglcname 
subprocess.call(script,shell=True)
#print script

# Estimate radrat and radsum based on depth
radrat = math.sqrt(depth)
radsum = (1.0+radrat)*(math.sin(PI*duration/period))/(1.0+math.sqrt(depth))   # Assumes b = 0
incl = (180.0/PI)*math.acos(0.5*(math.sin(PI*duration/period))/(1.0+math.sqrt(depth)))  # Assumes b = 0.5, which should be the average b value

# Get LD based on teff, logg, and metallicty
ldout = GET_LDS(teff,logg,metal)
info = ldout.split()
lda1 = ldb1 = 0.0  # Setting 0.0 for first claret term turns it into sing law
lda2 = ldb2 = float(info[0])
lda3 = ldb3 = float(info[1])
lda4 = ldb4 = float(info[2])



# Use some default initial values
task = 3
intring = 1.0
#incl = 89.95
massrat = -0.001
ecosw = 0.0
esinw = 0.0
gravda = 1.0
gravdb = 1.0
surfbr = 0.0
tlight = 0.0
ldlawa = "4par"
ldlawb = "4par"
#lda1 = 0.0
#lda2 = 0.3
#lda3 = 0.0
#lda4 = 0.0
#ldb1 = 0.0
#ldb2 = 0.3
#ldb3 = 0.0
#ldb4 = 0.0
refla = 0.0
reflb = 0.0
pshift = 0.0
lscale = 0.0
fit_radsum = 1
fit_radrat = 1
fit_incl = 1
fit_massrat = 0
fit_ecosw = 0
fit_esinw = 0
fit_gravda = 0
fit_gravdb = 0
fit_surfbr = 0
fit_tlight = 0
fit_lda1 = 0
fit_lda2 = 0
fit_lda3 = 0
fit_lda4 = 0
fit_ldb1 = 0
fit_ldb2 = 0
fit_ldb3 = 0
fit_ldb4 = 0
fit_refla = 0
fit_reflb = 0
fit_pshift = 0
fit_lscale = 1
fit_period = 1
fit_epoch = 1




# Only use data within 2.5 transit durations of time of central transit
script = "awk 'BEGIN{epoch=" + user_epoch + "; period=" + user_period + "; duration=" + user_duration + ";} {if(NR==1) while(epoch>$1) epoch-=period; if(($1-epoch)/period - int(($1-epoch)/period) < 2.5*duration/period || ($1-epoch)/period - int(($1-epoch)/period) > 1.0 - 2.5*duration/period) print($0)}' " + maglcname + " > " + clippedlcname
subprocess.call(script,shell=True)


# Run it the first time on the clipped data
#jktrunout = RUN_JKTEBOP(3,0,clippedlcname)
#print jktrunout
RUN_JKTEBOP(3,0,clippedlcname)

# Read in the fit parameters
READ_JKTEBOP()


# Scrapping next 6 commented lines I think - do modelshift first by filling in the baseline, then can start more complex fitting
# Run again, starting the fit params to the clipped data, now on the full time series
#RUN_JKTEBOP(0,maglcname)
# Re-read jkteop results
#READ_JKTEBOP()
# Make plots
#PLOT()



# Fill in baseline to do modshift



# Run modshift
domod = 0
if(domod==1):
    modshiftcmdout = RUN_MODSHIFT()

    # Read the modshift output back in to variables
    info = modshiftcmdout.split()
    mod_sig_pri = float(info[1])
    mod_sig_sec = float(info[2])
    mod_sig_ter = float(info[3])
    mod_sig_pos = float(info[4])
    mod_sig_fa = float(info[5])
    mod_Fred = float(info[6])
    mod_ph_pri = float(info[7])
    mod_ph_sec = float(info[8])
    mod_ph_ter = float(info[9])
    mod_ph_pos = float(info[10])
    mod_secdepth = float(info[11])
    mod_secdeptherr = float(info[12])




    # Now estaimte secondary params based on modshift output
    ecosw = (PI/2)*(mod_ph_sec - 0.5)  # http://www.astro.keele.ac.uk/jkt/pubs/Southworth-proc-PasDeDeux.pdf
    #ecosw = 0.0
    surfbr = mod_sig_sec/mod_sig_pri
    ldb2=0.3
    esinw = 0.001
    #esinw = (dur_sec - dur_pri)/(dur_sec + dur_pri)   # http://www.astro.keele.ac.uk/jkt/pubs/Southworth-proc-PasDeDeux.pdf   OR   http://arxiv.org/pdf/1201.1388.pdf
    fit_ecosw = 1
    fit_esinw = 1
    fit_surfbr = 1
    #fit_lda2 = 1
    #fit_ldb2 = 1

# And run JKTEBOP again, this time fitting for the secondary - I thought I could use the clean lc from the modshift test, but that's phased, and I don't want phased.
# May have to do my own clean. Mode 4 of jktebop will reject outliers, but we don't want to reject all systematics, and it takes a long time.
#jktrunout = RUN_JKTEBOP(4,0,maglcname)
#print jktrunout

#RUN_JKTEBOP(3,0,maglcname)


# Make plots
#PLOT(lcoutfile)


program_end_time = time.time()
print("--- %s second program runtime ---" % (program_end_time - program_start_time))





#print intring
#print incl
#print massrat
#print ecosw
#print esinw
#print gravda
#print gravdb
#print surfbr
#print tlight
#print lda1
#print ldb1
#print lda2
#print ldb2
#print refla
#print reflb
#print pshift
#print lscale
