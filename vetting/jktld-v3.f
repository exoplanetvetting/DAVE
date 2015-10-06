      PROGRAM JKTLD                    ! John Southworth  (jkt~astro.keele.ac.uk)
                                       ! Astrophysics Group, Keele University, UK
!================================================================================
! This program performs bilinear interpolation to determine limb darkening coeff-
! icients for the given effective temperature  (Teff)  and surface gravity (logg)
! from published grids of coefficients.   Depending on the number of arguments on
! the command line, JKTLD either creates an output file with a wide variety of LD
! coefficients or returns specific ones to standard output.
!================================================================================
! Linear limb darkening law:  I(mu)/I(0) = 1 - u*(1 - mu)         mu = cos(theta)
! Quadratic LD law:           I(mu)/I(0) = 1 - q1*(1 - mu) - q2*(1 - mu)^2
! Logarithmic LD law:         I(mu)/I(0) = 1 - l1*(1 - mu) - l2*mu*ln(mu)
! Square-root LD law:         I(mu)/I(0) = 1 - s1*(1 - mu) - s2*(1 - sqrt[mu])
! Exponential LD law:         I(mu)/I(0) = 1 - e1*(1 - mu) - e2/(1 - exp[mu])
! Sing three-parameter law:   I(mu)/I(0) = 1 - sum_(i=2..4) ti*[1 - mu^(i/2)]
! Claret four-parameter law:  I(mu)/I(0) = 1 - sum_(i=1..4) f_i*[1 - mu^(i/2)]
!================================================================================
! Source references for the limb darkening (LD) coefficients tabulations used:
! VanHam93a:  Van Hamme (1993AJ....106.2096V)
! DiazCo95a:  Diaz-Cordoves et al (1995A+AS..110..329D)
! Claret95a:  Claret et al (1995A+AS..114..247C)
! Claret00a:  Claret (2000A+A...363.1081C)
! Claret00p:  Claret (2000A+A...363.1081C)
! ClarHa03a:  Claret & Hauschildt (2003A+A...412..241C)
! Claret04a:  Claret (2004A+A...428.1001C)
! Claret04p:  Claret (2004A+A...428.1001C)
! Sing2010a:  Sing (2010A+A...510A..21S)
! DiazCo95a and Claret95a are complementary so are included together from now on.
!--------------------------------------------------------------------------------
! Model atmospheres used by each source.   In each case many combinations of Teff
! and logg are not available.  Similarly, coefficients for different [M/H] values
! are only available for the standard Vmicro=2km/s and coefficients for different
! Vmicro values are only available for the standard [M/H]=0.0dex for conciseness.
! SOURCE      Mod.atm.    Teffs         logg        [M/H]     Vmicro
! VanHam93a:  ATLAS9   3500 - 50000   0.0 - 5.0      0.0      n/a
! DiazCo95a:  ATLAS9   3500 - 50000   0.0 - 5.0      0.0      n/a
! Claret00a:  ATLAS9   3500 - 50000   0.0 - 5.0  -5.0 - +0.5  0,1,2,4,8
! Claret00p:  Phoenix  2000 - 9800    3.5 - 5.0      0.0      2
! ClarHa03a:  Phoenix  5000 - 10000   3.5 - 5.5      0.0      2
! Claret04a:  ATLAS9   3500 - 50000   0.0 - 5.0  -5.0 - +1.0  0,1,2,4,8
! Claret04p:  Phoenix  2000 - 9800    3.5 - 5.0      0.0      2
! Sing2010a:  ATLAS9   3500 - 40000   0.0 - 5.0  -5.0 - +1.0  n/a
!--------------------------------------------------------------------------------
! Available passbands and limb darkening laws:
! VanHam93a:  lin      log sqrt                   uvby UBV RcIc RjIj JKLMN
! DiazCo95a:  lin quad     sqrt                   uvby UBV
! Claret95a:  lin quad     sqrt                            RcIc JHK
! Claret00a:  lin quad log sqrt     4par          uvby UBV RcIc JHK
! Claret00p:  lin quad log sqrt     4par          uvby UBV RcIc JHK
! ClarHa03a:  lin quad log sqrt exp 4par          uvby UBV RcIc JHK
! Claret04a:  lin quad log sqrt     4par          ugriz
! Claret04p:  lin quad log sqrt     4par          ugriz
! Sing2010a:  lin quad              4par 3par     CoRoT Kepler
!================================================================================
! The LD tabulations must be regularly gridded in Teff, logg, Vmicro and [M/H].
! The actual values of these quantities are not important as long as those values
! are used consistently throughout the tabulations. The format for each line of a
! LD tabulations file is Teff Logg [M/H] Vmicro Filter Coefficients. Example:
!  3500 4.50  0.00  2  Kp  0.6842  0.4389  0.3036  1.5330 -1.0091  0.2420  0.4956
!  3750 4.50  0.00  2  Kp  0.6650  0.4017  0.3258  1.7807 -1.5071  0.4890  0.4790
!  4000 4.50  0.00  2  Kp  0.6888  0.5079  0.2239  1.6669 -1.4738  0.5729  0.5478
!  4250 4.50  0.00  2  Kp  0.7215  0.6408  0.0999  1.2877 -0.9263  0.4009  0.6928
!  4500 4.50  0.00  2  Kp  0.7163  0.6483  0.0842  1.1412 -0.6724  0.2793  0.7393
!================================================================================
! Warning numbers:  -1000 means the requested Teff or logg is out of range
!                  9.9999 no coeffs exist for that law, filter, [M/H] or Vmicro
!                   -2000 the specified LD table doesn't have the specified law
!================================================================================
! Installation:
! It is a monolithic program in the FORTRAN77 language, so is easy to compile and
! run. It was originally compiled with Gnu G77 but now I use the ifort compiler.
! It also compiles and runs successfully using gfortran and g95.
! Apart from the program you will also need to the LD tabulations. Put them into
! a directory somewhere on your computer. Then put the full file names and paths
! into the "data FILES" array (line ~128 below) and compile the program.
!================================================================================
! History:
!                 Started as a subroutine in JKTEBOP doing linear LD coeffs only
!                 Modified for nonlinear LD coefficients as a subroutine in JKTWD
! V1: 2006 Oct:   Became a stand-alone program with many more options
! V2: 2007/01/31: Much tidying up and fixed a bug with the VanHam93a LD coeffs
! V3: 2010/12/05: Added in the Sing2010a LD coefficients and tidied up further
!================================================================================
      implicit none
      real*8 TEFF,LOGG               ! Teff and log(g) to interpolate to
      real*8 MOH,VMICRO              ! [M/H] and Vmicro values to select
      character*1 LAW*1              ! Type of LD law (u, q, l or s)
      integer TABLE                  ! LD coefficient table to use
      character*2 FILTER             ! Single filter name or designation
      character*30 OUTFILE           ! Output file (if wanted)
      integer NARG                   ! Number of command-line arguments
      integer i,j,k,ERROR            ! Loop counters and error flag
      real*8 INITVAL                 ! Initialize arrays to this value

      integer NFILTER                ! Number of filter designations to do
      parameter ( NFILTER=25 )       ! bo uvby UBVRIJHKLMN RcIc ugriz C K
      character*2 FILTERS(NFILTER)   ! 2-character filter designations
      character*11 FILTERN(NFILTER)  ! Longer filter names
      real*8 VANHAM93(NFILTER,8)     ! VanHam93a: u Q(u) l1 l2 Q(l) s1 s2 Q(s)
      real*8 DIAZCO95(NFILTER,5)     ! DiazCo95a: u q1 q2 s1 s2
      real*8 CLARET00A(NFILTER,11)   ! Claret00a: u q1 q2 s1 s2 l1 l2 f1 f2 f3 f4
      real*8 CLARET00P(NFILTER,11)   ! Claret00p: u q1 q2 s1 s2 l1 l2 f1 f2 f3 f4
      real*8 CLARHA03(NFILTER,13)    ! ClarHa03p: u q1q2 s1s2 l1l2 e1e2 f1f2f3f4
      real*8 CLARET04A(NFILTER,11)   ! Claret04a: u q1 q2 s1 s2 l1 l2 f1 f2 f3 f4
      real*8 CLARET04P(NFILTER,11)   ! Claret04a: u q1 q2 s1 s2 l1 l2 f1 f2 f3 f4
      real*8 SING2010(NFILTER,10)    ! Sing2010a: u q1 q2 t1 t2 t3 f1 f2 f3 f4
      real*8 OUTARRAY(NFILTER,13)    ! Used when only outputting a few LDC

      integer NFILES                 ! Number of LD coefficient tabulations
      parameter ( NFILES=8 )         ! Number of LD coefficient tabulations
      character*50 FILES(NFILES)     ! Contains names of the files to do
      integer NCOL(NFILES)           ! Number of columns of LD coeff in each file

      data NCOL / 8,5,11,11,13,11,11,10 /

      data FILTERS / "bo","uS","vS","bS","yS","UJ","BJ","VJ",
     &               "RC","RJ","IC","IJ","JJ","HJ","KJ","LJ","MJ","NJ",
     &               "Su","Sg","Sr","Si","Sz","Co","Kp"/

      data FILTERN /     "bolometric ","Stromgren u","Stromgren v",
     &     "Stromgren b","Stromgren y","Johnson U  ","Johnson B  ",
     &     "Johnson V  ","Cousins R  ","Johnson R  ","Cousins I  ",
     &     "Johnson I  ","Johnson J  ","Johnson H  ","Johnson K  ",
     &     "Johnson L  ","Johnson M  ","Johnson N  ","Sloan u    ",
     &     "Sloan g    ","Sloan r    ","Sloan i    ","Sloan z    ",
     &     "CoRoT white","Kepler     "  /

      data FILES / "JKTLD-vanhamme.dat"    ,
     &             "JKTLD-diazcordoves.dat",
     &             "JKTLD-claret2000a.dat" ,
     &             "JKTLD-claret2000p.dat" ,
     &             "JKTLD-claret2003.dat"  ,
     &             "JKTLD-claret2004a.dat" ,
     &             "JKTLD-claret2004p.dat" ,
     &             "JKTLD-sing2010.dat"    /

      INITVAL = 99.9999d0
      do i = 1,NFILTER
        do j = 1,13
          OUTARRAY(i,j)  = INITVAL
          CLARHA03(i,j)  = INITVAL
          if ( j <= 8 ) VANHAM93(i,j) = INITVAL
          if ( j <= 5 ) DIAZCO95(i,j) = INITVAL
          if ( j <= 10 ) SING2010(i,j) = INITVAL
          if ( j <= 11 ) then
            CLARET00A(i,j) = INITVAL
            CLARET00P(i,j) = INITVAL
            CLARET04A(i,j) = INITVAL
            CLARET04P(i,j) = INITVAL
          end if
        end do
      end do

            ! Now read in all of the command-line arguments. If the
            ! wrong number were given then output helpful info instead.

      CALL READIN (NARG,TEFF,LOGG,MOH,VMICRO,LAW,TABLE,FILTER,OUTFILE)

            ! If only a few coefficients are wanted to standard output
            ! then calculate all coefficients for the LD file requested,
            ! find the index referring to the correct filter designation
            ! and output the results dependent on column number.

      if ( NARG == 7 ) then
        CALL DOFILE (FILES(TABLE),NCOL(TABLE),NFILTER,FILTERS,
     &                                    TEFF,LOGG,MOH,VMICRO,OUTARRAY)
        do i = 1,NFILTER
          if ( FILTER == FILTERS(i) ) exit
        end do
        if ( TABLE == 1 ) then
          if ( LAW == "u" ) then
            write (*,100) OUTARRAY(i,1)
          else if ( LAW == "l" ) then
            write (*,100) OUTARRAY(i,3),OUTARRAY(i,4)
          else if ( LAW == "s" ) then
            write (*,100) OUTARRAY(i,6),OUTARRAY(i,7)
          else
            write (*,*) " -2000   -2000"
          end if
        else if ( TABLE == 2 ) then
          if ( LAW == "u" ) then
            write (*,100) OUTARRAY(i,1)
          else if ( LAW == "q" ) then
            write (*,100) OUTARRAY(i,2),OUTARRAY(i,3)
          else if ( LAW == "s" ) then
            write (*,100) OUTARRAY(i,4),OUTARRAY(i,5)
          else
            write (*,*) " -2000   -2000"
          end if
        else if ( TABLE == 8 ) then
          if ( LAW == "u" ) then
            write (*,100) OUTARRAY(i,1)
          else if ( LAW == "q" ) then
            write (*,100) OUTARRAY(i,2),OUTARRAY(i,3)
          else if ( LAW == "t" ) then
            write (*,100) OUTARRAY(i,4),OUTARRAY(i,5),OUTARRAY(i,6)
          else if ( LAW == "f" ) then
            write (*,100) (OUTARRAY(i,j),j=7,10)
          else
            write (*,*) " -2000   -2000"
          end if
        else
          if ( LAW == "u" ) write (*,100) OUTARRAY(i,1)
          if ( LAW == "q" ) write (*,100) OUTARRAY(i,2),OUTARRAY(i,3)
          if ( LAW == "l" ) write (*,100) OUTARRAY(i,6),OUTARRAY(i,7)
          if ( LAW == "s" ) write (*,100) OUTARRAY(i,4),OUTARRAY(i,5)
          if ( LAW == "e" ) write (*,100) OUTARRAY(i,8),OUTARRAY(i,9)
          if ( LAW == "f" ) then
            if ( TABLE == 5 ) then
              write (*,100) (OUTARRAY(i,j),j=10,13)
            else
              write (*,100) (OUTARRAY(i,j),j=8,11)
            end if
          end if
        end if
        STOP
      end if
100   FORMAT (4(F10.4,1X))

            ! If a full set of coefficients are wanted then calculate
            ! the results for every LD file and output them to file.

      write (*,'(A24,A50)') "Calculating results for ",FILES(1)
      CALL DOFILE (FILES(1),NCOL(1),NFILTER,FILTERS,TEFF,LOGG,
     &                                              MOH,VMICRO,VANHAM93)
      write (*,'(A24,A50)') "Calculating results for ",FILES(2)
      CALL DOFILE (FILES(2),NCOL(2),NFILTER,FILTERS,TEFF,LOGG,
     &                                              MOH,VMICRO,DIAZCO95)
      write (*,'(A24,A50)') "Calculating results for ",FILES(3)
      CALL DOFILE (FILES(3),NCOL(3),NFILTER,FILTERS,TEFF,LOGG,
     &                                             MOH,VMICRO,CLARET00A)
      write (*,'(A24,A50)') "Calculating results for ",FILES(4)
      CALL DOFILE (FILES(4),NCOL(4),NFILTER,FILTERS,TEFF,LOGG,
     &                                             MOH,VMICRO,CLARET00P)
      write (*,'(A24,A50)') "Calculating results for ",FILES(5)
      CALL DOFILE (FILES(5),NCOL(5),NFILTER,FILTERS,TEFF,LOGG,
     &                                              MOH,VMICRO,CLARHA03)
      write (*,'(A24,A50)') "Calculating results for ",FILES(6)
      CALL DOFILE (FILES(6),NCOL(6),NFILTER,FILTERS,TEFF,LOGG,
     &                                             MOH,VMICRO,CLARET04A)
      write (*,'(A24,A50)') "Calculating results for ",FILES(7)
      CALL DOFILE (FILES(7),NCOL(7),NFILTER,FILTERS,TEFF,LOGG,
     &                                             MOH,VMICRO,CLARET04P)
      write (*,'(A24,A50)') "Calculating results for ",FILES(8)
      CALL DOFILE (FILES(8),NCOL(8),NFILTER,FILTERS,TEFF,LOGG,
     &                                              MOH,VMICRO,SING2010)

      CALL OUTTOFILE (OUTFILE,TEFF,LOGG,MOH,VMICRO,NFILTER,FILTERS,
     &                FILTERN,VANHAM93,DIAZCO95,CLARET00A,CLARET00P,
     &                CLARHA03,CLARET04A,CLARET04P,SING2010)

      END PROGRAM JKTLD
!=======================================================================
      SUBROUTINE READIN (NARG,TEFF,LOGG,MOH,VMICRO,LAW,TABLE,FILTER,
     &                                                          OUTFILE)
            ! If ther are command-line arguments it reads them in and
            ! checks them. If not, then it outputs useful information.
      implicit none
      integer NARG                  ! OUT: Number command-line arguments
      real*8 TEFF,LOGG              ! OUT: Teff,log(g) to interpolate to
      real*8 MOH,VMICRO             ! OUT: [M/H],Vmicro values to select
      character*1 LAW*1             ! OUT: Type of LD law (u, q, l or s)
      integer TABLE                 ! OUT: LD coefficient table  to  use
      character*2 FILTER            ! OUT: Single filter designation
      character*30 OUTFILE          ! OUT: Output file (if wanted)
      character*10 INPUT            ! LOCAL: Help read command-line args
      integer i,j,k,m,ERROR         ! LOCAL: Loop counters and errorflag

      NARG = iargc()

      if ( NARG /= 5 .and. NARG /= 7 ) then
        write (*,*) " "
        write (*,'(A38,A42)') "JKTLD                John Southworth  ",
     &                    "(Keele University)   jkt~astro.keele.ac.uk"
        write (*,'(A38,A42)') "This code outputs interpolated limb da",
     &                    "rkening coefficients for a given Teff,logg"
        write (*,'(A38,A42)') "for several LD laws, and using tabulat",
     &                    "ed coefficients from several publications."
        write (*,'(A38,A42)') "--------------------------------------",
     &                    "------------------------------------------"
        write (*,'(A38,A42)') "To output a large number of coefficien",
     &                    "ts for a given Teff, logg, [M/H], Vmicro: "
        write (*,'(A38,A42)') "  Usage:  jktld  <Teff>  <logg>  <M/H>",
     &                    "  <Vmicro>  <outfile>                     "
        write (*,'(A38,A42)') "--------------------------------------",
     &                    "------------------------------------------"
        write (*,'(A38,A42)') "To output coefficients for one Teff, l",
     &                    "ogg, [M/H], Vmicro, law, passband, source:"
        write (*,'(A38,A42)') "  Usage:  jktld  <Teff>  <logg>  <M/H>",
     &                    "  <Vmicro>  <law>  <table>  <filter>      "
        write (*,'(A38,A42)') "<law>  =  'u' for linear, 'q' for quad",
     &                    "ratic, 'l' for logarithmic, 's' for square"
        write (*,'(A38,A42)') "          root,  'e' for exponential, ",
     &                    " 't' for Sing 3-par,  'f' for Claret 4-par"
        write (*,'(A38,A42)') "<table>  =  1 VanHam93, 2 DiazCo95+Cla",
     &                    "ret95, 3 Claret00 ATLAS, 4 Claret00 Phoenx"
        write (*,'(A38,A42)') "            5 ClarHa03, 6 Claret04 ATL",
     &                    "AS, 7 for Claret04 Phoenix, 8 for Sing2010"
        write (*,'(A38,A42)') "<f> = bo uS vS bS yS UJ BJ VJ RC RJ IC",
     &                    " IJ JJ HJ HK LJ MJ NJ Su Sg Sr Si Sz Co Kp"
        write (*,*) " "
        STOP

      else
        CALL GETARG (1,INPUT)
        read (INPUT,*,iostat=ERROR) TEFF
        if ( ERROR /= 0 ) then
          write (*,*) " "
          write(*,'(A45)')"### ERROR reading Teff from the command line"
          write (*,*) " "
          STOP
        end if

        CALL GETARG (2,INPUT)
        read (INPUT,*,iostat=ERROR) LOGG
        if ( ERROR /= 0 ) then
          write (*,*) " "
          write(*,'(A45)')"### ERROR reading logg from the command line"
          write (*,*) " "
          STOP
        end if

        CALL GETARG (3,INPUT)
        read (INPUT,*,iostat=ERROR) MOH
        if ( ERROR /= 0 ) then
          write (*,*) " "
          write(*,'(A45)')"### ERROR reading [M/H] from command line   "
          write (*,*) " "
          STOP
        end if

        CALL GETARG (4,INPUT)
        read (INPUT,*,iostat=ERROR) VMICRO
        if ( ERROR /= 0 ) then
          write (*,*) " "
          write(*,'(A45)')"### ERROR reading Vmicro from  command line "
          write (*,*) " "
          STOP
        end if


        if ( NARG == 5 ) then
          CALL GETARG (5,OUTFILE)

        else if ( NARG == 7 ) then

          CALL GETARG (5,INPUT)
          read (INPUT,*,iostat=ERROR) LAW
          if ( ERROR /= 0 ) then
            write (*,*) " "
            write (*,'(A35,A45)') "### ERROR reading the type of LD la",
     &                   "w to use from the command line.             "
            write (*,*) " "
            STOP
          else if ( LAW /= "u" .and. LAW /= "q" .and. LAW /= "l"
     &              .and. LAW /= "s" .and. LAW /= "e" .and. LAW /= "f"
     &              .and. LAW /= "t" ) then
            write (*,*) " "
            write (*,'(A35,A45)') "### ERROR: type of LD law must be o",
     &                  "ne of the following: u q l s e t f           "
            write (*,*) " "
            STOP
          end if

          CALL GETARG (6,INPUT)
          read (INPUT,*,iostat=ERROR) TABLE
          if ( ERROR /= 0 ) then
            write (*,*) " "
            write (*,'(A35,A45)') "### ERROR reading the number of the",
     &                   " LD coefficient table from the command line."
            write (*,*) " "
            STOP
          else if ( TABLE < 1 .or. TABLE > 8 ) then
            write (*,*) " "
            write (*,'(A35,A45)') "### ERROR: the number of the LD coe",
     &                   "fficient table must be between 1 and 8      "
            write (*,*) " "
            STOP
          end if

          CALL GETARG (7,INPUT)
          read (INPUT,*,iostat=ERROR) FILTER
          if ( ERROR /= 0 ) then
            write (*,*) " "
            write (*,'(A35,A45)') "### ERROR reading the passband from",
     &                   " the command line                           "
            write (*,*) " "
            STOP
          end if

        end if

      end if

      END SUBROUTINE READIN
!=======================================================================
      SUBROUTINE OUTTOFILE (OUTFILE,TEFF,LOGG,MOH,VMICRO,NFILTER,
     &                      FILTERS,FILTERN,VANHAM93,DIAZCO95,CLARET00A,
     &                  CLARET00P,CLARHA03,CLARET04A,CLARET04P,SING2010)
      implicit none
      character*30 OUTFILE          ! IN: Name of file to write to
      real*8 TEFF,LOGG              ! IN: Teff and log(g) used
      real*8 MOH,VMICRO             ! IN: [M/H] and Vmicro values used
      integer i,j,k,m,ERROR         ! LOCAL: Loop counters and errorflag

      integer NFILTER               ! IN: Number of filters
      character*2 FILTERS(NFILTER)  ! IN: Short names of the filters
      character*11 FILTERN(NFILTER) ! IN: Longer names of the filters
      real*8 VANHAM93(NFILTER,8)    ! IN: VanHam93: u uQ l1 l2 cdQ s1 s2 efQ
      real*8 DIAZCO95(NFILTER,5)    ! IN: DiazCo95:        u q1 q2 s1 s2
      real*8 CLARET00A(NFILTER,7)   ! IN: Claret00(ATL):   u q1 q2 s1 s2 l1 l2
      real*8 CLARET00P(NFILTER,7)   ! IN: Claret00(Phe):   u q1 q2 s1 s2 l1 l2
      real*8 CLARHA03(NFILTER,9)    ! IN: ClarHa03(Phe):   u q1 q2 s1 s2 l1 l2 gh
      real*8 CLARET04A(NFILTER,7)   ! IN: Claret04(ATL):   u q1 q2 s1 s2 l1 l2
      real*8 CLARET04P(NFILTER,7)   ! IN: Claret04(Phe):   u q1 q2 s1 s2 l1 l2
      real*8 SING2010(NFILTER,10)   ! IN: Sing2010a:     u q1 q2 3par 4par

      open (50,file=OUTFILE,status="new",iostat=ERROR)
      if ( ERROR /= 0 ) then
        write (*,'(A50,A30)')
     &      "### ERROR: could not open the output file, name:  ",OUTFILE
        write (*,*) " "
        STOP
      end if

      write (50,'(A40,A40)') "========================================",
     &                       "========================================"
      write (50,'(A40,A40)') "Output from JKTLD                      J",
     &                       "ohn Southworth   (jkt~astro.keele.ac.uk)"
      write (50,'(A40,A40)') "Limb darkening coefficients bilinearly i",
     &                       "nterpolated from published LD tables for"
      write (50,'(A7,F7.1,A13,F6.4,A11,A10,F4.1,A13,F3.1,A5)')
     &                "Teff = ",TEFF," K    logg = ",LOGG," [cm/s^2]  ",
     &                  "  [M/H] = ",MOH,"    Vmicro = ",VMICRO," km/s"
      write (50,'(A40,A40)') "----------------------------------------",
     &                       "----------------------------------------"
      write (50,'(A40,A40)') "'9.9999' indicates unavailable result   ",
     &                       "  '*******' means parameter out of range"
      write (50,'(A40,A40)') "Note: the exponential, three-par and fou",
     &                       "r-par laws are not included below. Their"
      write (50,'(A40,A40)') "coefficients can be obtained using the o",
     &                       "ption to output results to the terminal."
      write (50,'(A40,A40)') "========================================",
     &                       "========================================"
      write (50,'(A40,A40)') "Linear limb darkening law:  I(mu)/I(0) =",
     &                       " 1 - u*(1 - mu)          mu = cos(theta)"
      write (50,'(A40,A40)') "Quadratic LD law:           I(mu)/I(0) =",
     &                       " 1 - q1*(1 - mu) - q2*(1 - mu)^2        "
      write (50,'(A40,A40)') "Logarithmic LD law:         I(mu)/I(0) =",
     &                       " 1 - l1*(1 - mu) - l2*mu*ln(mu)         "
      write (50,'(A40,A40)') "Square-root LD law:         I(mu)/I(0) =",
     &                       " 1 - s1*(1 - mu) - s2*(1 - sqrt[mu])    "
      write (50,'(A40,A40)') "Exponential LD law:         I(mu)/I(0) =",
     &                       " 1 - e1*(1 - mu) - e2/(1 - exp[mu])     "
      write (50,'(A40,A40)') "Sing three-parameter law:   I(mu)/I(0) =",
     &                       " 1 - sum_(i=2..4) ti*[1 - mu^(i/2)]     "
      write (50,'(A40,A40)') "Claret four-parameter law:  I(mu)/I(0) =",
     &                       " 1 - sum_(i=1..4) f_i*[1 - mu^(i/2)]    "
      write (50,'(A40,A40)') "----------------------------------------",
     &                       "----------------------------------------"
      write (50,'(A40,A40)') "Source references for the LD coefficient",
     &                       "s:                                      "
      write (50,'(A40,A40)') "VanHam93a:  Van Hamme (1993AJ....106.209",
     &                       "6V)                                     "
      write (50,'(A40,A40)') "DiazCo95a:  Diaz-Cordoves+ (1995A+AS..11",
     &                       "0..329D), Claret+ (1995A+AS..114..247C) "
      write (50,'(A40,A40)') "Claret00a:  Claret (2000A+A...363.1081C)",
     &                       "                                        "
      write (50,'(A40,A40)') "Claret00p:  Claret (2000A+A...363.1081C)",
     &                       "                                        "
      write (50,'(A40,A40)') "ClarHa03a:  Claret & Hauschildt (2003A+A",
     &                       "...412..241C)                           "
      write (50,'(A40,A40)') "Claret04a:  Claret (2004A+A...428.1001C)",
     &                       "                                        "
      write (50,'(A40,A40)') "Claret04p:  Claret (2004A+A...428.1001C)",
     &                       "                                        "
      write (50,'(A40,A40)') "Sing2010a:  Sing (2010A+A...510A..21S)  ",
     &                       "                                        "
      write (50,'(A40,A40)') "----------------------------------------",
     &                       "----------------------------------------"
      write (50,'(A40,A40)') "Model atmospheres used by each source:  ",
     &                       "                                        "
      write (50,'(A40,A40)') "Source      Mod.atm.    Teffs         lo",
     &                       "gg        [M/H]     Vmicro              "
      write (50,'(A40,A40)') "VanHam93a:  ATLAS9   3500 - 50000   0.0 ",
     &                       "- 5.0      0.0      n/a                 "
      write (50,'(A40,A40)') "DiazCo95a:  ATLAS9   3500 - 50000   0.0 ",
     &                       "- 5.0      0.0      n/a                 "
      write (50,'(A40,A40)') "Claret00a:  ATLAS9   3500 - 50000   0.0 ",
     &                       "- 5.0  -5.0 - +0.5  0,1,2,4,8           "
      write (50,'(A40,A40)') "Claret00p:  Phoenix  2000 - 9800    3.5 ",
     &                       "- 5.0      0.0      2                   "
      write (50,'(A40,A40)') "ClarHa03a:  Phoenix  5000 - 10000   3.5 ",
     &                       "- 5.5      0.0      2                   "
      write (50,'(A40,A40)') "Claret04a:  ATLAS9   3500 - 50000   0.0 ",
     &                       "- 5.0  -5.0 - +1.0  0,1,2,4,8           "
      write (50,'(A40,A40)') "Claret04p:  Phoenix  2000 - 9800    3.5 ",
     &                       "- 5.0      0.0      2                   "
      write (50,'(A40,A40)') "Sing2010a:  ATLAS9   3500 - 40000   0.0 ",
     &                       "- 5.0  -5.0 - +1.0  n/a                 "
      write (50,'(A40,A40)') "----------------------------------------",
     &                       "----------------------------------------"
      write (50,'(A40,A40)') "Available passbands and limb darkening l",
     &                       "aws:                                    "
      write (50,'(A40,A40)') "VanHam93a:  lin      log sqrt           ",
     &                       "        uvby UBV RcIc RjIj JKLMN        "
      write (50,'(A40,A40)') "DiazCo95a:  lin quad     sqrt           ",
     &                       "        uvby UBV                        "
      write (50,'(A40,A40)') "Claret95a:  lin quad     sqrt           ",
     &                       "                 RcIc JHK               "
      write (50,'(A40,A40)') "Claret00a:  lin quad log sqrt     4par  ",
     &                       "        uvby UBV RcIc JHK               "
      write (50,'(A40,A40)') "Claret00p:  lin quad log sqrt     4par  ",
     &                       "        uvby UBV RcIc JHK               "
      write (50,'(A40,A40)') "ClarHa03a:  lin quad log sqrt exp 4par  ",
     &                       "        uvby UBV RcIc JHK               "
      write (50,'(A40,A40)') "Claret04a:  lin quad log sqrt     4par  ",
     &                       "        ugriz                           "
      write (50,'(A40,A40)') "Claret04p:  lin quad log sqrt     4par  ",
     &                       "        ugriz                           "
      write (50,'(A40,A40)') "Sing2010a:  lin quad              4par 3",
     &                       "par     CoRoT Kepler                    "
      write (50,'(A40,A40)') "========================================",
     &                       "========================================"
      write (50,*) " "

      if ( TEFF > 5800.0d0 .and. TEFF < 7000.0d0 ) then
        write (50,'(A38,A42)') "WARNING: there are no coefficients bet",
     &                     "ween 5800 and 7000K in the Claret (2000,  "
        write (50,'(A38,A42)') "2004) Phoenix tables. The gap is inter",
     &                     "polated over but results may be inaccurate"
        write (50,*) " "
      end if


      write (50,'(A40,A40)') "General overview for linear, log, sqrt l",
     &                       "aws with quality indicators (VanHam93a):"
      write (50,'(A40,A40)') "Passband          u    quality     l1   ",
     &                       "  l2    quality     s1     s2    quality"
      do i = 1,18
        if ( FILTERS(i) /= "HJ" )
     &      write (50,'(A2,1X,A11,1X,F7.4,1X,F7.4,2(1X,3(1X,F7.4)))')
     &                       FILTERS(i),FILTERN(i),(VANHAM93(i,j),j=1,8)
      end do
      write (50,*) " "

      write (50,'(A40,A40)') "Limb darkening coefficients for the SDSS",
     &                       " passbands  (Claret 2004, ATLAS9):      "
      write (50,'(A40,A40)') "Passband              u         q1      ",
     &                       "q2        l1      l1        s1      s2  "
      do i = 19,23
            write (50,'(A2,2X,A11,4X,F7.4,3(3X,F7.4,1X,F7.4))')
     &                     FILTERS(i),FILTERN(i),(CLARET04A(i,j),j=1,3),
     &                    (CLARET04A(i,j),j=6,7),(CLARET04A(i,j),j=4,5)
      end do

      write (50,'(A40,A40)') "Limb darkening coefficients for the SDSS",
     &                       " passbands (Claret 2004, Phoenix):      "
      write (50,'(A40,A40)') "Passband              u         q1      ",
     &                       "q2        l1      l2        s1      s2  "
      do i = 19,23
            write (50,'(A2,2X,A11,4X,F7.4,3(3X,F7.4,1X,F7.4))')
     &                     FILTERS(i),FILTERN(i),(CLARET04P(i,j),j=1,3),
     &                    (CLARET04P(i,j),j=6,7),(CLARET04P(i,j),j=4,5)
      end do
      write (50,*) " "

      write (50,'(A40,A40)') "Limb darkening coefficients for the CoRo",
     &                       "T and Kepler passbands (from Sing 2010):"
      write (50,'(A40,A40)') "(The CoRoT colour passbands are differen",
     &                       "t for every star so cannot be predicted)"
      write (50,'(A40,A40)') "Passband              u         q1      ",
     &                       "q2        t1      t2      t3            "
      do i = 24,25
            write (50,'(A2,2X,A11,4X,F7.4,2(3X,F7.4,1X,F7.4),1X,F7.4)')
     &                     FILTERS(i), FILTERN(i), SING2010(i,1),
     &                    (SING2010(i,j),j=2,3), (SING2010(i,j),j=4,6)
      end do
      write (50,*) " "

      write (50,'(A44)') "Linear limb darkening coefficients          "
      write (50,'(A44)') "----------------------------------          "
      write (50,'(A40,A40)') "       VanHam93    DiazCo95    Claret00a",
     &                       "   Claret00p   ClarHa03                 "
      do i = 1,18
        write (50,'(A2,7(5X,F7.4))') FILTERS(i),VANHAM93(i,1),
     &         DIAZCO95(i,1),CLARET00A(i,1),CLARET00P(i,1),CLARHA03(i,1)
      end do
      write (50,*) " "

      write (50,'(A44)') "Quadratic limb darkening coefficients       "
      write (50,'(A44)') "-------------------------------------       "
      write (50,'(A40,A40)') "         DiazCl93          Claret00a    ",
     &                       "     Claret00p         ClarHa03         "
      write (50,'(A40,A40)') "        q1      q2        q1      q2    ",
     &                       "    q1      q2        q1      q2        "
      do i = 1,18
        if ( FILTERS(i) /= "LJ" .and. FILTERS(i) /= "MJ" .and.
     &       FILTERS(i) /= "NJ" .and. FILTERS(i) /= "RJ" .and.
     &       FILTERS(i) /= "IJ")    write (50,'(A2,8(3X,F7.4,1X,F7.4))')
     &          FILTERS(i),(DIAZCO95(i,j),j=2,3),(CLARET00A(i,j),j=2,3),
     &                      (CLARET00P(i,j),j=2,3),(CLARHA03(i,j),j=2,3)
      end do
      write (50,*) " "

      write (50,'(A44)') "Logarithmic limb darkening coefficients     "
      write (50,'(A44)') "---------------------------------------     "
      write (50,'(A40,A40)') "         VanHam93          Claret00a    ",
     &                       "     Claret00p         ClarHa03         "
      write (50,'(A40,A40)') "        l1      l2        l1      l2    ",
     &                       "    l1      l2        l1      l2        "
      do i = 1,18
        write (50,'(A2,8(3X,F7.4,1X,F7.4))')                 FILTERS(i),
     &                     (VANHAM93(i,j),j=3,4),(CLARET00A(i,j),j=6,7),
     &                     (CLARET00P(i,j),j=6,7),(CLARHA03(i,j),j=6,7)
      end do
      write (50,*) " "

      write (50,'(A44)') "Square-root limb darkening coefficients     "
      write (50,'(A44)') "---------------------------------------     "
      write (50,'(A40,A52)') "         VanHam93          DiazCo95     ",
     &           "     Claret00a         Claret00p         ClarHa03   "
      write (50,'(A40,A52)') "        s1      s2        s1      s2    ",
     &           "    s1      s2        s1      s2        s1      s2  "
      do i = 1,18
        write (50,'(A2,8(3X,F7.4,1X,F7.4))')                 FILTERS(i),
     &           (VANHAM93(i,j),j=6,7),(DIAZCO95(i,j),j=4,5),(CLARET00A
     &         (i,j),j=4,5),(CLARET00P(i,j),j=4,5),(CLARHA03(i,j),j=4,5)
      end do
      write (50,*) " "

!       write (50,'(A44)') "Exponential limb darkening coefficients     "
!       write (50,'(A44)') "---------------------------------------     "
!       write (50,'(A44)') "         ClarHa03                           "
!       write (50,'(A44)') "        e1      e2                          "
!       do i = 1,18
!         if ( FILTERS(i) /= "bo" .and. FILTERS(i) /= "RJ" .and.
!      &       FILTERS(i) /= "IJ" .and. FILTERS(i) /= "LJ" .and.
!      &       FILTERS(i) /= "MJ" .and. FILTERS(i) /= "NJ" )        write
!      &    (50,'(A2,8(3X,F7.4,1X,F7.4))')FILTERS(i),(CLARHA03(i,j),j=8,9)
!       end do
!       write (50,*) " "

      write (50,'(A40,A40)') "========================================",
     &                       "========================================"
      write (50,*) " "
      close (50)

      END SUBROUTINE OUTTOFILE
!=======================================================================
!=======================================================================
      SUBROUTINE DOFILE (FILE,NCOL,NFILTER,FILTERS,TEFF,LOGG,MOH,VMICRO,
     &                                                           OUTPUT)
      implicit none
      character*50 FILE             ! IN: File containing LD coeffs
      integer NCOL                  ! IN: Number of columns of LD coeffs
      integer NFILTER               ! IN: Number of filters to do
      character*2 FILTERS(NFILTER)  ! IN: Names of the filters
      real*8 TEFF,LOGG              ! IN: Teff, log(g) to interpolate to
      real*8 MOH,VMICRO             ! IN: [M/H], Vmicro values to select
      real*8 OUTPUT(NFILTER,NCOL)   ! OUT: interpolated LD coefficients
      integer NDATA                 ! LOCAL: Number of file lines read
      real*8 INDEP(4,200000)        ! LOCAL: Input Teff,logg,[M/H],Vmicr
      character*2 BAND(200000)      ! LOCAL: Incoming filter designat'ns
      real*8 DEP(13,200000)         ! LOCAL: Incoming LD coefficients
      real*8 TABARRAY(15,100000)    ! LOCAL: All coeffs for interpolat'n
      real*8 BILINEAR               ! FUNCTION do bilinear interpolation
      real*8 VSMALL                 ! LOCAL: small number
      integer i,j,k,m,ERROR         ! LOCAL: Loop counters and errorflag

      VSMALL = 1.0d-6

      open (51,file=FILE,status="old",iostat=ERROR)
      if ( ERROR /= 0 ) then
        write (*,'(A32,A50)') "### ERROR: problem opening file ",FILE
        write (*,*) " "
        STOP
      end if

      do j = 1,200000
        read (51,*,iostat=ERROR) (INDEP(k,j),k=1,4),BAND(j),
     &                           (DEP(k,j),k=1,NCOL)
        if ( ERROR /= 0 ) exit
      end do
      NDATA = j - 1
      close (51)

      do i = 1,NFILTER
        m = 0
        do j = 1,NDATA
          if ( BAND(j) == FILTERS(i) .and. abs(MOH-INDEP(3,j)) < VSMALL
     &         .and. abs(VMICRO-INDEP(4,j)) < VSMALL )  then
            m = m + 1
            TABARRAY(1,m) = INDEP(1,j)
            TABARRAY(2,m) = INDEP(2,j)
            do k = 1,NCOL
              TABARRAY(k+2,m) = DEP(k,j)
            end do
          end if
        end do

        if ( m > 5 ) then
          do k = 1,NCOL
            OUTPUT(i,k) = BILINEAR(TABARRAY,m,TEFF,LOGG,k+2)
          end do
        end if
      end do

      END SUBROUTINE DOFILE
!=======================================================================
      DOUBLEPRECISION FUNCTION BILINEAR (ARRAY,NDATA,IN1,IN2,COL)
            ! This performs bilinear interpolation on ARRAY to find the
            ! interpolated value for IN1 and IN2. ARRAY is assumed to
            ! have the two discrete distributions in the first two
            ! columns and then 8 further columns, one of which
            ! (given by COL) contains the results for interpolation.
      implicit none
      integer NDATA                 ! IN:  length of array ARRAY
      real*8 ARRAY(15,100000)       ! IN:  array containing tabular data
      real*8 IN1,IN2                ! IN:  values to interpolate towards
      integer COL                   ! IN:  column containing quantities
      real*8 LOW1,HIGH1,LOW2,HIGH2  ! LOCAL: adjacent values to IN1, IN2
      real*8 COEFFS(4),WT(4)        ! LOCAL: coeffs to sum using weights
      integer i,j,ERROR             ! LOCAL: loop counter and error flag

      data COEFFS /-1.0d6,-1.0d6,-1.0d6,-1.0d6/
      LOW1 = -1.0d6
      HIGH1 = 1.0d6
      LOW2 = -1.0d6
      HIGH2 = 1.0d6

            ! First go through ARRAY to find the values which bracket
            ! the two values to interpolate to (IN1 and IN2). Avoid
            ! getting the same values for LOW and HIGH (which would
            ! occur if input value was exactly the same as a tabulated
            ! one) by careful choice of inequalities.

      do i = 1,NDATA
        if (ARRAY(1,i) <  IN1)  LOW1  = max(LOW1, ARRAY(1,i))
        if (ARRAY(1,i) >= IN1)  HIGH1 = min(HIGH1,ARRAY(1,i))
        if (ARRAY(2,i) <  IN2)  LOW2  = max(LOW2, ARRAY(2,i))
        if (ARRAY(2,i) >= IN2)  HIGH2 = min(HIGH2,ARRAY(2,i))
      end do

      if (LOW1<-1.0d5.or.HIGH1>1.0d5.or.LOW2<-1.0d5.or.HIGH2>1.0d5) then
        BILINEAR = -1000.0
        return
      end if

            ! Now find the coeffs which are valid for the four values
            ! which surround the wanted point in 2D space.

      do i = 1,NDATA
        if      ( ARRAY(1,i) == LOW1  .and. ARRAY(2,i) == LOW2 ) then
          COEFFS(1) = ARRAY(COL,i)
        else if ( ARRAY(1,i) == HIGH1  .and. ARRAY(2,i) == LOW2) then
          COEFFS(2) = ARRAY(COL,i)
        else if ( ARRAY(1,i) == LOW1 .and. ARRAY(2,i) == HIGH2 ) then
          COEFFS(3) = ARRAY(COL,i)
        else if ( ARRAY(1,i) == HIGH1 .and. ARRAY(2,i) == HIGH2) then
          COEFFS(4) = ARRAY(COL,i)
        end if
      end do

            ! Now calculate the weights for the four values, which
            ! depend on how close they are to the wanted point.
            ! Total weights should add up to unity, so check this.

      WT(1) = (HIGH1-IN1) / (HIGH1-LOW1) * (HIGH2-IN2) / (HIGH2-LOW2)
      WT(2) = (IN1-LOW1)  / (HIGH1-LOW1) * (HIGH2-IN2) / (HIGH2-LOW2)
      WT(3) = (HIGH1-IN1) / (HIGH1-LOW1) * (IN2-LOW2)  / (HIGH2-LOW2)
      WT(4) = (IN1-LOW1)  / (HIGH1-LOW1) * (IN2-LOW2)  / (HIGH2-LOW2)

      BILINEAR = 0.0d0
      do i = 1,4
        BILINEAR = BILINEAR + WT(i)*COEFFS(i)
      end do

!       print*,"what ",in1,in2,col
!       print*,"hilo ",low1,high1,low2,high2
!       print*,"co   ",coeffs
!       print*,"wt   ",wt
!       print*,"out  ",BILINEAR,ndata

      END FUNCTION BILINEAR
!=======================================================================
!=======================================================================
