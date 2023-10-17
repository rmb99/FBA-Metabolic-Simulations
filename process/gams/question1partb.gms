********************Part b ****************
*** Specifying root dirctory of the project ***
$INLINECOM /*  */

$set myroot iAF1260/iAF1260


options
        limrow = 10000
        limcol = 10000
        optCR = 1E-9
        optCA = 0.0
        iterlim = 100000
        decimals = 8
        reslim = 100000
*work = 50000000
        sysout = off
        solprint = on;
        

*****************************
** Defenitions **
*      i                        = set of all metabolite
*      j                        =set of all reactions
*      offaeroglucose(j)        = set of reactions that should be turned off under aerobic glucose conditions
*      offanaero(j)             = set of reactions that should be turned off under aerobic glucose conditions
*      exchange(j)              = set of reactions defined in the minimal M9 media for uptake

*****************************

Sets

i
$include "%myroot%_cmp.txt"

j
$include "%myroot%_rxnnames.txt"

offaeroglucose(j)
$include "%myroot%_offaeroglu.txt"

offanaero(j)
$include "%myroot%_offanaero.txt"

exchange(j)
$include "%myroot%_source_M9.txt"
;

Parameters

S(i,j)
$include "%myroot%_sij.txt"

rxntype(j)
$include "%myroot%_rxntype.txt"

LB(j) , UB(j)
;


Variables

v(j)
z
;


***************************** DUAL VARIABLE DEFINITIONS **************************************
*
*       y(i)            = Dual variable associated with the stoichiometric constraints
*       y_minus(j)      = Dual variable associated with the lower bounds of reactions
*       y_plus(j)       = Dual variable associated with the upper bounds of reactions
*       w               = Dual objective value
*
*****************************************************************************************
Variables
y(i)
y_minus(j)
y_plus(j)
w
;

***************************** DUAL EQUATION DEFINITIONS **************************************
*
*       DualObj         = Objective function for the dual problem
*       DualStoic(j)    = Constraints derived from the stoichiometric matrix
*
*****************************************************************************************
Equations
DualObj
DualStoic(j)
;
Parameters
biomass_lo, biomass_up;

LB('Ec_biomass_iAF1260_core_59p81M') = 0;
UB('Ec_biomass_iAF1260_core_59p81M') = 0;

biomass_lo = v('Ec_biomass_iAF1260_WT_59p81M').lo;
biomass_up = v('Ec_biomass_iAF1260_WT_59p81M').up;
* DEFINING THE DUAL CONSTRAINTS AND OBJECTIVE FUNCTION *************************************************************
DualObj..           w =e= sum(j, LB(j)*y_minus(j) + UB(j)*y_plus(j));

DualStoic(j)..      sum(i, S(i,j)*y(i)) =g= (v('Ec_biomass_iAF1260_WT_59p81M').lo + v('Ec_biomass_iAF1260_WT_59p81M').up)/2;

y_minus.lo(j) = 0;
y_plus.lo(j) = 0;

*****************************************************************************************

************************ DECLARING THE DUAL MODEL with constraints****************************
Model dualModel
/
DualObj
DualStoic
/;

dualModel.optfile = 1;

*******************SOLVING THE Dual Linear Programming (LP) problem***************************
Solve dualModel using lp minimizing w ;

*******************EXPORTING THE DUAL RESULTS IN A TXT FILE***********************************
file ff_dual /DualBiomass.txt/;
put ff_dual;
put "The min value of the dual problem is : " w.l:0:8//;

*** loop through each metabolite i one at a time and print the name (with .tl extension)
*** and the level value of the dual variable (with .l extension)
loop(i,
    put i.tl:0:30, "    ", y.l(i):0:8/;
);
putclose ff_dual;
*****************************************************************************************
