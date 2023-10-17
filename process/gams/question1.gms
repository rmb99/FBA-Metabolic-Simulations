** Question one **

** In this question we aim to find the dimentionality of subspace spanning metablic fluxes consistent with
** maximum biomass formation under aerobic glucose condition and compare it with anerobic condition for
** iAF1260 GSM model of E-coli.!

***********************************
** Output files directory **



**********************************

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
*      offanaero(j)             = set of reactions that should be turned off under anerobic glucose conditions
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

***********************

Variables

v(j)
z
;

Equations
obj
Stoic
;


Obj..           z =e= v('Ec_biomass_iAF1260_WT_59p81M');
Stoic(i)..      sum(j, S(i,j) * v(j) )  =e=  0 ;


************************************************************************************

***** DEFINING THE UPPER AND LOWER BOUNDS FOR ALL THE REACTIONS *************************
scalar vmax /1000/;

LB(j)$(rxntype(j) = 0) = 0;
UB(j)$(rxntype(j) = 0) = vmax;
LB(j)$(rxntype(j) = 2) = 0;
UB(j)$(rxntype(j) = 2) = 0;
LB(j)$(rxntype(j) = 1) = -vmax;
UB(j)$(rxntype(j) = 1) = vmax;
LB(j)$(rxntype(j) = 4) = 0;
UB(j)$(rxntype(j) = 4) = vmax;
LB(j)$(exchange(j)) = -vmax;


LB('Ec_biomass_iAF1260_core_59p81M') = 0;
UB('Ec_biomass_iAF1260_core_59p81M') = 0;

*Reset bounds and constraints for aerobic conditions
LB(j)$(offaeroglucose(j)) = 0;
UB(j)$(offaeroglucose(j)) = 0;

*** Set lower bound of glucose uptake to -10 under aerobic conditions
LB('EX_glc(e)') = -10;
*** Set lower bound of oxygen uptake to -20 under aerobic conditions
LB('EX_o2(e)') = -20;
LB('ATPM') = 8.39;
UB('ATPM') = 8.39;



*** SETTING LOWER AND UPPER BOUNDS FOR THE FLUXES (variable bounds)
v.lo(j) = LB(j);
v.up(j) = UB(j);

v.fx('LDH_D_f') = 0;
v.fx('LDH_D2') = 0;
v.fx('L-LACD2') = 0;
v.fx('L-LACD3') = 0;
v.fx('PFL') = 0;


Model maxbio
/
Obj
Stoic
/;

maxbio.optfile = 1;

* Solve the model under aerobic glucose conditions
Solve maxbio using lp maximizing z;

* Record the dimensionality of the solution under aerobic glucose conditions
 
Scalar ActiveFluxesCountAerobic;
ActiveFluxesCountAerobic = sum(j, ord(j) * (v.l(j) > 1e-6));  

* Display the dimensionality results under aerobic glucose conditions
Display "Dimensionality of Active Fluxes under Aerobic Glucose Conditions:", ActiveFluxesCountAerobic;

*** Set lower bound of glucose uptake to -10 under anaerobic conditions
LB('EX_glc(e)') = -10;
*** Set lower bound of oxygen uptake to 0 under anaerobic conditions
LB('EX_o2(e)') = 0;

* Solve the model under anaerobic conditions
Solve maxbio using lp maximizing z;

* Record the maximum biomass value
Scalar PrimalMaxBiomass;
PrimalMaxBiomass = z.l;

* Record the dimensionality of the solution under anaerobic conditions
Scalar ActiveFluxesCountAnaerobic;
ActiveFluxesCountAnaerobic = sum(j, ord(j) * (v.l(j) > 1e-6));  
* Display the dimensionality results under anaerobic conditions
Display "Dimensionality of Active Fluxes under Anaerobic Conditions:", ActiveFluxesCountAnaerobic;

