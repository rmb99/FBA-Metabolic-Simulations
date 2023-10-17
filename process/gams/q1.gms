$inlineCom # # 
Sets
    i / ... /  # Define your set of metabolites here #
    j / ... /;  # Define your set of reactions here 

Parameters
    S(i,j) / ... / # Define the stoichiometry matrix (metabolite-reaction coefficients)
    LB(j) / ... /  # Define lower bounds for reactions
    UB(j) / ... /  # Define upper bounds for reactions
    RegulatoryConstraints(j) / ... /;  # Define regulatory constraints

Variables
    Flux(j);  # Define flux variables for reactions

Equations
    BiomassObjective,  # Objective to maximize biomass
    RegulatoryConstraintsEquations(j);  # Constraints to turn off reactions

BiomassObjective..  # Maximize biomass objective function
    sum(j, S('biomass_reaction',j) * Flux(j)) =e= 1;  # Adjust 'biomass_reaction' to match the actual biomass reaction name

* Apply flux bounds
Flux.L(j) = LB(j);
Flux.U(j) = UB(j);

* Apply regulatory constraints (turn off reactions)
RegulatoryConstraintsEquations(j).. Flux(j) =l= RegulatoryConstraints(j);

* Define the model
Model BiomassMaximization /BiomassObjective, RegulatoryConstraintsEquations/;

* Solve under aerobic glucose conditions
LB(j)$(offaeroglucose(j)) = 0;  $Turn off reactions for aerobic glucose
UB(j)$(offaeroglucose(j)) = 0;
RegulatoryConstraints(j)$(offaeroglucose(j)) = 0;
Solve BiomassMaximization using MIP maximizing Flux('biomass_reaction');

* Record the dimensionality of the solution
Variable ActiveFluxesCount;
ActiveFluxesCount = sum(j, ord(j) * (Flux(j) > 1e-6));  # Count non-zero fluxes

* Display the dimensionality results
Display "Dimensionality of Active Fluxes under Aerobic Glucose Conditions:", ActiveFluxesCount;

* Reset bounds and constraints for anaerobic conditions
LB(j)$(offanaero(j)) = 0;  # Turn off reactions for anaerobic conditions
UB(j)$(offanaero(j)) = 0;
RegulatoryConstraints(j)$(offanaero(j)) = 0;
Solve BiomassMaximization using MIP maximizing Flux('biomass_reaction');

* Record the dimensionality of the solution
ActiveFluxesCount = sum(j, ord(j) * (Flux(j) > 1e-6));  # Count non-zero fluxes

* Display the dimensionality results for anaerobic conditions
Display "Dimensionality of Active Fluxes under Anaerobic Conditions:", ActiveFluxesCount;