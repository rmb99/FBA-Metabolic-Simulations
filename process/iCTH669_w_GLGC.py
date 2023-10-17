from utils import load_model
from cobra.flux_analysis import flux_variability_analysis


def fba_analysis(path):
    # Load the model
    model = load_model(path)

    # Define a list of conditions to test, e.g., different cellulose uptake rates
    uptake_conditions = [3.0, 3.50, 4.0]

    # Open a file to write the results
    with open("../results/fba_results.txt", "w") as file:
        # Perform FBA under different conditions
        for condition in uptake_conditions:
            # Set the uptake rate
            model.reactions.EXCH_cellb_e.lower_bound = -condition  # the uptake reactions typically have negative fluxes

            # Optimize the model
            solution = model.optimize()

            # Write the results to the file
            file.write(f"Results for cellulose uptake rate of {condition}:\n")
            file.write(f"Objective value (growth rate): {solution.objective_value}\n")
            file.write("Fluxes:\n")
            for reaction_id, flux in solution.fluxes.items():
                file.write(f"{reaction_id}: {flux}\n")
            file.write("\n")


def fva_analysis(path):
    # Load the model
    model = load_model(path)

    # Define a list of conditions to test, e.g., different cellulose uptake rates
    uptake_conditions = [3.0, 3.50, 4.0]

    # Open a file to write the results
    with open("../results/fva_results.txt", "w") as file:
        # Perform FVA under different conditions
        for condition in uptake_conditions:
            # Set the uptake rate
            model.reactions.EXCH_cellb_e.lower_bound = -condition  # the uptake reactions typically have negative fluxes

            # Perform FVA
            fva_result = flux_variability_analysis(model)

            # Write the results to the file
            file.write(f"Results for cellulose uptake rate of {condition}:\n")
            for reaction_id, fluxes in fva_result.iterrows():
                file.write(f"{reaction_id}: Minimum Flux: {fluxes['minimum']}, Maximum Flux: {fluxes['maximum']}\n")
            file.write("\n")


if __name__ == '__main__':
    model_path = "../models/iCTH669_w_GLGC.sbml"
    fba_analysis(model_path)
    fva_analysis(model_path)