import cobra
import matplotlib.pyplot as plt

from utils import load_model


def analyse_fva(model_path):
    model = load_model(model_path)

    # Set the environmental conditions by changing the bounds of the respective reactions For aerobic conditions,
    # the oxygen uptake rate should be non-zero (e.g., -20), and for anaerobic, it should be zero.
    model.reactions.get_by_id("EX_o2_e").bounds = (-20, 0)  # example for aerobic conditions

    # Change to anaerobic conditions
    model.reactions.get_by_id("EX_o2_e").bounds = (0, 0)
    anaerobic_solution = model.optimize()

    biomass_values = []
    succinate_production = []
    biomass_reaction = model.reactions.get_by_id("BIOMASS_Ec_iAF1260_core_59p81M")
    succinate_reaction = model.reactions.get_by_id("EX_succ_e")

    for percentage in range(0, 110, 10):  # from 0 to 100% in increments of 10%
        biomass_constraint = anaerobic_solution.objective_value * (percentage / 100.0)
        biomass_values.append(biomass_constraint)  # append the current biomass constraint to biomass_values
        model.reactions.get_by_id(biomass_reaction.id).bounds = (biomass_constraint, biomass_constraint)

        # Perform FVA
        fva_result = cobra.flux_analysis.flux_variability_analysis(
            model, reaction_list=[succinate_reaction.id], loopless=True
        )

        # Access the FVA result using proper DataFrame indexing
        min_flux = fva_result.at[succinate_reaction.id, 'minimum']
        max_flux = fva_result.at[succinate_reaction.id, 'maximum']

        succinate_production.append((min_flux, max_flux))

    # Now, biomass_values should have the same length as min_succinate and max_succinate
    min_succinate, max_succinate = zip(*succinate_production)

    # Proceed with plotting
    plt.plot(biomass_values, min_succinate, label='Minimum Succinate Production')
    plt.plot(biomass_values, max_succinate, label='Maximum Succinate Production')
    plt.xlabel('Biomass Constraint')
    plt.ylabel('Succinate Production Flux')
    plt.legend()
    plt.savefig('../results/iAF1260/succinate_production.png', format='png', dpi=300)
    plt.show()

    plt.close()


if __name__ == '__main__':
    path = "../models/iAF1260.xml"
    analyse_fva(path)
