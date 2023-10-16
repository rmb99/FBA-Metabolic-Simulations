# -*- coding: utf-8 -*-

import cobra
import numpy as np


def calculate_aero_rank():
    # Load the iAF1260 model, using xml format of iAF1260 model from ref[16]
    model = cobra.io.read_sbml_model("models/iAF1260.xml")

    # Set the model to aerobic conditions and maximize biomass:
    model.reactions.get_by_id("EX_o2_e").lower_bound = -1000

    biomass_reaction = model.reactions.get_by_id("BIOMASS_Ec_iAF1260_core_59p81M")
    model.objective = biomass_reaction
    solution_aerobic = model.optimize()

    # Calculate null space dimensionality under aerobic condition
    stoichiometry_matrix_aerobic = cobra.util.create_stoichiometric_matrix(model, array_type="dense")
    _, s1, _ = np.linalg.svd(stoichiometry_matrix_aerobic)
    rank_aerobic = np.sum(s1 > 1e-6)
    nullity_aerobic = stoichiometry_matrix_aerobic.shape[1] - rank_aerobic
    return nullity_aerobic


def calculate_anaero_rank():
    model_anaero = cobra.io.read_sbml_model("models/iAF1260.xml")

    #  Set the model to anaerobic conditions and maximize biomass
    model_anaero.reactions.get_by_id("EX_o2_e").lower_bound = 0
    model_anaero.reactions.get_by_id("EX_o2_e").upper_bound = 0

    solution_anaerobic = model_anaero.optimize()

    # calculate null space dimentionality under anaerobic condition
    stoich_matrix_anaerobic = cobra.util.create_stoichiometric_matrix(model_anaero, array_type="dense")
    _, s, _ = np.linalg.svd(stoich_matrix_anaerobic)
    rank_anaerobic = np.sum(s > 1e-6)
    nullity_anaerobic = stoich_matrix_anaerobic.shape[1] - rank_anaerobic
    return nullity_anaerobic


if __name__ == '__main__':
    rank_nullity_aerobic = calculate_aero_rank()
    rank_nullity_anerobic = calculate_anaero_rank()

    print(f"Dimensionality of the subspace under aerobic conditions: {rank_nullity_aerobic}")
    print(f"Dimensionality of the subspace under anaerobic conditions: {rank_nullity_anerobic}")