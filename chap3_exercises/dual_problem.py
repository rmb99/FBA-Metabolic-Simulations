import numpy as np
from scipy.optimize import linprog
import cobra

from utils import load_model


def calculate_primal_problem(model_path):
    model = load_model(model_path)
    # Set environmental conditions, e.g., anaerobic
    model.reactions.get_by_id("EX_o2_e").lower_bound = -1

    # Optimize the primal problem (Maximize biomass)
    biomass_reaction = model.reactions.get_by_id("BIOMASS_Ec_iAF1260_core_59p81M")
    model.objective = biomass_reaction
    primal_solution = model.optimize()

    return primal_solution


def calculate_dual_solution(model_path):
    model = load_model(model_path)

    # Create the stoichiometric matrix
    stoichiometric_matrix = cobra.util.array.create_stoichiometric_matrix(model)

    # The right-hand side of the equality constraints, b, is zero for steady-state FBA (S*v = 0)
    b = np.zeros(len(model.metabolites))

    # The c vector for the primal problem is the coefficients of the objective function (usually 1 for the biomass
    # reaction and 0 for others)
    biomass_reaction = model.reactions.get_by_id("BIOMASS_Ec_iAF1260_core_59p81M")
    c = np.zeros(len(model.reactions))
    c[model.reactions.index(biomass_reaction)] = 1

    # For the dual problem, the objective coefficients are the b vector from the primal
    c_dual = b

    # The inequality constraints in the dual are derived from the stoichiometric matrix and the c vector from the primal
    a_ub_dual = np.transpose(stoichiometric_matrix)
    b_ub_dual = c

    # The signs are flipped because the primal is a maximization problem, and the dual is minimization
    c_dual = -c_dual
    b_ub_dual = -b_ub_dual

    # Solve the dual problem
    res = linprog(c_dual, A_ub=a_ub_dual, b_ub=b_ub_dual, method='highs')

    # Check if the solutions are close enough (considering some tolerance due to numerical optimization precision)
    primal_solution = calculate_primal_problem(model_path)
    tolerance = 1e-6
    if np.isclose(primal_solution.objective_value, -res.fun, atol=tolerance):
        print(f"Dual solution: {-res.fun}")
        print("Primal and dual solutions match!")
    else:
        print(f"Dual solution: {-res.fun}")
        print("Primal and dual solutions do not match within the specified tolerance.")


if __name__ == '__main__':
    path = "../models/iAF1260.xml"
    calculate_dual_solution(path)
