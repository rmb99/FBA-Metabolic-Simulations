
class Fba:
    def __init__(self, model, biomass_reaction):
        self.model = model
        self.biomass_reaction = biomass_reaction

    def maximize(self):
        self.model.objective = self.biomass_reaction
        primal_solution = self.model.optimize()
        return primal_solution
