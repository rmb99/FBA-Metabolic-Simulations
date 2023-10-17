
from escher import Builder
from utils import load_model


def preform_visualization(model_path):
    model = load_model(model_path)
    solution = model.optimize()
    model.summary()
    # Build an Escher map with the solution
    builder = Builder(
        model=model,
        map_name='e_coli_core.Core metabolism',
        reaction_data=solution.fluxes
    )

    # Save the Escher map as an HTML file
    builder.save_html('../results/iAF1260_metabolic.html')


if __name__ == '__main__':
    models_path = '../models/iAF1260.xml'
    preform_visualization(models_path)
