import cobra


def load_model(model_path):
    # Load the model
    model = cobra.io.read_sbml_model(model_path)
    return model
