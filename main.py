"""
Bhagya Wijeratne
15219

22/2/2023

Network based candidate gene predictor - main script

"""
import networkx as nx

import script

# defining the main method to call and implement methods
if __name__ == '__main__':
    print("Creating a Graph for the protein interaction file.....")
    # storing the networkx graph object in a variable
    protein_graph = script.Network.networkToGraph("string_interactions_short_prototype2.tsv")
    # converting to a file such that cytoscape can open ( to confirm creation of network graph)
    print("Creating gml file such that cytoscape can open....")
    nx.write_gml(protein_graph, "proteinNetwork.gml")
    print('\n')

    print("=================Basic Information===================")
    print("Creating a Network object for the provided Protein Interaction file.......")
    # create a Network object to call the rest of the functions
    network1 = script.Network("string_interactions_short_prototype2.tsv")
    # storing the number of proteins and interactions in two variables, calling the method defined in script
    n_proteins, n_interactions = network1.networkInfo()
    # printing the output to console
    print("Number of proteins in the network:", n_proteins)
    print("Number of interactions in the network:", n_interactions)
    print('\n')

    print("=================Using Majority Voting Algorithm===================")
    # calling the method to predict the candidate genes for the single function - transcription regulatory activity
    print(
        "Predicting candidate genes for single function, transcription regulatory activity using Majority Voting...... ")
    outputMVSingle = network1.gene_predict_single_function("GO_annotations_test_single_all.txt")
    print(outputMVSingle)
    print('\n')

    print("\x1B[4m" + "Predicting candidate genes for three functions,using Majority Voting...... " + "\x1B[0m")
    outputMVMultiple = network1.gene_predict_multiple_function("enrichment_Function.tsv")
    print (outputMVMultiple)
    print('\n')

    print(
        "\x1B[4m" + "Predicting most accurate function for unknown proteins,using Majority Voting Scores...... " + "\x1B[0m")
    outputMVFunction = network1.function_predict_multiple_function("enrichment_Function.tsv")
    print(outputMVFunction)
    print('\n')

    print("=================Using Hishigaki Algorithm===================")
    # calling the method to predict the candidate genes for the single function - transcription regulatory activity
    print(
        "Predicting candidate genes for single function, transcription regulatory activity using Hishigaki Algorithm...... ")
    outputHSSingle = network1.gene_predict_single_function_Hishigaki("GO_annotations_test_single_all.txt")
    print(outputHSSingle)
    print('\n')

    print("\x1B[4m" + "Predicting candidate genes for three functions,using Hishigaki Algorithm...... " + "\x1B[0m")
    outputHSMultiple = network1.gene_predict_multiple_function_Hishigaki("enrichment_Function.tsv")
    print(outputHSMultiple)
    print('\n')

    print(
        "\x1B[4m" + "Predicting most accurate function for unknown proteins,using Hishigaki Scores...... " + "\x1B[0m")
    outputHSFunction = network1.function_predict_multiple_function_Hishigaki("enrichment_Function.tsv")
    print(outputHSFunction)
    print('\n')
    print("End of program...")
