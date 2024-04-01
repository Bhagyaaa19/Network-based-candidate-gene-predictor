"""
Bhagya Wijeratne
15219

22/2/2023

Network based candidate gene predictor

"""

# importing the required packages
import networkx as nx
from collections import OrderedDict
import sys
from io import StringIO
from contextlib import redirect_stdout


# Define the capture_output decorator
def capture_output(func):
    def wrapper(*args, **kwargs):
        # Create a StringIO object to capture output
        output_buffer = StringIO()
        with redirect_stdout(output_buffer):
            # Call the original function with arguments and keyword arguments
            func(*args, **kwargs)
        # Retrieve the captured output as a string
        output_value = output_buffer.getvalue()
        return output_value

    return wrapper


class Network:
    n_networks = 0  # number of networks made using the constructor

    def __init__(self, network_file):
        self.hs_dict_multiple = None
        self.known_proteins_in_network = None
        self.known_multiple = None
        self.unknown_multiple = None
        self.all_proteins_set = None
        self.known_set = None
        self.unknown_genes = None
        self.protein_graph = self.networkToGraph(network_file)
        self.n_proteins, self.n_interactions = self.networkInfo()
        Network.n_networks += 1

    # ========================================================================================================
    # Method 1 - ●	A method to create and return a NetworkX graph object when a PPI network file from the STRING database
    # is given as an input in txt format.

    # Input: PPI network file from the STRING database in txt format
    # Output: NetworkX graph object from the data
    # -------------------------------------------------------------------
    @staticmethod  # static method used such that the methods can be called with or without network object
    def networkToGraph(network_file):  # input the file name as a string

        # creating an empty graph first for the imported protein network
        protein_graph = nx.Graph()

        with open(network_file, "r") as file:  # opening the file in read mode
            # skip the first line as it contains headers
            next(file)

            # make a network connecting each node 1 and node 2
            for line in file:  # iterating through the file line by line
                columns = line.split('\t')  # split the lines into columns using the tab seperation
                # extracting the required node 1 and node 2 information
                # converting to uppercase to avoid mismatches
                node1 = columns[0].upper()
                node2 = columns[1].upper()

                # adding edges to the graph
                protein_graph.add_edge(node1, node2)

        # returning the network graph object

        return protein_graph

    # =======================================================================================================================

    # ●	A method to return the number of proteins and interactions in a given PPI network.
    # Input: The protein network graph
    # Output: The number of proteins and interactions
    # ----------------------------------------------------------------------------------------------------------------------
    def networkInfo(self):  # input the networkx graph object as an input
        n_proteins = self.protein_graph.number_of_nodes()
        n_interactions = self.protein_graph.number_of_edges()
        return n_proteins, n_interactions

    # ========================================================================================================================

    # extra method to Extract known proteins from the seed protein list for the particular function.
    # also returns all proteins set, known genes set and unknown genes set
    def extractSeedProteinsSingleFunction(self, known_proteins_file):

        # 1.	Make a set with the names of all the proteins ( all_proteins_set)
        self.all_proteins_set = set(list(self.protein_graph.nodes))

        # Extract known proteins from the seed protein list for the particular function.
        # opening the file
        with open(known_proteins_file, "r") as known_proteins_file:
            self.known_set = set()  # creating an empty set to store the known gene names
            for line in known_proteins_file:
                columns = line.split('\t')  # splitting the lines into columns based on tab space

                # extracting the required protein information from the file
                # to prevent case mismatches, all the data are turned into upper case
                protein = columns[0].upper()
                # assigning the elements to the set
                self.known_set.add(protein)

        # Unknown_genes = all_proteins_set – known_set
        self.unknown_genes = self.all_proteins_set.difference(self.known_set)
        # print("known,",
        #       self.known_set)
        return self.all_proteins_set, self.known_set, self.unknown_genes

    # ===================================================================================================================
    # extra Method to count the number of neighbours known for the particular function ( Single function)
    def countNeighboursSingle(self, gene, known_set):
        neighbors = nx.neighbors(self.protein_graph, gene)  # Returns an iterator over all neighbors of protein
        known_neighbors = list()  # making an empty list to store the known neighbors of the particular protein
        for neighbor in neighbors:  # iterating through all the neighbors
            if neighbor in known_set:  # checking whether the protein is a known one
                known_neighbors.append(neighbor)  # if so append
                neighborCountForFunction = len(
                    known_neighbors)  # setting value of the protein as the number of known proteins interacting

        return neighborCountForFunction

    # ====================================================================================================================
    # ●	A method to predict candidate genes for a single function using the majority voting algorithm.

    # Input:
    # The network graph object containing PPI network information for a particular protein
    # GO annotations for seed protein list for a single function – transcription regulatory activity related genes in Arabidopsis thaliana

    # Output:
    # Predicted candidate genes for a particular function
    # -----------------------------------------------------------------------------------------------------------------------
    # Process:
    @capture_output
    def gene_predict_single_function(self, known_proteins_file):

        # Calling extra method to Extract known proteins from the seed protein list for the particular function.
        # also returns all proteins set, known genes set and unknown genes set 
        self.extractSeedProteinsSingleFunction(known_proteins_file)
        # 1.Makes a set with the names of all the proteins ( all_proteins_set)
        # 2.Extracts known proteins from the seed protein list for the particular function.- known_set
        # 3. Unknown_genes = all_proteins_set – known_set

        # 4.	Make a dictionary to store the unknown protein names with their respective MVS ( mvs_dict)
        mvs_dict = {}
        for gene in self.unknown_genes:
            # a.	Mvs[protein] = number of edges to known proteins in the network ( number of neighbors known to the function)
            # calling method to count the number of neighbours known for the particular function
            mvs_dict[gene] = self.countNeighboursSingle(
                gene, self.known_set)  # setting value of the protein as the number of known proteins interacting

        # 6.	Ordering the mvs_dict according to the value in the descending order of MVS ( orderedMVS)
        ordered_mvs = OrderedDict(reversed(sorted(mvs_dict.items(), key=lambda x: x[1])))

        # 7.	Defining a threshold MVS, above which the candidate genes are chosen.
        # print(known_set)
        # print(all_proteins_set)
        # print(unknown_genes)
        # print("Majority scores:", ordered_mvs)
        candidate_genes = list()  # list to store the predicted candidate genes
        threshold = 1
        for gene, value in ordered_mvs.items():  # selecting genes above a particular threshold value only
            if value > threshold:
                candidate_genes.append(gene)
        print(
            "The candidate genes for a single function using Majority Voting Algorithm– transcription regulatory activity ")
        print(candidate_genes)  # 8.Output the names of the candidate genes
        return

    # =====================================================================================================================

    # Extra method extract the functions and the relevant proteins from the multiple function seed protein list
    # returns all_proteins_set, known_multiple dictionary of for each function as keys and proteins as values, and unknown_multiple dictionary with functions as keys and respective unknown proteins as the values
    def extractSeedProteinsMultipleFunction(self, known_multiple_function_proteins_file):
        # Make a set with the names of all the proteins ( all_proteins_set)
        self.all_proteins_set = set(list(self.protein_graph.nodes))
        self.unknown_multiple = dict()  # dictionary to store unknown proteins

        # Extract known proteins from the seed protein list for each function ( three different sets).
        with open(known_multiple_function_proteins_file, "r") as proteins_file:
            self.known_multiple = dict()  # creating an empty dictionary to store the known gene names with their functions

            next(proteins_file)  # skip the first line
            for line in proteins_file:
                columns = line.split('\t')  # splitting the lines into columns based on tab space

                # extracting the required functions and the required protein information from the file
                # to prevent case mismatches, all the data are turned into upper case
                function = columns[1].upper()
                proteins = columns[7].strip().upper()
                # assigning the elements to the dictionary
                self.known_multiple[function] = list()  # make an empty list for each key in the dictionary
                self.unknown_multiple[function] = list()  # make a empty list for each key function in the dictionary

                for protein in proteins.split(","):
                    self.known_multiple[function].append(
                        protein)  # assigning the proteins to the dictionary with the respective
                    # function

        # print(self.known_multiple)

        for protein in self.all_proteins_set:
            for functionx, known_proteins in self.known_multiple.items():
                if protein not in known_proteins:
                    self.unknown_multiple[functionx].append(
                        protein)  # assigning the unknown proteins to each function, one at a time

        # print("known", self.known_multiple.items())
        # print(self.unknown_multiple.items())

        return self.all_proteins_set, self.known_multiple, self.unknown_multiple

    # ======================================================================================================================

    # ●	A method to predict candidate genes for multiple functions using the majority voting algorithm iteratively, one function at a time.
    #
    # Input:
    # The network graph object containing PPI network information for a particular protein
    # GO annotations for multiple seed protein list for a three different functions
    #  – transcription regulatory activity related genes in Arabidopsis thaliana

    # Output:
    # Predicted candidate genes for a multiple functions, iteratively one function at a time
    # ------------------------------------------------------------------------------------------------------------------------
    @capture_output
    def gene_predict_multiple_function(self, known_multiple_function_proteins_file):
        # Calling Extra method extract the functions and the relevant proteins from the multiple function seed protein list
        # returns all_proteins_set, known_multiple dictionary of for each function as keys and proteins as values,
        # and unknown_multiple dictionary with functions as keys and respective unknown proteins as the values
        self.extractSeedProteinsMultipleFunction(known_multiple_function_proteins_file)

        # 7.	Make a dictionary to store the unknown protein names with their respective MVS for each function
        # ( mvs_dict_function1, mvs_dict_function2, mvs_dict_function3)
        mvs_dict = {}

        for function, proteins in self.unknown_multiple.items():  # to iterate through the items in the unknown proteins
            mvs_dict[function] = {}
            for protein in proteins:
                # a.	mvs_dict[function][protein] = number of edges to known proteins in the network
                count = self.countNeighboursSingle(protein, self.known_multiple[function])
                mvs_dict[function][protein] = count  # setting value of the protein as the number of
                # known proteins interacting
        # print("MVS dict", mvs_dict)

        # 6.Ordering the mvs_dict according to the value in the descending order of MVS ( orderedMVS) – done for all three functions
        ordered_mvs = {k: dict(sorted(v.items(), key=lambda item: item[1], reverse=True)) for k, v in mvs_dict.items()}

        # 7.	Defining a threshold MVS (T = 1), above which the candidate genes are chosen. - done for all three functions

        # print("ordered MVS:", ordered_mvs)
        threshold = 1
        candidate_genes_multiple = dict()  # list to store the predicted candidate genes
        for function, protein_scores in ordered_mvs.items():
            candidate_genes_multiple[function] = list()
            for gene, score in protein_scores.items():
                # print(function, gene, score)
                if score > threshold:
                    candidate_genes_multiple[function].append(gene)
        # 8.	Output the names of the candidate genes separately for all three functions

        print("The candidate genes for multiple functions")
        print(candidate_genes_multiple)
        return

    # ==========================================================================================================================
    """
    •	A method to predict the most accurate function for unknown proteins by considering all seed lists at one time using the majority voting algorithm.
    
    Input: 
    The network graph object containing PPI network information for a particular protein
    GO annotations for multiple seed protein list for three different functions
     – transcription regulatory activity related genes in Arabidopsis thaliana
    
    Output: 
    Predicted most accurate function for unknown proteins based on the majority score
 
===========================================================================================================================
    """

    # Apply the decorator to the method
    @capture_output
    def function_predict_multiple_function(self, known_multiple_function_proteins_file):
        # Calling Extra method extract the functions and the relevant proteins from the multiple function seed protein list
        # returns all_proteins_set, known_multiple dictionary of for each function as keys and proteins as values,
        # and unknown_multiple dictionary with functions as keys and respective unknown proteins as the values
        self.extractSeedProteinsMultipleFunction(known_multiple_function_proteins_file)

        # 7.	Make a dictionary to store the unknown protein names with their respective MVS for each function
        # ( mvs_dict_function1, mvs_dict_function2, mvs_dict_function3)
        mvs_dict = {}

        for function, proteins in self.unknown_multiple.items():  # to iterate through the items in the unknown proteins
            mvs_dict[function] = {}
            for protein in proteins:
                # a.	mvs_dict[function][protein] = number of edges to known proteins in the network
                count = self.countNeighboursSingle(protein, self.known_multiple[function])
                mvs_dict[function][protein] = count  # setting value of the protein as the number of
                # known proteins interacting
        # print("MVS dict", mvs_dict)

        # 6.Ordering the mvs_dict according to the value in the descending order of MVS ( orderedMVS) – done for all three functions
        ordered_mvs = {k: dict(sorted(v.items(), key=lambda item: item[1], reverse=True)) for k, v in mvs_dict.items()}
        # ordering is done to minimise the times the functions have to be updated based on the scores in highest_scores.

        highest_scores = {}  # 7. dictionary to store each unknown protein and the functions for it based on the highest MVS score

        for function, proteins in ordered_mvs.items():  # 8. iterate the functions and the proteins in the ordered dictionary
            for protein, score in proteins.items():  # accessing the proteins and their respective HS
                if protein not in highest_scores or score > highest_scores[protein][0][1]:
                    # if the protein is not already assigned to the highest scores do it,
                    # if assigned, then check if the score for the current function is higher
                    highest_scores[protein] = [[function, score]]  # if so, update the new function and score

                elif score == highest_scores[protein][0][1]:
                    # in case where there is the same multiple HS score for multiple functions, we have to output all functions
                    highest_scores[protein].append([function, score])

        # print("Highest scores:", highest_scores)

        for protein, scores in highest_scores.items():
            print(f"Protein: {protein}, Functions and Respective Majority Voting Scores:")
            for function, score in scores:
                print(f"    Function: {function}, Majority Voting Score: {score}")
        return

    """
    ==============================================================================================================================================
    
    ●	A method to predict candidate genes for a single function using the Hishigaki algorithm.
    
    Input: 
    The network graph object containing PPI network information for a particular protein
    GO annotations for seed protein list for a single function – transcription regulatory activity related genes in Arabidopsis thaliana
    
    Output: 
    Predicted candidate genes for a particular function
    -----------------------------------------------------------------------------------------------------------------------------------
    """

    # -----------------------------------------------------------------------------------------------------------------------
    # Process:
    @capture_output
    def gene_predict_single_function_Hishigaki(self, known_proteins_file):
        # Calling extra method to Extract known proteins from the seed protein list for the particular function.
        # also returns all proteins set, known genes set and unknown genes set
        self.extractSeedProteinsSingleFunction(known_proteins_file)
        # 1.Makes a set with the names of all the proteins ( all_proteins_set)
        # 2.Extracts known proteins from the seed protein list for the particular function.- known_set
        # 3. Unknown_genes = all_proteins_set – known_set

        # updating the list with the known proteins belonging to the network
        known_proteins_in_network = list()  # list to store the known proteins in the network
        for protein in self.known_set:
            if protein in self.all_proteins_set:
                known_proteins_in_network.append(protein)

        freq_function = len(known_proteins_in_network) / len(
            self.all_proteins_set)  # frequency of the protein of the given function in the entire network

        # 4.	Make a dictionary to store the unknown protein names with their Hishigaki scores ( hs_dict)
        hs_dict = {}
        for gene in self.unknown_genes:
            nf = 0  # initializing variable
            # a.	nf = number of proteins with the given function in the neighborhood
            neighbor_count = 0  # counting the number of neighbours
            neighbors = nx.neighbors(self.protein_graph, gene)  # Returns an iterator over all neighbors of protein
            known_neighbors = list()  # making an empty list to store the known neighbors of the particular protein
            for neighbor in neighbors:  # iterating through all the neighbors
                neighbor_count += 1  # incrementing the number of neighbors
                if neighbor in self.known_set:  # checking whether the protein is a known one
                    known_neighbors.append(neighbor)  # if so append
            nf = len(known_neighbors)
            # print("nf", gene, nf)

            ef = freq_function * neighbor_count  # expected frequency
            hs_score = ((nf - ef) ** 2) / ef  # calculating hs
            hs_dict[gene] = hs_score
        # print("ef",gene, ef)
        # print(hs_dict)

        # 6.	Ordering the hs_dict according to the value in the descending order of HS ( orderedHS)
        orderedHS = OrderedDict(reversed(sorted(hs_dict.items(), key=lambda x: x[1])))
        # print("Hishigaki scores Ordered:", orderedHS)

        # 7.	Defining a threshold HS, above which the candidate genes are chosen.
        threshold = 0.5
        candidate_genes = list()  # list to store the predicted candidate genes
        for gene, value in orderedHS.items():
            if value > threshold:
                candidate_genes.append(gene)
        print(
            "The candidate genes for a single function using Hishigaki Algorithm– transcription regulatory activity  ")
        print(candidate_genes)  # 8.Output the names of the candidate genes
        return

    # ==================================================================================================================================
    # Extra method to calculate the Hishigaki scores
    # this extra method was used to reduce code repetition
    def calculateHishigakiScores(self):
        # finding the number of known proteins for each function belonging to the network
        self.known_proteins_in_network = {}  # initializing a dictionary to store the counts

        for function1, known_proteins in self.known_multiple.items():
            self.known_proteins_in_network[function1] = 0
            for protein1 in known_proteins:
                if protein1 in self.all_proteins_set:
                    self.known_proteins_in_network[function1] += 1
        # print("known proteins in network", self.known_proteins_in_network)

        # 7.Make a dictionary to store the unknown protein names with their respective Hishigaki scores for each function
        self.hs_dict_multiple = {}
        neighbor_count = {}  # dictionary to count the number of neighbours
        # nf = number of proteins with the given function in the neighborhood for a particular function
        nf = {}  # initializing dictionary to store the number of proteins with the given function in the neighborhood for a particular function
        freq_function = {}  # initializing dictionary to store the frequencies of functions
        ef = {}  # initializing dictionary to store the expected frequency

        for function, proteins in self.unknown_multiple.items():  # to iterate through the items in the unknown proteins
            self.hs_dict_multiple[
                function] = dict()  # making a dictionary for each function to store the HS for each protein
            neighbor_count[
                function] = dict()  # making nested dictionaries for every function to store the neighbor counts for each proteins
            nf[function] = dict()  # making nested dictionaries for each function to store nf values for each proteins
            # calculating hishigaki score for each protein
            for protein2 in proteins:
                neighbor_count[function][
                    protein2] = 0  # initializing the variable counting the number of neighbors for each protein
                neighbor_count[function][protein2] = 0
                neighbors = nx.neighbors(self.protein_graph,
                                         protein2)  # Returns an iterator over all neighbors of protein
                known_neighbors = list()  # making an empty list to store the known neighbors of the particular protein
                for neighbor in neighbors:  # iterating through all the neighbors
                    neighbor_count[function][protein2] += 1  # incrementing the number of neighbors
                    if neighbor in self.known_multiple[
                        function]:  # checking whether the protein is a known one for that function
                        known_neighbors.append(neighbor)  # if so append
                nf[function][protein2] = len(known_neighbors)

        # print("neighbour count", neighbor_count)
        # print("nf", nf)

        # calculating ef, the expected frequency of the function = frequency of proteins with given function in the entire network * number of neighbors ( neighbor_count)
        for function, proteins in self.unknown_multiple.items():
            # 8. frequency of the protein of the given function in the entire network
            freq_function[function] = self.known_proteins_in_network[function] / len(self.all_proteins_set)
            ef[function] = dict()  # making a nested dictionary to store the ef for each protein for each function

            for protein in proteins:
                ef[function][protein] = freq_function[function] * neighbor_count[function][
                    protein]  # expected frequency

                hs_score = ((nf[function][protein] - ef[function][protein]) ** 2) / ef[function][
                    protein]  # calculating hs
                self.hs_dict_multiple[function][protein] = hs_score  # calculating the hs
        # print("ef", ef)
        # print('Hishigaki Scores unordered: ', hs_dict_multiple)

    # =======================================================================================================================

    """
====================================================================================================================================
    ●	A method to predict candidate genes for multiple functions using the Hishigaki algorithm iteratively, one function at a time.
    
    Input: 
    The network graph object containing PPI network information for a particular protein
    GO annotations for multiple seed protein list for a three different functions
     – transcription regulatory activity related genes in Arabidopsis thaliana
    
    Output: 
    Predicted candidate genes for multiple functions, iteratively one function at a time using Hishigaki Algorthm
==========================================================================================================================  
    """

    @capture_output
    def gene_predict_multiple_function_Hishigaki(self, known_multiple_function_proteins_file):
        # Calling Extra method extract the functions and the relevant proteins from the multiple function seed protein list
        # returns all_proteins_set, known_multiple dictionary of for each function as keys and proteins as values,
        # and unknown_multiple dictionary with functions as keys and respective unknown proteins as the values
        self.extractSeedProteinsMultipleFunction(known_multiple_function_proteins_file)

        # print("Known proteins:", self.known_multiple)
        # print("unknown proteins for each function:", self.unknown_multiple.items())

        # ==========================================Calculating the Hishigaki scores ===============================

        # calling methods to calculate the Hishigaki Scores for each function for each protein
        self.calculateHishigakiScores()  # makes hs_dict_multiple[function][protein]:HS

        # 6.Ordering the hs_dict according to the value in the descending order of MVS ( orderedMVS) – done for all three functions
        orderedHS_multiple = {k: dict(sorted(v.items(), key=lambda item: item[1], reverse=True)) for k, v in
                              self.hs_dict_multiple.items()}
        # print("Hishigaki scores for multiple function predictions:", orderedHS_multiple)

        # 7.	Defining a threshold HS, above which the candidate genes are chosen. - done for all three functions
        threshold = 0.5
        candidate_genes_multiple_HS = dict()  # list to store the predicted candidate genes
        for function3, protein_scores in orderedHS_multiple.items():
            candidate_genes_multiple_HS[function3] = list()
            for gene, score in protein_scores.items():
                # print(function, gene, score)
                if score > threshold:
                    candidate_genes_multiple_HS[function3].append(gene)
        # 8.	Output the names of the candidate genes separately for all three functions

        print("The candidate genes for multiple functions using Hishigaki Algorithm")
        print(candidate_genes_multiple_HS)
        return

    """
    •	A method to predict the most accurate function for unknown proteins by considering all seed lists at one time using the majority voting algorithm.
    
    Input: 
    The network graph object containing PPI network information for a particular protein
    GO annotations for multiple seed protein list for three different functions
     – transcription regulatory activity related genes in Arabidopsis thaliana
    
    Output: 
    Predicted most accurate function for unknown proteins based on the Hiahigaki score
    
    """

    # Apply the decorator to the method
    @capture_output
    def function_predict_multiple_function_Hishigaki(self, known_multiple_function_proteins_file):
        # Calling Extra method extract the functions and the relevant proteins from the multiple function seed protein list
        # returns all_proteins_set, known_multiple dictionary of for each function as keys and proteins as values,
        # and unknown_multiple dictionary with functions as keys and respective unknown proteins as the values
        self.extractSeedProteinsMultipleFunction(known_multiple_function_proteins_file)

        # print("Known proteins:", self.known_multiple)
        # print("unknown proteins for each function:", self.unknown_multiple.items())

        # ==========================================Calculating the Hishigaki scores ===============================

        # calling methods to calculate the Hishigaki Scores for each function for each protein
        self.calculateHishigakiScores()  # makes hs_dict_multiple[function][protein]:HS

        # 6.Ordering the hs_dict according to the value in the descending order of MVS ( orderedMVS) – done for all three functions
        orderedHS_multiple = {k: dict(sorted(v.items(), key=lambda item: item[1], reverse=True)) for k, v in
                              self.hs_dict_multiple.items()}
        # print("Hishigaki scores after ordering:", orderedHS_multiple)
        # ordering is done to minimise the times the functions have to be updated based on the scores in highest_scores.

        highest_scores = {}  # 7. dictionary to store each unknown protein and the functions for it based on the highest Hishigaki score

        for function, proteins in orderedHS_multiple.items():  # 8. iterate the functions and the proteins in the ordered dictionary
            for protein, score in proteins.items():  # accessing the proteins and their respective HS
                if protein not in highest_scores or score > highest_scores[protein][0][1]:
                    # if the protein is not already assigned to the highest scores do it,
                    # if assigned, then check if the score for the current function is higher
                    highest_scores[protein] = [[function, score]]  # if so, update the new function and score

                elif score == highest_scores[protein][0][1]:
                    # in case where there is the same multiple HS score for multiple functions, we have to output all functions
                    highest_scores[protein].append([function, score])

        # print("Highest scores:", highest_scores)

        for protein, scores in highest_scores.items():
            print(f"Protein: {protein}, Functions and Respective Hishigaki Scores:")
            for function, score in scores:
                print(f"    Function: {function}, Hishigaki Score: {score}")

        # 9. Output the names of the unknown proteins with the relevant functions and MVS

        return
