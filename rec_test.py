# -*- coding: utf-8 -*-

from __future__ import print_function, division
import argparse
import os
import functions
import pickle


classic_threshold = {28: -38.81, 39: -58.45}

if __name__ == "__main__":
    #Arguments
    parser = argparse.ArgumentParser(
            description=("Use this script to test a set of sequences"))
    parser.add_argument('-t', '--test', required=True,
            help=("Test file. All sequences should have the same length"
            " and it should be 28 or 39."))
    parser.add_argument('-o', '--output', required=True,
            help='Name of the output file.')
    parser.add_argument('-m', '--model', required=False, help=("Name of"
    "the folder with the model."))
    parser.add_argument('--algorithm', required=False, help=("The "
    "algorith used to classify the RSS. By default all are used."),
        choices=('RIC', 'rRIC', 'CROSSED', 'REC'),
        default='REC')

    

    args = parser.parse_args()
    algorithm = args.algorithm
    #Getting rec folder
    base_folder = os.path.dirname(os.path.realpath(__file__))

    #Check length of the sequence
    length = functions.check_length(args.test)

    #Get model data
    if args.model == None:
        model = os.path.join(os.path.join(base_folder, "Data"), 
                "base_model_" + str(length))
    else:
        model = args.model
    data = pickle.load(open(os.path.join(model, "data")))



    #Tranform data
    if algorithm == "CROSSED":
        crossed_test = functions.convert_crossed(filename = args.test,
                length = length, output = "test_crossed_" + args.output,
                features = data['crossed_features'], positive = True,
                folder = "", extra = data['extra'])
    else:
        if algorithm == "REC":
            crossed_test = functions.convert_crossed(
                    filename = args.test, length = length, 
                    output = "test_crossed_" + args.output, 
                    features = data['rec_crossed_features'], 
                    positive = True, folder = "", extra = data['extra'])
        ric_test = functions.convert_ric(filename = args.test, 
                length = length, output = "test_ric_" + args.output,
                folder = "")

    #Apply algorithms
    if algorithm == "CROSSED":
        crossed_output =  functions.run_crossed(filename=crossed_test,
                model=os.path.join(model, data['crossed_model']),
                output="test_crossed_results_" + args.output,
                base_folder=base_folder, folder="")
    else:
        if algorithm != "rRIC":
            ric_output = functions.run_ric(test_file = ric_test,
                    length = length, base_folder = base_folder, 
                    training_file = os.path.join(os.path.join(
                        base_folder, "Data"), "classic_mouse_" + str(
                        length-16) + ".fasta"), output = 
                        "test_ric_results_" + args.output, folder="")
        if algorithm != "RIC":
            rric_output = functions.run_ric(test_file = ric_test, 
                    length = length, base_folder = base_folder, 
                    training_file = os.path.join(model, 
                        data['ric_template']), output = 
                    "test_rric_results_" + args.output, folder = "")
        if algorithm == "REC":
            crossed_output =  functions.run_crossed(
                    filename = crossed_test, 
                    model = os.path.join(model, 
                        data['rec_crossed_model']), 
                    output = "test_crossed_results_" + args.output, 
                    base_folder = base_folder, folder="")

    #Read scores
    if algorithm == "CROSSED" or algorithm == "REC":
        crossed_scores = functions.read_crossed( 
                filename = crossed_output, original = args.test, 
                length = length)
        ric_scores = None
        rric_scores = None
    else:
        crossed_scores = None
    if algorithm != "CROSSED":
        if algorithm != "rRIC":
            ric_scores = functions.read_ric(ric_output)
        else:
            ric_scores = None
        if algorithm != "RIC":
            rric_scores = functions.read_ric(rric_output)
        else:
            rric_scores = None

    #Classify sequences
    result = functions.classify(algorithm = algorithm, 
            crossed_scores = crossed_scores, ric_scores = ric_scores,
            rric_scores = rric_scores, data = data, length = length)

    #Save result
    functions.save(output = args.output, result = result, 
            algorithm = algorithm, original = args.test, 
            ric_scores = ric_scores, rric_scores = rric_scores, 
            crossed_scores = crossed_scores)
    try:
        os.remove("test_crossed_" +args.output)
    except:
        pass
    try:
        os.remove("test_ric_" + args.output)
    except:
        pass
    try:
        os.remove("test_ric_results_" + args.output)
    except:
        pass
    try:
        os.remove("test_rric_results_" + args.output)
    except:
        pass





        

    
# vim: ai sts=4 et sw=
