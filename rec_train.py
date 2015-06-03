# -*- coding: utf-8 -*-

from __future__ import print_function, division
import argparse
import os
import functions
import shutil
import pickle

if __name__ == "__main__":
    #Arguments
    parser = argparse.ArgumentParser(
            description=("Use this script to get a new model of the "
                "rss. The length of the sequences should be 28 or "
                "39."))
    parser.add_argument('-p', '--positive', required=True,
            help='Positive fasta file.')
    parser.add_argument('-n', '--negative', required=True,
            help='Negative fasta file.')
    parser.add_argument('-o', '--output', required=True,
            help='Name of the folder with the model.')
    parser.add_argument('-f', '--features', required=False,
            help='Feature file to use.')
    parser.add_argument('-e', '--extra', required=False,
            help='For RSS 28 the extra nucleotide to use in the end.',
            default='T')

    #Get all the arguments 
    args = parser.parse_args()

    #Getting rec folder
    base_folder = os.path.dirname(os.path.realpath(__file__))

    #Check length of the sequence
    length = functions.check_length(args.positive)
    

    #Define output unique element
    unique = args.output

    #Get extra
    extra = args.extra

    #Check if output folder exists and remove
    if os.path.exists(unique):
        shutil.rmtree(unique)
    os.mkdir(unique)

    #Convert input files to RIC format
    ric_positives = functions.convert_ric(
            filename=args.positive, length=length,
            output="train_ric_positives.fasta", folder=unique)
    ric_negatives = functions.convert_ric(
            filename=args.negative, length=length,
            output="train_ric_negatives.fasta", folder=unique)

    #Aply RIC to the test files
    ric_positives_output = functions.run_ric(test_file=ric_positives,
            length=length, base_folder=base_folder,
            training_file=ric_positives, output="test_ric_positives",
            folder=unique)
    ric_negatives_output = functions.run_ric(test_file=ric_negatives,
            length=length, base_folder=base_folder,
            training_file=ric_positives, output="test_ric_negatives",
            folder=unique)
    os.remove(ric_negatives)
    #Read scores
    ric_positive_scores = functions.read_ric(
            filename=ric_positives_output)
    os.remove(ric_positives_output)
    ric_negative_scores = functions.read_ric(
            filename=ric_negatives_output)
    os.remove(ric_negatives_output)
    ric_value_fdr = functions.get_threshold_unique(
            positive = ric_positive_scores,
            negative = ric_negative_scores,
            optim = "FDR")
    ric_value_sen = functions.get_threshold_unique(
            positive = ric_positive_scores,
            negative = ric_negative_scores,
            optim = "SEN")

    
    #Do CRoSSeD
    #Check if feature file was provided and if it's ok
    if args.features != None:
        features = functions.read_feature_file(args.features)
    else:
        features = functions.read_feature_file(os.path.join(
            base_folder, "Data/crossed_features_" + str(length) + 
            ".txt"))
    #Convert input files to CRoSSeD format
    crossed_positives = functions.convert_crossed(
            filename=args.positive,
            length=length,
            output="train_crossed_positives.fasta",
            features=features, positive=True, folder=unique, 
            extra=extra)
    crossed_negatives = functions.convert_crossed(
            filename=args.negative,
            length=length,
            output="train_crossed_negatives.fasta",
            features=features, positive=False, folder=unique, 
            extra=extra)
    #Train CRoSSeD model
    model = functions.train_crossed(positive=crossed_positives,
            negative=crossed_negatives, model="CRoSSeD_model",
            unique=unique, base_folder=base_folder)
    #Aply CRoSSeD to input data
    crossed_positives_output = functions.run_crossed(
            filename=crossed_positives, model=model,
            output="test_crossed_positives", base_folder=base_folder,
            folder=unique)
    crossed_negatives_output = functions.run_crossed(
            filename=crossed_negatives, model=model,
            output="test_crossed_negatives", base_folder=base_folder,
            folder=unique)
    #Read scores
    crossed_positive_scores = functions.read_crossed(
            filename=crossed_positives_output, original=args.positive,
            length=length)
    crossed_negative_scores = functions.read_crossed(
            filename=crossed_negatives_output, original=args.negative,
            length=length)
    crossed_value_fdr = functions.get_threshold_unique(
            positive = crossed_positive_scores, 
            negative = crossed_negative_scores, 
            optim = "FDR")
    crossed_value_sen = functions.get_threshold_unique(
            positive = crossed_positive_scores, 
            negative = crossed_negative_scores, 
            optim = "SEN")

    """
    #Do REC CRoSSeD
    #Check if feature file was provided and if it's ok
    if args.features != None:
        rec_features = functions.read_feature_file(args.features)
    else:
        rec_features = functions.read_feature_file(os.path.join(
            base_folder, "Data/rec_features_" + str(length) + 
            ".txt"))
    #Convert input files to CRoSSeD format
    crossed_positives = functions.convert_crossed(
            filename=args.positive,
            length=length,
            output="train_crossed_positives.fasta",
            features=rec_features, positive=True, folder=unique, 
            extra=extra)
    crossed_negatives = functions.convert_crossed(
            filename=args.negative,
            length=length,
            output="train_crossed_negatives.fasta",
            features=rec_features, positive=False, folder=unique, 
            extra=extra)
    #Train CRoSSeD model
    rec_model = functions.train_crossed(positive=crossed_positives,
            negative=crossed_negatives, model="rec_CRoSSeD_model",
            unique=unique, base_folder=base_folder)
    #Aply CRoSSeD to input data
    crossed_positives_output = functions.run_crossed(
            filename=crossed_positives, model=rec_model,
            output="test_crossed_positives_rec", 
            base_folder=base_folder,
            folder=unique)
    crossed_negatives_output = functions.run_crossed(
            filename=crossed_negatives, model=rec_model,
            output="test_crossed_negatives_rec", 
            base_folder=base_folder,
            folder=unique)

    #Read scores
    crossed_positive_scores = functions.read_crossed(
            filename=crossed_positives_output, original=args.positive,
            length=length)
    crossed_negative_scores = functions.read_crossed(
            filename=crossed_negatives_output, original=args.negative,
            length=length)
    """
    rec_model = model
    rec_features = features
   

    
    #Combine RIC and CROSSED scores to store in a dictionary.

    positive_scores = functions.combine_scores(filename=args.positive,
            ric=ric_positive_scores, crossed=crossed_positive_scores)
    negative_scores = functions.combine_scores(filename=args.negative,
            ric=ric_negative_scores, crossed=crossed_negative_scores)

    #Compute thresholds
    rec_ric_value_fdr, rec_crossed_value_fdr = functions.get_threshold(
                positive_scores, negative_scores, "BP",
                "FDR")
    rec_ric_value_sen, rec_crossed_value_sen = functions.get_threshold(
                positive_scores, negative_scores, "BP",
                "SEN")


    #Save all the data in the model folder
    data = {
            'ric_template': "train_ric_positives.fasta",
            'rec_crossed_model': rec_model, 
            'rec_crossed_features': rec_features,'extra': extra, 
            'crossed_model':model, 'crossed_features':features}
    data['ric_value_fdr'] = ric_value_fdr
    data['ric_value_sen'] = ric_value_sen
    data['crossed_value_fdr'] = crossed_value_fdr
    data['crossed_value_sen'] = crossed_value_sen
    data['rec_ric_value_fdr'] = rec_ric_value_fdr
    data['rec_ric_value_sen'] = rec_ric_value_sen
    data['rec_crossed_value_fdr'] = rec_crossed_value_fdr
    data['rec_crossed_value_sen'] = rec_crossed_value_sen
    data['max_dist'] = {
            'rric':max(ric_positive_scores) - min(ric_negative_scores),
            'crossed':max(crossed_positive_scores) - min(
                crossed_negative_scores)}
    data['min_values'] = {'rric':min(ric_negative_scores), 
            'crossed': min(crossed_negative_scores)}
    pickle.dump(data, open(os.path.join(unique, "data"), 'w'))





#a vim: ai sts=4 et sw=4
