# -*- coding: utf-8 -*-
from __future__ import print_function, division
from Bio import SeqIO
import os
import subprocess
import math


classic_threshold = {28: -38.81, 39: -58.45}

def read_feature_file(filename):
    """
    Read a feature file and return a list of dictionaries
    to convert sequences in the numerical features.
    Format:
    ID\tDinucleotide\t...\tFeature Name
    """
    features = []
    try:
        with open(filename) as fh:
            for line in fh:
                temp = line.strip().split("\t")
                if temp[0] == "ID":
                    header = temp
                else:
                    feature = dict()
                    for x in xrange(2, len(temp)):
                        feature[header[x].upper()] = temp[x]
                    features.append(feature)
        return features
    except:
        print("Error with the feature file")
        raise Exception


def check_length(filename):
    #Check if the length of the sequences in the file is unique and 
    #is 28 or 39
    try:
        for record in SeqIO.parse(filename, 'fasta'):
            length = len(record.seq)
            if length != 28 and length != 39:
                print("The sequences do not have 28 or 39 of length")
                raise Exception
            else:
                try:
                    if length != size:
                        print("Sequences don't have the same length")
                        raise Exception
                except NameError:
                    size = length
        return size
    except:
        print("Problem reading the positive file. Maybe not fasta?")
        raise Exception


def convert_ric(filename, length, folder, output):
    #Convert file to RIC format
    out = open(os.path.join(folder, output), 'w')
    for record in SeqIO.parse(filename, 'fasta'):
        if len(record.seq) == length:
            record.seq = record.seq.lower()
            SeqIO.write(record, out, 'fasta')
        else:
            print("Sequence " + record.id + 
                    " do not have the rigth length")
            raise Exception
    out.close()
    return os.path.abspath(os.path.join(folder, output))


def run_ric(test_file, length, base_folder, training_file, output, 
        folder):
    #Run RIC useing taining_file as input for training RIC and
    #test_file as input for RIC to test
    output = os.path.join(folder, output)
    if os.path.exists(output):
        os.remove(output)
    script_file = os.path.join(os.path.join(base_folder, "Scripts"),
            "ComputeRIC") + str(length-16) + ".pl"
    try:
        process = subprocess.Popen("perl " + script_file + " " + test_file +
                " " + output + " " + training_file, shell=True).wait()
        if process == 0:
            return output
        else:
            print("Error with RIC computation")
            raise Exception
    except:
        print("Error calling RIC")
        raise Exception


def convert_crossed(filename, length, output, features, positive, 
        folder, extra = "T"):
    out = output
    #Convert a fasta file into CRoSSeD input format
    output = open(os.path.join(folder, out), 'w')
    for record in SeqIO.parse(filename, 'fasta'):
        seq = record.seq.upper().tostring()
        if len(seq) == length:
            if len(seq) == 28:
                seq = seq + extra
            if positive:
                output.write("1\n")
            else:
                output.write("0\n")
            output.write(seq[0])
            for feature in features:
                output.write("\t+")
            output.write("\n")
            for pos in xrange(1, len(seq)):
                output.write(seq[pos])
                for feature in features:
                    output.write("\t" + feature[seq[pos-1:pos + 1]])
                output.write("\n")
            output.write("\n")
        else:
            print("Sequence " + record.id + 
                    " do not have the rigth length")
            raise Exception
    output.close()
    return os.path.join(folder, out)


def train_crossed(positive, negative, model, unique, base_folder):
    #Train a new CRoSSeD model
    model = os.path.join(unique, model)
    training_file = os.path.join(unique, "train_crossed")
    output = open(training_file, 'w')
    with open(positive) as fh:
        for line in fh:
            output.write(line)
    with open(negative) as fh:
        for line in fh:
            output.write(line)
    output.close()
    script_file = os.path.join(os.path.join(base_folder, "Scripts"),
            "crossed_train.pl")
    try:
        process = subprocess.Popen("perl " + script_file + 
                " -corr 0 -tf " + training_file + " -mf " + model + 
                " -thr 0.999999 > temp",
                shell=True).wait()
        if process == 0:
            os.remove(training_file)
            os.remove("temp")
            return os.path.abspath(model)
        else:
            print("Error with CRoSSeD computation")
            raise Exception
    except:
        print("Error calling CRoSSeD train")
        raise Exception


def run_crossed(filename, model, output, base_folder, folder):
    #Apply crossed to filename
    output = os.path.join(folder, output)
    script_file = os.path.join(os.path.join(base_folder, "Scripts"),
            "crossed_test.pl")
    try:
        process = subprocess.Popen("perl " + script_file + 
                " -corr 0 -tf " + filename + " -mf " + model + " -out " 
                + output + " > temp ",
                shell=True).wait()
        if process == 0:
            os.remove("temp")
            return output
        else:
            print("Error with CRoSSeD computation")
            raise Exception
    except:
        print("Error calling CRoSSeD test")
        raise Exception


def read_ric(filename):
    scores = []
    with open(filename) as fh:
        for line in fh:
            name, seq, score = line.strip().split("\t")
            scores.append(float(score))
    return scores


def read_crossed(filename, original, length):
    scores = []
    data = dict()
    with open(filename, 'r') as fh:
        fh.readline()
        for record in SeqIO.parse(original, 'fasta'):
            rss = record.seq.upper().tostring()
            if rss in data:
                scores.append(float(fh.readline().strip().split("\t"
                    )[-1].split("/")[1].strip()))
                data.pop(rss)
            else:
                while rss not in data:
                    line = fh.readline()
                    seq = "".join(line.strip().
                            split("\t")[0:length]).replace("|", "")
                    score = float(line.strip().split("\t")[-1].split(
                        "/")[1].strip())
                    data[seq] = score
                scores.append(data[rss])
                data.pop(rss)
    os.remove(filename)
    return scores


def combine_scores(**kwargs):
    filename = kwargs['filename']
    ric = kwargs['ric']
    crossed = kwargs['crossed']
    scores = dict()
    pos = 0
    for record in SeqIO.parse(filename, 'fasta'):
        scores[record.id] = [ric[pos], crossed[pos]]
        pos += 1
    return scores


def get_overlap(pos, neg):
    overlap = []
    for score in pos:
        if score <= neg[-1]:
            overlap.append(score)
        else:
            overlap.append(score)
            break
    return overlap


def compare(**kwargs):

    ric = kwargs['ric']
    crossed = kwargs['crossed']
    ric_score = kwargs['ric_score']
    crossed_score = kwargs['crossed_score']
    
    if ric_score >= ric and crossed_score >= crossed:
        return True
    else:
        return False

def get_matrix(**kwargs):
    matrix = {'tp': 0, 'tn': 0, 'fp': 0, 'fn': 0}
    positive = kwargs['positive']
    method = kwargs['method']
    for rss in positive:
        pos = compare(method=method, ric_score=positive[rss][0],
                crossed_score=positive[rss][1], ric=kwargs['ric'],
                crossed=kwargs['crossed'])
        if pos:
            matrix['tp'] += 1
        else:
            matrix['fn'] += 1
    negative = kwargs['negative']
    for rss in negative:
        pos = compare(method=method, ric_score=negative[rss][0],
                crossed_score=negative[rss][1], ric=kwargs['ric'],
                crossed=kwargs['crossed'])
        if pos:
            matrix['fp'] += 1
        else:
            matrix['tn'] += 1
    return matrix


def mcc(matrix):
    tp = matrix['tp']
    tn = matrix['tn']
    fp = matrix['fp']
    fn = matrix['fn']
    if (tp + fp) != 0 and (tp + fn) != 0 and (tn + fp) != 0 and (tn + fn) != 0:
        return (((tp * tn) - (fp * fn)) / (math.sqrt((tp + fp) * (tp + fn) * (
            tn + fp) * (tn + fn))))
    else:
        return 0


def fdr(matrix):
    try:
        return 1 - (matrix['fp'] / (matrix['fp'] + matrix['tp']))
    except:
        return 0


def sen(matrix):
    try:
        return matrix['tp'] / (matrix['tp'] + matrix['fn'])
    except:
        return 0


statistics = {'MCC': mcc, 'FDR': fdr, 'SEN': sen}
default_order = ['FDR', 'MCC', 'SEN']


def get_matrix_rec(**kwargs):
    matrix = {'tp': 0, 'tn': 0, 'fp': 0, 'fn': 0}
    positive = kwargs['positive']
    for score in positive:
        if score >= kwargs['rec']:
            matrix['tp'] += 1
        else:
            matrix['fn'] += 1
    negative = kwargs['negative']
    for score in negative:
        if score >= kwargs['rec']:
            matrix['fp'] += 1
        else:
            matrix['tn'] += 1
    return matrix


def get_threshold_unique(positive, negative, optim):
    if optim == "FDR":
        max_neg = max(negative)
        positive.sort()
        for value in positive:
            if value > max_neg:
                return value
    elif optim == "SEN":
        return min(positive)
def get_false_positives(ric_value, crossed_value, negative):
    fp = 0
    for rss in negative:
        ric = negative[rss][0]
        crossed = negative[rss][1]
        if ric >= ric_value and crossed >= crossed_value:
            fp += 1
    return fp

def get_true_positives(ric_value, crossed_value, positive):
    tp = 0
    for rss in positive:
        ric = positive[rss][0]
        crossed = positive[rss][1]
        if ric >= ric_value and crossed >= crossed_value:
            tp += 1
    return tp


def get_threshold(positive, negative, method, optim):
    ric_pos = list(zip(*positive.values())[0])

    crossed_pos = list(zip(*positive.values())[1])

    best_fp = 10000
    best_tp = 0
    if optim == "FDR":
        for ric_value in ric_pos:
            for crossed_value in crossed_pos:
                fp = get_false_positives(ric_value, crossed_value, 
                        negative)
                tp = get_true_positives(ric_value, crossed_value, 
                        positive)
                if fp < best_fp or (fp == best_fp and tp > best_tp):
                    best_fp = fp
                    best_tp = tp
                    best_ric = ric_value
                    best_crossed = crossed_value
        return best_ric, best_crossed
    elif optim == "SEN":
        for ric_value in ric_pos:
            for crossed_value in crossed_pos:
                fp = get_false_positives(ric_value, crossed_value, 
                        negative)
                tp = get_true_positives(ric_value, crossed_value, 
                        positive)
                if tp > best_tp or (tp == best_tp and fp < best_fp):
                    best_fp = fp
                    best_tp = tp
                    best_ric = ric_value
                    best_crossed = crossed_value
        return best_ric, best_crossed    
                




def classify(algorithm, crossed_scores, ric_scores, rric_scores, data,
        length):
    #Generate a dict of list with the YES/NO value for sequence
    result = dict()
    if algorithm == "RIC" or algorithm == "REC":
        temp = []
        for score in ric_scores:
            if score >= classic_threshold[length]:
                temp.append(True)
            else:
                temp.append(False)
        result['RIC'] = temp
    if algorithm == "rRIC":
        temp_fdr = []
        temp_sen = []
        for score in rric_scores:
            if score >= data['ric_value_sen']:
                temp_sen.append(True)
                if score >= data['ric_value_fdr']:
                    temp_fdr.append(True)
                else:
                    temp_fdr.append(False)
            else:
                temp_sen.append(False)
                temp_fdr.append(False)
        result['rRIC_fdr'] = temp_fdr
        result['rRIC_sen'] = temp_sen
    if algorithm == "CROSSED":
        temp_fdr = []
        temp_sen = []
        for score in crossed_scores:
            if score >= data['crossed_value_sen']:
                temp_sen.append(True)
                if score >= data['crossed_value_fdr']:
                    temp_fdr.append(True)
                else:
                    temp_fdr.append(False)
            else:
                temp_sen.append(False)
                temp_fdr.append(False)
        result['CROSSED_fdr'] = temp_fdr
        result['CROSSED_sen'] = temp_sen
    if algorithm == "REC":
        temp_fdr = []
        temp_sen = []
        temp_rec = []
        for pos in xrange(len(rric_scores)):
            r = rric_scores[pos]
            c = crossed_scores[pos]
            if r >= data['rec_ric_value_sen'] and c >= data[
                    'rec_crossed_value_sen']:
                temp_sen.append(True)
            else:
                temp_sen.append(False)
            if r >= data['rec_ric_value_fdr'] and c >= data[
                    'rec_crossed_value_fdr']:
                temp_fdr.append(True)
            else:
                temp_fdr.append(False)
            dr = r - data['min_values']['rric']
            rt = dr / data['max_dist']['rric']
            dc = c - data['min_values']['crossed']
            ct = dc / data['max_dist']['crossed']
            temp_rec.append((rt + ct)/2)
        result['REC_fdr'] = temp_fdr
        result['REC_sen'] = temp_sen
        result['REC'] = temp_rec
    return result



   

def save(output, result, algorithm, original, ric_scores, rric_scores,
        crossed_scores):
    output = open(output, 'w')
    pos = 0 
    if algorithm == "REC":
        output.write("Sequence ID\tRIC score\trRIC score\tCRoSSeD "
                "score\tREC score\tPASS/FAIL RIC\tPASS/FAIL REC "
                "stringent\tPASS/FAIL REC sensitive\n")
    elif algorithm == "RIC":
        output.write("Sequence ID\tScore\tPASS/FAIL\n")
    else:
        output.write("Sequence ID\tScore\tPASS/FAIL stringent\t "
        "PASS/FAIL sensitive\n")
    for record in SeqIO.parse(original, 'fasta'):
        id = record.id
        output.write(id + "\t")
        if algorithm == "REC":
            output.write(str(ric_scores[pos]) + "\t")
            output.write(str(rric_scores[pos]) + "\t")
            output.write(str(crossed_scores[pos]) + "\t")
            output.write(str(result['REC'][pos]) + "\t")
            if result['RIC'][pos]:
                output.write("1\t")
            else:
                output.write("0\t")
            if result['REC_fdr'][pos]:
                output.write("1\t")
            else:
                output.write("0\t")
            if result['REC_sen'][pos]:
                output.write("1\n")
            else:
                output.write("0\n")
            
        elif algorithm == "RIC":
            output.write(str(ric_scores[pos]) + "\t")
            if result['RIC'][pos]:
                output.write("1\n")
            else:
                output.write("0\n")
        elif algorithm == "rRIC":
            output.write(str(rric_scores[pos]) + "\t")
            if result['rRIC_fdr'][pos]:
                output.write("1\t")
            else:
                output.write("0\t")
            if result['rRIC_sen'][pos]:
                output.write("1\n")
            else:
                output.write("0\n")
        elif algorithm == "CROSSED":
            output.write(str(crossed_scores[pos]) + "\t")
            if result['CROSSED_fdr'][pos]:
                output.write("1\t")
            else:
                output.write("0\t")
            if result['CROSSED_sen'][pos]:
                output.write("1\n")
            else:
                output.write("0\n")
            
        pos += 1
    
    output.close()


if __name__ == "__main__":
    pass


# vim: ai sts=4 et sw=4
