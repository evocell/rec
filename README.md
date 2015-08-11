#########
#PURPOSE#
#########
Bioinformatics tool to test 28-bp and 39-bp DNA sequences for there ability to be targets of the V(D)J recombinase (RAG) and to map potential recombination signal sequences (RSS) in all jawed vertebrates. This tool also enables to produce classifiers with alternative training sets.

#############
#INSTALATION#
#############
Copy the folder on your computer.

Help available in the read-me file.

##############
#REQUIREMENTS#
##############
python-2.7

biopython-1.57

REC uses 2 other tools in order to classity sequences as RSS or not. 

The specific requirements for these tools should also be meet.

 * RIC (1): 
  - perl-5.14
 * CRoSSeD (2): 
  - Only tested under Linux system.
  - Perl-5.8
  - CRF++ toolkit ( http://crfpp.sourceforge.net/)

#######
#USAGE#
#######
To test a set of sequences you should use the script rec_test.py

Options:

-h --help -> shows the help message for this script

-t --test -> Input file to be tested. All sequences should have length equal to either 28 or 39bp, and all sequences of the same file should have the same length.

-o --output -> output file to write the result. Note that if a file with the same name already exists it will be replaced.

-m --model -> If you have trained the REC with another model, give the path to that model. Otherwise our default trainig set and model will be used.

--algorithm -> You can select which algorithm will be used to test you sequences. By default the REC option is selected and in the output all algorithms (RIC, iRIC, CRoSSeD and REC) will be represented. If you select a specific algorithm, only this one will be in the output file.



example:

python rec_test.py -t input.fasta -o output.txt

with new model

python rec_test.py -t input.fasta -o output.txt -m model_folder 



OUTPUT file:

The default order of data in the output file is:

Sequence ID -> The id you provide for each sequence in the input fasta file

RIC score -> The score with the original RIC

rRIC score -> The score with our retrained RIC

CRoSSeD score -> The score of CRoSSeD

REC score -> The REC score (bare in mind that this is only indicative and not used in the decision of calling the sequence an RSS or not, as both rRIC and CRoSSeD need to pass the thresholds.)

PASS/FAIL RIC -> 1 if the sequence passes the original RIC score thershold

The REC classification (next 2 columns) is defined based on the rRIC score and the CRoSSeD score of each sequence. We have calculated 2 different thresholds. The stringent threshold is the best for genomic searches as it minimizes the number of wrongly identified sequences as RSS. If you desire to look for a specific region and find all sequences that could work as RSS you can look at the sensitive classification (bare in mind that you will have an higher probability to observe false positives in this classification).

PASS/FAIL REC stringent -> 1 if the sequence passes our stringent threshold.

PASS/FAIL REC sensitive -> 1 if the sequence passes our sensitive threshold.



TO DESIGN A NEW MODEL

You can use your own dataset to create a new model following the same model we used. To this purpose, use the script rec_train.py



References:

(1) Cowell LG, Davila M, Kepler TB, Kelsoe G. Identification and utilization of arbitrary correlations in models of recombination signal sequences. Genome Biology. 2002;3(12):research0072.1-research0072.20.

(2) Meysman P, Dang TH, Laukens K, et al. Use of structural DNA properties for the prediction of transcription-factor binding sites in Escherichia coli. Nucleic Acids Research. 2011;39(2):e6. doi:10.1093/nar/gkq1071.

