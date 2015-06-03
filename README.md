#########
#PURPOSE#
#########
REC allows you to test a set of sequences and know if they can work as recombination signal sequences (RSS) recognized by RAG.
You can also develop a new classifier with your on training set.

#############
#INSTALATION#
#############
Just copy the folder for you computer and it's ready to use.
##############
#REQUIREMENTS#
##############
python-2.7
biopython-1.57
REC uses 2 other tools in order to classity sequences as RSS or not.
The specific requriments for these tools should also be meet.
RIC - perl-5.14
CRoSSeD - Only tested under Linux system.
        - Perl-5.8
        - CRF++ toolkit ( http://crfpp.sourceforge.net/)
#######
#USAGE#
#######
To test the a set of sequences you should use the script rec_test.py
Options:
-h --help -> show the help message for this script
-t --test -> Input file to be tested. All sequences should have length equal to either 28 or 39.
-o --output -> output file to write the result. If it exist it will be replaced.
-m --model -> If you have train REC with another model give the path to that model where. Otherwise our default trainig set and model will be used.
--algorithm -> You can select wich algorithm will be used to test you sequences. By default the REC option is selected and in the output all algorithms will be represented. If you select any of the others only the selected algorithm will be in the output.

OUTPUT
The default output or with REC option:
Sequence ID -> The id you provid for each sequence
RIC score -> The score with the original RIC
rRIC score -> The score with our retrain RIC
CRoSSeD score -> The score of CRoSSeD
REC score -> The REC score (bare in mind that this is only indicative and not used in the decision of calling the sequence an RSS or not.)
PASS/FAIL RIC -> 1 if the sequence pass the original RIC score thershold 
The REC classification (next 2 columns) is defined based on the rRIC score and the CRoSSeD score of each sequences. We have calculated 2 different thresholds. The stringent threshold is the best for genomic searches as it minimizes the number of wrongly indentified sequences as RSS. If you desire to look for a specific region and find all sequences that could work as RSS you can look at the sensitive classification (bare in mind that you will have an higher probability of observe false positives in this classification).
PASS/FAIL REC stringent -> 1 if the sequence pass our stringent threshold.
PASS/FAIL REC sensitive -> 1 if the sequence pass our sensitive threshold.



TO DESIGN A NEW MODEL
You can use your own dataset to create a new model following the same model we used.
Use the script rec_train.py.
