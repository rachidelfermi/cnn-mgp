# Load packages
import os, argparse
from tqdm import tqdm
import numpy as np
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import tensorflow as tf
tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)
from tensorflow.keras.models import load_model

def validate_args(args):
    # Validate input file location
    if not os.path.isfile(args.fileName):
        print('The specified input file does not exist (i.e., "' + args.fileName + ')".')
        print('Make sure you entered this correctly and try again.')
        quit()
    # Validate parameters
    if args.minOrfLenght < 80:
        print("minProLen must be >= 80. Fix your input and try again.")
        quit()
    if args.unresolvedCodon != 0 and args.unresolvedCodon != 1:
        print("unresolvedCodon must be 1 for (ATG, CTG, GTG, TTG) or 0 for (ATG). Fix your input and try again.")
        quit()
    # Validate output file location
    if args.outputFileName != None:
        print(os.path.isfile(args.outputFileName),args.outputFileName)
        if os.path.isfile(args.outputFileName):
            print('There is already a file at "' + args.outputFileName + '".')
            print('Either specify a new file name, delete this older file')
            quit()
        elif args.sequenceType.lower() == 'nucl':
            args.nuclOutName = None
            args.coordOutName= None
            prefixNucName = args.outputFileName.rsplit('.', maxsplit=1)
            if len(prefixNucName) == 1:
                prefixNucName.append('.fasta')
            args.nuclOutName = prefixNucName[0] + '_nucl.' + prefixNucName[1]
            args.coordOutName = prefixNucName[0] + '_CNN-MGP_coord.coords'
            if os.path.isfile(args.nuclOutName) :
                print(
                    'There is already a file "' + args.nuclOutName + '". Either specify a new file name, delete this file')
                quit()
            if os.path.isfile(args.coordOutName):
                print(
                    'There is already a file  "' + coordOutName + '". Either specify a new file name or delete this older file')
                quit()

        elif args.sequenceType.lower() == 'prot':
            args.protOutName = None
            args.coordOutName = None
            prefixprotName = args.outputFileName.rsplit('.', maxsplit=1)
            if len(prefixprotName) == 1:
                prefixprotName.append('.fasta')
            args.protOutName = prefixprotName[0] + '_prot.' + prefixprotName[1]
            args.coordOutName = prefixprotName[0] + '_CNN-MGP_coord.coords'
            if os.path.isfile(args.protOutName) :
                print(
                    'There is already a file  "' + args.protOutName + '". Either specify a new file name, delete this file')
                quit()
            if os.path.isfile(args.coordOutName):
                print(
                        'There is already a file  "' + coordOutName + '". Either specify a new file name or delete this older file')
                quit()

        else:
            args.protOutName = None  # For consistency, since the below else statement generates these values
            args.nuclOutName = None
            args.coordOutName = None
            # Derive output file names
            outPrefix = args.outputFileName.rsplit('.', maxsplit=1)
            if len(outPrefix) == 1:
                outPrefix.append('.fasta')
            args.protOutName = outPrefix[0] + '_prot.' + outPrefix[1]
            args.nuclOutName = outPrefix[0] + '_nucl.' + outPrefix[1]
            args.coordOutName = outPrefix[0] + '_CNN-MGP_coord.coords'
            if os.path.isfile(args.protOutName) :
                print(
                    'There is already a file  "' + args.protOutName + '". Either specify a new file name or delete this file')
                quit()
            if os.path.isfile(args.nuclOutName) :
                print(
                    'There is already a file  "' + args.nuclOutName + '". Either specify a new file name or delete this older file')
                quit()
            if os.path.isfile(args.coordOutName):
                print(
                    'There is already a file  "' + coordOutName + '". Either specify a new file name or delete this older file')
                quit()

    return args
# Reqs
usage = """%(prog)s predict gene in Metagenomics fragments """
p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", "-input", dest="fileName",
                   help="Input fasta file name")
p.add_argument("-o", "-output", dest="outputFileName",
                   help="Output fasta file name")
# Opts
p.add_argument("-min", "--minimum", type=int, dest="minOrfLenght",
                   help="Minimum ORF length.any number lower than 80 well be override to 80, Default = 80.", default=80)

p.add_argument("-st", "--seqtype", dest="sequenceType", choices = ['prot', 'nucl', 'both', 'PROT', 'NUCL', 'BOTH'],
                   help="Specify the type of output you want to generate (i.e., protein translated gene, nucleotide CDS, or both). If you specify 'both', two outputs with '_prot' and '_nucl' suffix will be generated. Default == 'nucl'.", default="nucl")

p.add_argument("-u", "--unresolved", dest="unresolvedCodon", type=int,
                   help="choose 1 for  unresolved codons(ie, ATG, CTG, GTG, TTG).0 for codon start ATG.", default=1)
args = p.parse_args()
if args.fileName != None and args.outputFileName != None:
    args = validate_args(args)
else:
    print(
        "you didn't specify either the input file or the output file , shifting to interactive mode."
    )
fileName = args.fileName
outputFileName = args.outputFileName
minOrfLenght = args.minOrfLenght
sequenceType = args.sequenceType
unresolvedCodon = args.unresolvedCodon
coordOutName = None

if fileName == None or outputFileName == None:
    # Locate our file of interest
    print(
        'Enter the name of the fasta-formatted file you wish to extract ORFs from. This should be a nucleotide sequence file. Include file extension (e.g., ".fasta").')
    while True:
        try:
            fileName = input()
            print(fileName)
            print(os.path.isfile(fileName))
            if os.path.isfile(fileName) == False:
                raise Exception
            print('Fasta file located successfully')
            break
        except KeyboardInterrupt:
            quit()
        except:
            print(
                'Fasta file failed to load. If you misspelt the name, try again. If this script isn\'t in the same directory as the fasta file, move it there then try again.')
            continue
    print('')


    print(
        'Enter the name which you want the output fasta file to be called. Make sure not to use illegal characters (i.e. \\/:?"<>|).')
    while True:
        try:
            illegalCharacters = '\/:?"<>|;,'
            outputFileName = input()
            if outputFileName == '':
                print(
                    'You didn\'t name this file anything. You need to have at least one character in your output file name. Try again.')
                continue
            for character in illegalCharacters:
                if character in outputFileName:
                    raise Exception
            if os.path.isfile(outputFileName):
                print('This is already a file. Try again.')
                continue
            break
        except KeyboardInterrupt:
            quit()
        except:
            print(
                'You used an illegal character (i.e. \\/:?"<>|). Try to name your file without these characters again.')
            continue

    # Allow user to determine minimum protein length
    print(
        'Enter the minimum ORF length that you want to accept as legitimate ORF. A setting of 80 is recommended. 80 is taken if the user types any value < 80 ')
    while True:
        minOrfLenght = input()
        try:
            minOrfLenght = int(minOrfLenght)
            minOrfLenght=80 if minOrfLenght< 80 else minOrfLenght
            break
        except KeyboardInterrupt:
            quit()
        except:
            print('You seem to have not typed a number here. Try again.')
    print('')

    print(
        'Enter 1 for alternative start codons (ATG, TTG ,CTG, GTG) or 0 for only (ATG) start codon.')
    print('')
    while True:
        unresolvedCodon = input()
        try:
            if int(unresolvedCodon)==0 or int(unresolvedCodon)==1:
                unresolvedCodon=int(unresolvedCodon)
                break
            elif int(unresolvedCodon)!=0 and int(unresolvedCodon)!=1:
                print('choose 1 or 0. Try again.')
                continue
        except ValueError:
            print("The input was not a valid integer.")
        except:
            print('You seem to have not typed a number here. Try again.')

    # Allow user to determine what kind of output they want to have
    print(
        'Do you want to produce a translated protein ORF file, a nucleotide CDS file, or both? Enter \'prot\', \'nucl\', or \'both\'')
    print(
          'Note that \'prot\' and  \'both\' affect the speed of the script')
    print('')
    seqtypeChoices = ['prot', 'nucl', 'both', 'PROT', 'NUCL', 'BOTH']
    while True:
        try:
            sequenceType = input()
            if sequenceType.lower() not in seqtypeChoices:
                raise Exception
            break
        except KeyboardInterrupt:
            quit()
        except:
            print('You didn\'t type \'prot\', \'nucl\', or \'both\'. Try again.')
    print('')

# Check if we should be overwriting files / get our output names if sequenceType.lower() == 'both'

if outputFileName != None:
    if os.path.isfile(outputFileName):
        print('There is already a file at "' + outputFileName + '".')
        print('Either specify a new file name, delete this older file')
        quit()
    elif sequenceType.lower() == 'nucl':
        nuclOutName = None
        prefixNucName = outputFileName.rsplit('.', maxsplit=1)
        if len(prefixNucName) == 1:
            prefixNucName.append('.fasta')
        nuclOutName = prefixNucName[0] + '_nucl.' + prefixNucName[1]
        coordOutName = prefixNucName[0] + '_CNN-MGP_coord.coords'
        if os.path.isfile(nuclOutName):
            print(
                'There is already a file "' + nuclOutName + '". Either specify a new file name, delete this file')
            quit()
        if os.path.isfile(coordOutName):
            print(
                'There is already a file  "' + coordOutName + '". Either specify a new file name or delete this older file')
            quit()

    elif sequenceType.lower() == 'prot':
        protOutName = None
        prefixprotName = outputFileName.rsplit('.', maxsplit=1)
        if len(prefixprotName) == 1:
            prefixprotName.append('.fasta')
        protOutName = prefixprotName[0] + '_prot.' + prefixprotName[1]
        coordOutName = prefixprotName[0] + '_CNN-MGP_coord.coords'
        if os.path.isfile(protOutName):
            print(
                'There is already a file  "' + protOutName + '". Either specify a new file name, delete this file')
            quit()
        if os.path.isfile(coordOutName):
            print(
                'There is already a file  "' + coordOutName + '". Either specify a new file name or delete this older file')
            quit()

    else:
        protOutName = None  # For consistency, since the below else statement generates these values
        nuclOutName = None
        # Derive output file names
        outPrefix = outputFileName.rsplit('.', maxsplit=1)
        if len(outPrefix) == 1:
            outPrefix.append('.fasta')
        protOutName = outPrefix[0] + '_prot.' + outPrefix[1]
        nuclOutName = outPrefix[0] + '_nucl.' + outPrefix[1]
        coordOutName = outPrefix[0] + '_CNN-MGP_coord.coords'
        if os.path.isfile(protOutName):
            print(
                'There is already a file  "' + protOutName + '". Either specify a new file name or delete this file')
            quit()
        if os.path.isfile(nuclOutName):
            print(
                'There is already a file  "' + nuclOutName + '". Either specify a new file name or delete this older file')
            quit()
        if os.path.isfile(coordOutName):
            print(
                'There is already a file  "' + coordOutName + '". Either specify a new file name or delete this older file')
            quit()


# loading the  models
projectpath = str(os.getcwd())
import platform
separator= "\\" if platform.system()=='Windows' else '/'
model_Path = projectpath + separator+ "models"
# reading in splice junction input data and converting to required format
CGrange0 = load_model(model_Path + separator+'cgrange0fulldata.h5')
CGrange1 = load_model(model_Path + separator+'cgrange1fulldata.h5')
CGrange2 = load_model(model_Path + separator+'cgrange2fulldata.h5')
CGrange3 = load_model(model_Path + separator+'cgrange3fulldata.h5')
CGrange4 = load_model(model_Path + separator+'cgrange4fulldata.h5')
CGrange5 = load_model(model_Path + separator+'cgrange5fulldata.h5')
CGrange6 = load_model(model_Path + separator+ 'cgrange6fulldata.h5')
CGrange7 = load_model(model_Path + separator+'cgrange7fulldata.h5')
CGrange8 = load_model(model_Path + separator+'cgrange8fulldata.h5')
CGrange9 = load_model(model_Path + separator+'cgrange9fulldata.h5')
#all the fuction needed for thr preproccessing
def get_candidategenes(orfinput, cgrange):
    if cgrange <= 0.3657:
        class_pred_full = CGrange0.predict_classes(orfinput)
    if cgrange > 0.3657 and cgrange <= 0.4157:
        class_pred_full = CGrange1.predict_classes(orfinput)
    if cgrange > 0.4157 and cgrange <= 0.46:
        class_pred_full = CGrange2.predict_classes(orfinput)
    if cgrange > 0.46 and cgrange <= 0.5014:
        class_pred_full = CGrange3.predict_classes(orfinput)
    if cgrange > 0.5014 and cgrange <= 0.5428:
        class_pred_full = CGrange4.predict_classes(orfinput)
    if cgrange > 0.5428 and cgrange <= 0.5814:
        class_pred_full = CGrange5.predict_classes(orfinput)
    if cgrange > 0.5814 and cgrange <= 0.6185:
        class_pred_full = CGrange6.predict_classes(orfinput)
    if cgrange > 0.6185 and cgrange <= 0.65:
        class_pred_full = CGrange7.predict_classes(orfinput)
    if cgrange > 0.65 and cgrange <= 0.6828:
        class_pred_full = CGrange8.predict_classes(orfinput)
    if cgrange > 0.6828 and cgrange <= 1:
        class_pred_full = CGrange9.predict_classes(orfinput)
    return class_pred_full

#onehot encoding
def onehote(sequence):
    mapping = {"A": 0, "C": 1, "G": 2, "T": 3}
    seq2 = [mapping[i] for i in sequence]
    return np.eye(4)[seq2]

# padding
def padding(onehotenc):
    # make the lenght 700 so the padding becomes 700
    addition = [[0, 0, 0, 0] for i in range(0, 700 - len(onehotenc[0]), 1)]
    onehotenc[0] = np.concatenate([onehotenc[0], addition])
    return np.stack(tf.keras.preprocessing.sequence.pad_sequences(onehotenc, padding="post"))

#checking the dna
def check_Dna(orf):
    index = 0
    dna_checked = False
    codonA, codonT, codonC, codonG = 0, 0, 0, 0
    if len(orf) > minOrfLenght :
        for ind, codon in enumerate(orf):
            codonA += codon.count('A')
            codonT += codon.count('T')
            codonC += codon.count('C')
            codonG += codon.count('G')
            if codon.count('A') + codon.count('C') + codon.count('T') + codon.count('G') == 3:
                dna_checked = True
                index = ind + 1
            elif len(orf[:ind]) > minOrfLenght:
                dna_checked = True
                index = ind
                break
            else:
                break
    else:
        dna_checked = False
    if codonA > 1 and codonT > 1 and codonC > 1 and codonG > 1:
        return dna_checked, orf[:index]
    else:
        return False, []

def queryPosition(strt,lenth):
   return strt*3,strt*3+lenth*3
def strand(x):
    return "+" if x < 3 else "-"
translateaminoacides = {    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
                            'TGC':'C', 'TGT':'C', 'GAC':'D', 'GAT':'D',
                            'GAA':'E', 'GAG':'E', 'TTC':'F', 'TTT':'F',
                            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
                            'CAC':'H', 'CAT':'H', 'ATA':'I', 'ATC':'I',
                            'ATT':'I', 'AAA':'K', 'AAG':'K', 'CTT':'L',
                            'TTA':'L', 'TTG':'L',
                            'CTA':'L', 'CTC':'L', 'CTG':'L', 'ATG':'M',
                            'AAC':'N', 'AAT':'N', 'CCA':'P', 'CCC':'P',
                            'CCG':'P', 'CCT':'P', 'CAA':'Q', 'CAG':'Q',
                            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
                            'AGA':'R', 'AGG':'R', 'AGC':'S', 'AGT':'S',
                            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
                            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
                            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
                            'TGG':'W', 'TAC':'Y', 'TAT':'Y', 'TAA':'*',
                            'TAG':'*','TGA':'*'}

def create_OutputFiles(predictedGene):
    if sequenceType.lower() == 'nucl':
        nucORFsfile = open(nuclOutName, 'a+')
        nucORFsfile.write("".join(predictedGene) + "\n")
    elif sequenceType.lower() == 'prot':
        protORFsfile = open(protOutName, 'a+')
        for codon in range(0,len(predictedGene),3):
            protORFsfile.write(translateaminoacides[predictedGene[codon:codon+3]])
        protORFsfile.write("\n")
    else:
        protORFsfile = open(protOutName, 'a+')
        nucORFsfile = open(nuclOutName, 'a+')
        for codon in range(0, len(predictedGene), 3):
            protORFsfile.write(translateaminoacides[predictedGene[codon:codon + 3]])
        protORFsfile.write("\n")
        nucORFsfile.write("".join(predictedGene) + "\n")
startCodons=[]
stopCodons = ["TAA", "TAG", "TGA"]
if unresolvedCodon == 1:
    startCodons = ["ATG", "TTG", "GTG", "CTG"]
else:
    startCodons = ["ATG"]
orfscoords = open(coordOutName, 'a+')
# extracting the sequences from the multifasta file
sequences = []
descr = None
# here is the path of multifalsta file
with open(fileName) as file:
    line = file.readline()[:-1]  # always trim newline
    while line:
        if line[0] == '>':
            if descr:  # any sequence found yet?
                sequences.append((descr, seq))
            descr = str(line[1:].split('>'))
            seq = ''  # start a new sequence
        else:
            seq += line
        line = file.readline()[:-1]
    sequences.append((descr, seq))
print(len(sequences), " Fragment found !")

print('Scanning Fragments for coding ORFS')

for index, value in tqdm(
        enumerate(sequences)):  # extracting all orf in each fragment then predicting  coding from non-coding orfs
    Frames = []
    dna = value[1]  # extract the fragment
    description = value[0]
    Frames = []
    reversedna = []
    listOfOrfs = []
    listOfOrfID = []
    # create the positive frames
    # split the frames into codons for better performance
    Frames.append([dna[i:i + 3] for i in range(0, len(dna), 3)])
    Frames.append([dna[i:i + 3] for i in range(1, len(dna), 3)])
    Frames.append([dna[i:i + 3] for i in range(2, len(dna), 3)])
    # reverse compliment of the fragment
    reverse = {"A": "T", "C": "G", "T": "A", "G": "C"}
    for i in range(len(dna)):
        reversedna.append(reverse[dna[-i - 1]]) if dna[-i - 1] in reverse.keys() else reversedna.append(
            dna[-i - 1])  # if any contamination found we keep it for further more check
    reversedna = ''.join(reversedna)
    # create the negative frames
    Frames.append([reversedna[i:i + 3] for i in range(0, len(reversedna), 3)])
    Frames.append([reversedna[i:i + 3] for i in range(1, len(reversedna), 3)])
    Frames.append([reversedna[i:i + 3] for i in range(2, len(reversedna), 3)])
    #  our models predict genes for each cgrange with a different model
    CG_count = (dna.count('C') + dna.count('G')) / len(dna)


    def extract_IncompleteORFS_From_TheStart(stran, frame):
        # get the first tree indixes of stop codons
        lenth = len(
            frame)  # if no stop or start codon found we take the lenght of the frame(only start orfs, only stop orfs or orfs missing both start and stop)
        TAG = frame.index("TAG") if "TAG" in frame else lenth
        TAA = frame.index("TAA") if "TAA" in frame else lenth
        TGA = frame.index("TGA") if "TGA" in frame else lenth
        # get the first four indixes of start codons
        ATG = frame.index("ATG") if "ATG" in frame else lenth
        TTG = frame.index("TTG") if "TTG" in frame else lenth
        GTG = frame.index("GTG") if "GTG" in frame else lenth
        CTG = frame.index("CTG") if "CTG" in frame else lenth
        if unresolvedCodon == 1:
            minStart = min(ATG, TTG, GTG, CTG)
        else:
            minStart = ATG
        if min(TAG, TAA, TGA) < minStart:  # the first codon found is stop : ------------stop--
            goodDna, orf = check_Dna(frame[:min(TAG, TAA, TGA) + 1])  # check the orf lenght and contamination
            if goodDna:
                listOfOrfs.append(''.join(orf))
                qstart, qend = queryPosition(0, min(TAG, TAA, TGA))
                listOfOrfID.append([index * 2, strand(stran), qstart, qend, stran, "I"])
        elif min(TAG, TAA, TGA) == lenth and minStart < lenth:  # stop codon not found  but the start codons is, we take the first start codon :---start-----------------------
            goodDna, orf = check_Dna(frame[minStart:])  # check the orf lenght and contamination
            if goodDna:
                listOfOrfs.append(''.join(orf))
                qstart, qend = queryPosition(minStart, len(orf))
                listOfOrfID.append([index * 2, strand(stran), 1, qstart, qend, stran, "I"])
        elif min(TAG, TAA, TGA) == lenth and minStart== lenth:  # neather stop or start codon found: ------------------------
            goodDna, orf = check_Dna(check_Dna(frame[:]))
            if goodDna:
                listOfOrfs.append(''.join(orf))
                qstart, qend = queryPosition(0, len(orf))
                listOfOrfID.append([index * 2, strand(stran), 1, qstart, qend, stran, "I"])


    def extract_IncompleteORFS_From_TheEnd(stran, frame):
        lenth = len(frame)
        # searching for the last stop codon
        TAG = list(reversed(frame)).index("TAG") if "TAG" in frame else 0
        TAA = list(reversed(frame)).index("TAA") if "TAA" in frame else 0
        TGA = list(reversed(frame)).index("TGA") if "TGA" in frame else 0
        # getting the last stop codon index
        if min(TAG, TAA,TGA) != 0:  # if no codon stop found means eather only start or both codon missing, we already extacted that
            lastStopCodon = lenth - min(TAG, TAA, TGA) - 1  # getting the last index if stop codons
            # searching for anypossible  start codont after the last stop codon: --------------stop-----start-----------
            ATG = frame[lastStopCodon:].index("ATG") if "ATG" in frame[lastStopCodon:] else lenth
            TTG = frame[lastStopCodon:].index("TTG") if "TTG" in frame[lastStopCodon:] else lenth
            GTG = frame[lastStopCodon:].index("GTG") if "GTG" in frame[lastStopCodon:] else lenth
            CTG = frame[lastStopCodon:].index("CTG") if "CTG" in frame[lastStopCodon:] else lenth
            if unresolvedCodon == 1:
                minStart = min(ATG, TTG, GTG, CTG)
            else:
                minStart = ATG
            if minStart < len(frame):  # if no stop codon found  we exit the function otherwise we check the orf
                newStart = lastStopCodon + minStart  # getting the new start position
                goodDna, orf = check_Dna(frame[newStart:])  # checking the orf for lenght and contamination
                if goodDna:
                    listOfOrfs.append(''.join(orf))
                    qstart, qend = queryPosition(newStart, len(orf))
                    listOfOrfID.append([index * 2, strand(stran), 1, qstart, qend, stran, "I"])


    # the extraction of the all orfs
    for i in range(0, len(Frames), 1):
        # delete the contamination if first codon or the last codon to boost the performance for the checkDna function
        if len(Frames[i][-1]) < 3 or "N" in Frames[i][-1]:
            del Frames[i][-1]
        if "N" in Frames[i][0]:
            del Frames[i][0]
        # extract all the incomplete orfs
        extract_IncompleteORFS_From_TheEnd(i, Frames[i])
        extract_IncompleteORFS_From_TheStart(i, Frames[i])
        # extract  all the complete  orfs
        start = 0
        while start < len(Frames[i]):
            # get the index of the first start
            if Frames[i][start] in startCodons:
                # get the index of the first first stop
                for stop in range(start + 1, len(Frames[i]), 1):
                    if Frames[i][stop] in stopCodons:
                        goodDna, orf = check_Dna(Frames[i][start:stop + 1])  # check lenght of the dna, contamination
                        if goodDna:
                            listOfOrfs.append(''.join(orf))
                            qstart, qend = queryPosition(start, len(orf))
                            listOfOrfID.append([index * 2, strand(i), 1, qstart, qend, i, "C"])
                            start = stop  # skip into the index after stop  then repeat the search
                        break
            start += 1
    # after the extraction of all orfs for this fragment; it's time for preprocessing the orfs,
    # the models only accepts orfs with lenght 700*4
    # first  all orfs should be one hot encoded  the (lenght*4)
    # second we performe padding (adding vectors of zeros (0,0,0,0) until the lenght 700 is satisfied
    if len(listOfOrfs) != 0:  # did we extract any ORFS
        onehotencodedorf = []
        # performe onehot encoding
        for seq in listOfOrfs:
            onehotencodedorf.append(onehote(seq))
        # padd the lenght to 700
        inputOrfs = padding(onehotencodedorf)
        # each cgrange has a specific model
        # feed the orfs into the CNN-MGP alongside the CG content of the fragment
        candidategenes = get_candidategenes(inputOrfs, CG_count)
        # add the predicted  genes into the outputfile
        for i in range(0, len(candidategenes), 1):
            if candidategenes[i] == 1:
                create_OutputFiles(listOfOrfs[i])
                orfscoords.write(">" + str(listOfOrfID[i]).replace('[', '').replace(']', '') + "\n")

orfscoords.close()
print('Program completed successfully!')
