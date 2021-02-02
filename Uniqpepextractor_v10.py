import os
import os
import read_fasta_file
import tryptic_peptide
import sys
from datetime import date, datetime

cwd = os.getcwd()

today = date.today()
dt = today.strftime("%m%d%y")

now = datetime.now()
current_time = now.strftime("%H%M%S")

aa_known = {}
aa_unknown = {}
def amino_acids(AA_seq):
    aa_dicts = {'G': 'Glycine@Gly', 'A': 'Alanine@Ala', 'L': 'Leucine@Leu', 'M': 'Methionine@Met', 'F': 'Phenylalanine@Phe', 'W': 'Tryptophan@Trp', 'K': 'Lysine@Lys', 'Q': 'Glutamine@Gln',
            'E': 'Glutamic Acid@Glu', 'S': 'Serine@Ser', 'P': 'Proline@Pro', 'V': 'Valine@Val', 'I': 'Isoleucine@Ile', 'C': 'Cysteine@Cys', 'Y': 'Tyrosine@Tyr', 'H': 'Histidine@His',
            'R': 'Arginine@Arg', 'N': 'Asparagine@Asn', 'D': 'Aspartic Acid@Asp', 'T': 'Threonine@Thr', 'U' : 'Selenocysteine@Sec'}
    for a in AA_seq.upper():
        if a in aa_dicts:
            aa_known[a] = aa_dicts[a]
        else:
            aa_unknown[a] = ''                       

    for l, m in aa_unknown.items():
        #print (l)
        return l


   
def Unique_pep(infile, outfile, miss_cleave, min_len, max_len):
    tryptic_pep = {}
    list_fasta = {}
    fasta_file = [infile + '/' + list_file for list_file in os.listdir(infile) if os.path.isfile(infile + '/' + list_file) if list_file.split('.')[-1] == 'fasta']
    for f in fasta_file:
        list_fasta[f.split('/')[-1]] = f.split('/')[-1]
        fasta = read_fasta_file.read_fasta(f)
        for rows in fasta:
            seq = rows[1].rstrip()
            aa_count = amino_acids(seq)
            for iter_cleavage in range(int(miss_cleave) + 1):
                pep = tryptic_peptide.tryptic_peptide_trypsin(seq,iter_cleavage,int(min_len),int(max_len))
                for i in pep:
                #print (i, rows[0], protein_fasta_file)
                    if 'C' in i:
                        rm = i
                    elif 'M' in i:
                        rm = i
                    elif 'X' in i:
                        rm = i
                    elif 'Z' in i:
                        rm = i
                    else:
                        if i not in tryptic_pep:
                            tryptic_pep[i] = [rows[0] + '\t' + f.split('/')[-1]]

                        else:
                            tryptic_pep[i].append(rows[0] + '\t' + f.split('/')[-1])
                            
    print ('Protein digestion by trypsin is complete and ' + str(len(tryptic_pep)) + ' unique peptides are stored')
    print ('Number of known amino acids found: ' + str(len(aa_known)))
    print ('Number of unknown amino acids found: ' + str(len(aa_unknown)))

    for k, v in aa_known.items():
        print (k + '\t' + v.split('@')[0] + '\t' + v.split('@')[1])

    for l, m in aa_unknown.items():
        print (l)

    for k, v in tryptic_pep.items():
        if len(v) == 1:
            for j in v:
                if j.split('\t')[-1] in list_fasta:
                    for list_folder in os.listdir(outfile):
                        file = list_folder + '.fasta'
                        try:
                            if j.split('\t')[-1] in file:
                                path = outfile + '/' + list_folder
                                #print (path + '/' + list_fasta[j.split('\t')[-1]].rstrip('.fasta') + '_' + 'Unique_Peptides.txt')
                                #print (k + '\t' + j.split('\t')[0].split(' ')[0] + '\t' + j.split('\t')[0] + '\t' + str(len(k)) + '\t' + j.split('\t')[-1] + '\n')
                                write1 = open(path + '/' + list_fasta[j.split('\t')[-1]].rstrip('.fasta') + '_Unique_Peptides_' + dt + '.txt', 'a')
                                write1.write(k + '\t' + j.split('\t')[0].split(' ')[0] + '\t' + j.split('\t')[0] + '\t' + str(len(k)) + '\t' + j.split('\t')[-1] + '\n')

                                write1.close()
                        except:
                            pass
                            #print (list_folder)

        
if __name__== "__main__":
    Unique_pep(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])

#path = 'python Uniqpepextractor_v10.py D:/Skyline/NTMs/ALL_Mycobcaterium_Species/ D:/Skyline/NTMs/NTMs_In_silico_Peptides/ 0 7 25'
