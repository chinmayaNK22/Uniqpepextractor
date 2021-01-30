import os
import read_fasta_file
import tryptic_peptide
import sys

cwd = os.getcwd()

    
def Unique_pep(infile, outfile, miss_cleave, min_len, max_len):
    tryptic_pep = {}
    list_fasta = {}
    for list_file in os.listdir(infile):
        if os.path.isfile(infile + '/' + list_file):
            if list_file.split('.')[-1] == 'fasta':
                list_fasta[list_file] = list_file
                fasta = read_fasta_file.read_fasta(infile + '/' + list_file)
                for rows in fasta:
                    seq = rows[1].rstrip()
                    for iter_cleavage in range(miss_cleave + 1):
                        pep = tryptic_peptide.tryptic_peptide_trypsin(seq,int(iter_cleavage),int(min_len),int(max_len))
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
                            
    print ('Protein digestion by trypsin is complete and unique peptides are stored')

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
                                write1 = open(path + '/' + list_fasta[j.split('\t')[-1]].rstrip('.fasta') + '_' + 'Unique_Peptides.txt', 'a')
                                write1.write(k + '\t' + j.split('\t')[0].split(' ')[0] + '\t' + j.split('\t')[0] + '\t' + str(len(k)) + '\t' + j.split('\t')[-1] + '\n')

                                write1.close()
                        except:
                            pass
                            #print (list_folder)

        
if __name__== "__main__":
    Unique_pep(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])

#path = 'D:/Skyline/NTMs/ALL_Mycobcaterium_Species/Mycobacteroides_Mycolicibacter/test/ D:/Skyline/NTMs/NTMs_In_silico_Peptides/ 0 7 25'
