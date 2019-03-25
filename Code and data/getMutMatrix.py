# coding: utf-8
"""
This module contains util functions for preprocessing protein sequences 
to get mutation matrix.
"""

from Bio import SeqIO
from Bio.Align.Applications import MafftCommandline
from pathlib import Path
from random import randint

import os
import matlab.engine



def get_seqs_by_protein(subtype_path, prot_name, records, len_threshold):
    """
    This function separates sequences by protein.
    """
    
    target_sequences = []
    
    # allseq_file = subtype_path + 'fullGenomeSet_protein.fasta'
    # records = list(SeqIO.parse(allseq_file, 'fasta'))
    
    for r in records:
        des = r.description
        des_items = des.split('|')
        pnames = [i for i in des_items if i.startswith('Protein Name:')]
        if(len(pnames) == 1):
            pname_record = pnames[0]
            # print(pname_record)
            
            check_flag = 'Protein Name:'+prot_name
            if check_flag in pname_record:
                if len(r.seq) != len_threshold:
                    pass
                else: 
                    target_sequences.append(r)
        else:
            # print(des)
            pass
    
    protein_seqs_dir = subtype_path + '/protein_seqs/'
    if not os.path.exists(protein_seqs_dir):
        os.makedirs(protein_seqs_dir)
            
    SeqIO.write(target_sequences, protein_seqs_dir + prot_name+'.fasta', 'fasta')
    
    print("There are %i %s sequences with length %i" % (len(target_sequences), prot_name, int(len_threshold)))
    
    return True


# EXAMPLE ##################################################    
# h3n2_path = './H3N2/'
# allseq_file = h3n2_path + 'fullGenomeSet_protein.fasta'
# records = list(SeqIO.parse(allseq_file, 'fasta'))
# get_seqs_by_protein(h3n2_path, 'HA', records, 252)
############################################################


def main_get_seqs_by_protein(subtype_path, protein_theta_dict):
    
    # allseq_file = subtype_path + 'fullGenomeSet_complete_protein.fasta'
    
    allseq_file = subtype_path + 'fullGenomeSet_protein.fasta'
    records = list(SeqIO.parse(allseq_file, 'fasta'))
    print('There are %i protein sequences in total.' % len(records))
    
    for pname, theta in protein_theta_dict.items():
        # print(pname, theta)
        get_seqs_by_protein(subtype_path, pname, records, theta)
    
    return True


# Pre-defined constants: The path for the MAFFT program.
mafft_bat = 'D:\mafft-7.409-win64-signed\mafft-win\mafft.bat'

def alignSeqs(in_file, out_file, mafft_bat):
    
    cline = MafftCommandline(mafft_bat, input = in_file)
    # print(cline)
    
    [stdout, stderr] = cline()
    
    with open(out_file, 'w+') as handle:
        handle.write(stdout)
    with open('error.txt', 'w+') as handle:
        handle.write(stderr)
        
    return True




def main_alignSeqs(file_path, protein_list, mafft_bat):
    
    for pname in protein_list:
        print('Aligning %s sequences ...' % pname)
        alignSeqs(file_path + pname + '.fasta', file_path + pname.strip() + '_align.fasta', mafft_bat)
        
    return True





def get_seqs_by_year(subtype_file_path, prot_name):
    
    seq_file = subtype_file_path + 'protein_seqs/' + prot_name.strip() + '_align.fasta' #e.g. HA_align.fasta
    records = list(SeqIO.parse(seq_file, 'fasta'))
    
    for r in records:
        des = r.description
        seq = r.seq
        
        des_items = des.split('|')
        gbid = des_items[0]
        strain_name = ''
        for i in des_items:
            if i.startswith('Strain Name:'):
                strain_name = i
                break
            else:
                pass
        year = strain_name.split('/')[-1]
        
        year_dir = subtype_file_path + 'year_seqs/' + prot_name.strip() + '_year/'
        if not os.path.exists(year_dir):
            os.makedirs(year_dir)
            
        year_file = Path(year_dir + str(year) + '.fasta')
        
        if year_file.exists():
            with open(str(year_file), 'a+') as f:
                f.write('>%s|%s\n%s\n' % (gbid, strain_name, seq))
        else:
            # print(str(year))
            with open(str(year_file), 'w+') as f:
                f.write('>%s|%s\n%s\n' % (gbid, strain_name, seq))
    
    return True
    



def main_get_seqs_by_year(subtype_path, protein_list):
    
    for pname in protein_list:
        print(pname)
        get_seqs_by_year(subtype_path, pname)
        
    return True
    


# ! Remember to check sequences divided by year mannually. 
# Sometimes the strain name does not conform to the convention with YEAR.
# For example, strains beginning with 'H090-' were in 2009; strains beginning with 'SirirajICRC' were in 2010. Merge sequences into the corresponding year_files. 


def getMutMatrix(seq, nYears):
    mutMatrix = []
    for i in range(1, len(seq)):
        # print(i)
        if i%nYears == 0:
            continue
        s1 = seq[i-1]
        s2 = seq[i]
        
        flags = list(i[0]!=i[1] for i in zip(s1, s2))
        
        aRow = []
        for j in range(1, len(flags)+1):
            aRow.append(int(flags[j-1]))
        mutMatrix.append(aRow)
        
    return mutMatrix



def writeMutMatrix(filename, mutMatrix):
    with open(filename, 'w+') as f:
        for i in range(1, len(mutMatrix)+1):
            # print(i)
            Seq_t = mutMatrix[i-1]
            type(Seq_t)
            Seq = ','.join(str(v) for v in Seq_t)
            f.write('%s\n' % Seq)
    return True



def writeTempSeq(filename, seqs, year1, endYear, nloops):
    TempSeqFile = []
    TempSeq = []
    # print(seqs.keys())
    
    startYear = year1
    for y in range(year1, endYear+1):
        if y in seqs.keys():
            continue
        else:
            startYear = y+1
            
    for i in range(1, nloops+1):
        # print(i)
        for y in range(startYear, endYear+1):
            if len(seqs[y])>1:
                N = len(seqs[y])
                r = randint(1, N)
                TempSeq.append(seqs[y][r-1])
            elif len(seqs[y])==1:
                TempSeq.append(seqs[y][0])
    TempSeqFile.append(TempSeq)
    TempSeqFile.append(startYear)
    
    with open(filename, 'w+') as f:
        for i in range(1, len(TempSeq)+1):
            s = TempSeq[i-1]
            f.write('%s\n' % s)
    # print("WriteTempSeq over")
    
    return TempSeqFile
                    



def writeSeqs(filename, seqs, keyl, keyu):
    with open(filename, 'w+') as f:
        for key, value in seqs.items():
            seqNum = len(value)
            for i in range(1, seqNum+1):
                seq = value[i-1]
                s = str(key) + '\t' + seq
                f.write('%s\n' % s)
    return True


def main_getMutMatrix(subtype_path, pname, year1, year2, nloops):
    filepath = subtype_path + 'year_seqs/' + pname.strip() + '_year/'
    
    AllSeq = {}
    for y in range(year1, year2+1):
        filename = filepath+str(y)+'.fasta' 
        Seqs = []
        
        records = list(SeqIO.parse(filename, 'fasta'))
        for r in records:
            Seqs.append(r.seq)
        AllSeq[y] = Seqs
    
    writeSeqs(filepath + pname.strip() + '_allseq.fasta', AllSeq, year1, year2)
    TempSeqFile = writeTempSeq(filepath + pname.strip() + '_tempSeq.fasta', AllSeq, year1, year2, nloops)
    TempSeq = TempSeqFile[0]
    StartYear = TempSeqFile[1]
    
    nY = year2-StartYear + 1
    
    matrix = getMutMatrix(TempSeq, nY)
    
    matrix_dir = subtype_path + 'mutMatrix/'
    if not os.path.exists(matrix_dir):
        os.makedirs(matrix_dir)
    
    writeMutMatrix(matrix_dir + pname.strip() + '_mutMatrix.txt', matrix)
    
    return True



def call_mutMatrix2ARM(subtype_path, protein_names_list):
    eng = matlab.engine.start_matlab()
    mutMatrix_path = subtype_path + 'mutMatrix/'
    
    for pname in protein_names_list:
        print(pname)
        eng.mutMatrix2ARM(mutMatrix_path, pname, nargout = 0)
        
    return True



def remove_invalid_lines(input_csv, output_csv):
    # remove INVALID lines
    # including blank lines and lines with only one site. 
    
    output = open(output_csv, 'w+')
    with open(input_csv, 'r') as f:
        for line in f:
            if line.strip():
                if len(line.split(' '))>1:
                    # print(line.split(' '))
                    output.write(line)
            else:
                # print(line)
                pass
        output.close()
        
    return True



def main_remove_invalid_lines(file_path, protein_names_list):
    print('--------- %s ------------' % file_path)
    for pname in protein_names_list:
        # print(pname)
        input_csv = file_path + pname.strip() + '_forARM.csv'
        forARM_dir = file_path + 'forARM_cleaned/'
        if not os.path.exists(forARM_dir):
            os.makedirs(forARM_dir)
            
        output_csv = forARM_dir + pname.strip() + '_forARM.dat'
        remove_invalid_lines(input_csv, output_csv)
    return True

