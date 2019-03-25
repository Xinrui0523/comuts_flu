# coding: utf-8


import subprocess
import os
import csv
import pandas as pd



spmf_jar = './jar_SPMF/spmf.jar'

h3n2_path = './H3N2/'
ph1n1_path = './H1N1/pH1N1/'
eph1n1_path = './H1N1/epH1N1/'
b_path = './B/_b_all_lineages/'

# All proteins 
h3n2_pnames_list_ori = ['HA', 'M1', 'M2', 'NA', 'NP', 'NS1', 'NS2', 'PA', 'PA-N155', 'PA-N182', 'PB1', 'PB1-F2', 'PB1-N40', 'PB2']
ph1n1_pnames_list_ori = ['HA', 'M1', 'M2', 'NA', 'NP', 'NS1', 'NS2', 'PA', 'PA-N155', 'PA-N182', 'PA-X', 'PB1', 'PB1-N40', 'PB2']
eph1n1_pnames_list_ori = ['HA', 'M1', 'M2', 'NA', 'NP', 'NS1', 'NS2', 'PA', 'PA-N155', 'PA-N182', 'PA-X', 'PB1', 'PB1-F2', 'PB1-N40', 'PB2']
b_pnames_list_ori = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'NB', 'M1', 'BM2', 'NS1', 'NS2']


# Predefined proteins to analyze for each subtype
h3n2_pnames_list = ['HA_NA']
ph1n1_pnames_list = ['HA_NA']
eph1n1_pnames_list = ['HA_NA']
b_pnames_list = ['HA_NA']
####################################################


def convert2SPMF(in_file, out_file):
    
    tmp_writer = open(in_file + '.tmp', 'w+')
    
    with open(in_file, 'r') as f:
        lines = f.readlines()
        nlines = len(lines)
        # print(nlines)
        for line in lines:
            newline = ' '.join(line.split()).replace(' ', ',')
            tmp_writer.write(newline + '\n')
                    
    tmp_writer.close()
    
    subprocess.call(['java', '-jar', spmf_jar, 'run', 'Convert_a_sequence_database_to_SPMF_format', in_file + '.tmp', out_file, 'CSV_INTEGER', str(nlines)])
    
    os.remove(in_file + '.tmp')
    
    return True



def main_convert(subtype):
    subtype_path = ''
    pnames_list = []
    
    if subtype == 'h3n2':
        subtype_path = h3n2_path
        pnames_list = h3n2_pnames_list
        
    elif subtype == 'ph1n1':
        subtype_path = ph1n1_path
        pnames_list = ph1n1_pnames_list
    
    elif subtype == 'eph1n1':
        subtype_path = eph1n1_path
        pnames_list = eph1n1_pnames_list
    
    elif subtype == 'b':
        subtype_path = b_path
        pnames_list = b_pnames_list
    
    outfile_spmf_dir = subtype_path + 'spmf_file/'
    if not os.path.exists(outfile_spmf_dir):
        os.makedirs(outfile_spmf_dir)
            
    for p in pnames_list:
        in_file_path = subtype_path + 'mutMatrix/forARM_cleaned/'
        in_file = in_file_path + p + '_forARM.dat'
        out_file = outfile_spmf_dir + p + '_spmf.dat'
        convert2SPMF(in_file, out_file)
        
    return True



def ERMiner(in_file, out_file, minsup, minconf):
        
    subprocess.call(['java', '-jar', spmf_jar, 'run', 'ERMiner', in_file, out_file, minsup, minconf])
    
    return True


def main_ERMiner(subtype, minsup, minconf):
    subtype_path = ''
    pnames_list = []
    
    if subtype == 'h3n2':
        subtype_path = h3n2_path
        pnames_list = h3n2_pnames_list
        
    elif subtype == 'ph1n1':
        subtype_path = ph1n1_path
        pnames_list = ph1n1_pnames_list
    
    elif subtype == 'eph1n1':
        subtype_path = eph1n1_path
        pnames_list = eph1n1_pnames_list
    
    elif subtype == 'b':
        subtype_path = b_path
        pnames_list = b_pnames_list
    
    outfile_spmf_dir = subtype_path + 'spmf_file/ERMiner/'
    if not os.path.exists(outfile_spmf_dir):
        os.makedirs(outfile_spmf_dir)
            
    for p in pnames_list:
        in_file_path = subtype_path + 'spmf_file/'
        in_file = in_file_path + p + '_spmf.dat'
        out_file = outfile_spmf_dir + p + '_ERMiner.dat'
        ERMiner(in_file, out_file, minsup, minconf)
        
    return True



def TNS(in_file, out_file, k, minconf, delta):
    subprocess.call(['java', '-jar', spmf_jar, 'run', 'TNS', in_file, out_file, k, minconf, delta])
    return True



def main_TNS(subtype, k, minconf, delta):
    subtype_path = ''
    pnames_list = []
    
    if subtype == 'h3n2':
        subtype_path = h3n2_path
        pnames_list = h3n2_pnames_list
        
    elif subtype == 'ph1n1':
        subtype_path = ph1n1_path
        pnames_list = ph1n1_pnames_list
    
    elif subtype == 'eph1n1':
        subtype_path = eph1n1_path
        pnames_list = eph1n1_pnames_list
    
    elif subtype == 'b':
        subtype_path = b_path
        pnames_list = b_pnames_list
    
    outfile_spmf_dir = subtype_path + 'spmf_file/TNS/'
    if not os.path.exists(outfile_spmf_dir):
        os.makedirs(outfile_spmf_dir)
            
    for p in pnames_list:
        in_file_path = subtype_path + 'spmf_file/'
        in_file = in_file_path + p + '_spmf.dat'
        out_file = outfile_spmf_dir + p + '_TNS.dat'
        TNS(in_file, out_file, k, minconf, delta)
        
    return True



def SPAM_AGP(in_file, out_file, minsup):
    subprocess.call(['java', '-jar', spmf_jar, 'run', 'SPAM_AGP', in_file, out_file, minsup])
    return True


def main_SPAM(subtype, minsup):
    subtype_path = ''
    pnames_list = []
    
    if subtype == 'h3n2':
        subtype_path = h3n2_path
        pnames_list = h3n2_pnames_list
        
    elif subtype == 'ph1n1':
        subtype_path = ph1n1_path
        pnames_list = ph1n1_pnames_list
    
    elif subtype == 'eph1n1':
        subtype_path = eph1n1_path
        pnames_list = eph1n1_pnames_list
    
    elif subtype == 'b':
        subtype_path = b_path
        pnames_list = b_pnames_list
    
    outfile_spmf_dir = subtype_path + 'spmf_file/SPAM_AGP/'
    if not os.path.exists(outfile_spmf_dir):
        os.makedirs(outfile_spmf_dir)
            
    for p in pnames_list:
        in_file_path = subtype_path + 'spmf_file/'
        in_file = in_file_path + p + '_spmf.dat'
        out_file = outfile_spmf_dir + p + '_SPAM_AGP.dat'
        SPAM_AGP(in_file, out_file, minsup)
        
    return True


def TKS(in_file, out_file, k):
    subprocess.call(['java', '-jar', spmf_jar, 'run', 'TKS', in_file, out_file, k])
    return True



def main_TKS(subtype, k):
    subtype_path = ''
    pnames_list = []
    
    if subtype == 'h3n2':
        subtype_path = h3n2_path
        pnames_list = h3n2_pnames_list
        
    elif subtype == 'ph1n1':
        subtype_path = ph1n1_path
        pnames_list = ph1n1_pnames_list
    
    elif subtype == 'eph1n1':
        subtype_path = eph1n1_path
        pnames_list = eph1n1_pnames_list
    
    elif subtype == 'b':
        subtype_path = b_path
        pnames_list = b_pnames_list
    
    outfile_spmf_dir = subtype_path + 'spmf_file/TKS/'
    if not os.path.exists(outfile_spmf_dir):
        os.makedirs(outfile_spmf_dir)
            
    for p in pnames_list:
        in_file_path = subtype_path + 'spmf_file/'
        in_file = in_file_path + p + '_spmf.dat'
        out_file = outfile_spmf_dir + p + '_TKS.dat'
        TKS(in_file, out_file, k)
        
    return True



def parse_patterns(in_file, out_file):
    
    out_f = open(out_file, 'w+')
    fieldnames = ['sites', 'support']
    
    writer = csv.DictWriter(out_f, delimiter=',', lineterminator='\n', fieldnames = fieldnames)
    writer.writeheader()
    
    with open(in_file, 'r') as f:
        lines = f.readlines()
        for l in lines:
            components = l.split(' -1 ')
            
            sup = components[-1]
            sup_exact = sup.split(' ')[-1].rstrip()
            # print(sup_exact)
            
            sites = components[0:-1]
            # print(sites)
            if len(sites) == 1:
                pass
            else:
                sites_format = ','.join(sites)
                writer.writerow({'sites': sites_format, 'support': sup_exact})
    
    
    out_f.close()
    
    check_null = open(out_file, 'r')
    n = len(check_null.readlines())
    check_null.close()
    
    if n == 1:
        os.remove(out_file)
        
    return True
    


def main_parse_patterns(subtype, algo):
    
    subtype_path = ''
    pnames_list = []
    
    if subtype == 'h3n2':
        subtype_path = h3n2_path
        pnames_list = h3n2_pnames_list
        
    elif subtype == 'ph1n1':
        subtype_path = ph1n1_path
        pnames_list = ph1n1_pnames_list
    
    elif subtype == 'eph1n1':
        subtype_path = eph1n1_path
        pnames_list = eph1n1_pnames_list
    
    elif subtype == 'b':
        subtype_path = b_path
        pnames_list = b_pnames_list
        
    file_path = subtype_path + 'spmf_file/' + algo + '/'
    outfile_path = file_path + 'results/'
    if not os.path.exists(outfile_path):
        os.makedirs(outfile_path)
    
    for p in pnames_list:
        infile_name = file_path + p + '_' + algo + '.dat'
        outfile_name = outfile_path + p + '_' + algo + '.tab'
        parse_patterns(infile_name, outfile_name)
    
    return True


# Algo for detecting patterns:
# 1. SPAM_AGP
# 2. TKS

# Algo for detecting rules: 
# 1. ERMiner
# 2. TNS


def parse_rules_with_conf(in_file, out_file):
    
    out_f = open(out_file, 'w+')
    fieldnames = ['target', 'source', 'support', 'confidence']
    
    writer = csv.DictWriter(out_f, delimiter=',', lineterminator='\n', fieldnames = fieldnames)
    writer.writeheader()
    
    with open(in_file, 'r') as f:
        lines = f.readlines()
        for l in lines:
            components = l.split('#')

            rules = components[0]
            sup = components[1]
            conf = components[2].rstrip()
            
            r = rules.split(' ==> ')
            sources = r[0]
            target = r[1].rstrip()
            
            sup_exact = sup.split(' ')[1]
            conf_exact = conf.split(' ')[1]
            
            writer.writerow({'target': str(target), 'source': str(sources),  
                             'support': sup_exact, 'confidence': conf_exact})
    
    out_f.close()
    
    check_null = open(out_file, 'r')
    n = len(check_null.readlines())
    # print(n)
    check_null.close()
    
    if n == 1:
        os.remove(out_file)
        
    return True



def main_parse_rules_with_conf(subtype, algo):
    subtype_path = ''
    pnames_list = []
    
    if subtype == 'h3n2':
        subtype_path = h3n2_path
        pnames_list = h3n2_pnames_list
        
    elif subtype == 'ph1n1':
        subtype_path = ph1n1_path
        pnames_list = ph1n1_pnames_list
    
    elif subtype == 'eph1n1':
        subtype_path = eph1n1_path
        pnames_list = eph1n1_pnames_list
    
    elif subtype == 'b':
        subtype_path = b_path
        pnames_list = b_pnames_list
        
    file_path = subtype_path + 'spmf_file/' + algo + '/'
    outfile_path = file_path + 'results/'
    if not os.path.exists(outfile_path):
        os.makedirs(outfile_path)
    
    for p in pnames_list:
        infile_name = file_path + p + '_' + algo + '.dat'
        outfile_name = outfile_path + p + '_' + algo + '.csv'
        parse_rules_with_conf(infile_name, outfile_name)
    
    return True



def split_rules(in_file, out_file):
    
    out_f = open(out_file, 'w+')
    fieldnames = ['target', 'source', 'support']
    
    writer = csv.DictWriter(out_f, delimiter=',', lineterminator='\n', fieldnames = fieldnames)
    writer.writeheader()
    
    with open(in_file, 'r') as f:
        lines = f.readlines()
        for l in lines:
            components = l.split('#')
            
            rules = components[0]
            sup = components[1]
            conf = components[2].rstrip()
            
            r = rules.split(' ==> ')
            sources = r[0]
            target = r[1].rstrip()
            
            sup_exact = sup.split(' ')[1]
            conf_exact = conf.split(' ')[1]
            
            sites_list = sources.split(',')
            targets_list = target.split(',')
            
            for s in sites_list:
                for t in targets_list:
                    writer.writerow({'target': t, 'source': s, 'support': sup_exact})

            
    out_f.close()
    
    check_null = open(out_file, 'r')
    n = len(check_null.readlines())
    # print(n)
    check_null.close()
    
    if n == 1:
        os.remove(out_file)
        
    return True



def sort_csv(in_file, out_file):
    
    if os.path.isfile(in_file):
        f = open(in_file, 'r')
        records = csv.DictReader(f, delimiter=',', lineterminator='\n')
        sortedlist = sorted(records, key=lambda row:(row['target'],row['source']), reverse=False)
        f.close()
        
        with open(out_file, 'w') as f:
            fieldnames = ['target', 'source', 'support']
            writer = csv.DictWriter(f, delimiter=',', lineterminator='\n', fieldnames=fieldnames)
            writer.writeheader()
            for row in sortedlist:
                writer.writerow(row)
        os.remove(in_file)
        return True
    
    else:
        pass
        # print('No rules extracted ... ')
        


def remove_duplicates(in_file, out_file):
    
    if os.path.isfile(in_file):
        data = pd.read_csv(in_file)
        
        target_list = []
        source_list = []
        support_list = []
        
        # clean_data = {'target': target_list, 'source': source_list, 'support': support_list}
            
        i = 0 
        while i < len(data.index):
            # print('i = ' + str(i))
            
            target = data.iloc[i]['target']
            source = data.iloc[i]['source']
            support = data.iloc[i]['support']
            
            j = i + 1
            while j < len(data.index):
                
                ref_target = data.iloc[j]['target']
                ref_source = data.iloc[j]['source']
                ref_support = data.iloc[j]['support']
                
                flag1 = (ref_target == target)
                flag2 = (ref_source == source)
                
                if flag1 == True and flag2 == True:
                    # support = support + ref_support
                    support = max([support, ref_support])
                    j = j + 1
                    
                else: 
                    # print(support)
                    break
            
            i = j
            target_list.append(target)
            source_list.append(source)
            support_list.append(support)
            
            
        clean_data = {'target': target_list, 'source': source_list, 'support': support_list}
        clean_df = pd.DataFrame(clean_data)
        output = clean_df[['target', 'source', 'support']]
        output.to_csv(out_file, encoding='utf-8', index=False)
        
        os.remove(in_file)
        return True
        
    else:
        pass
        # print('No sorted file ---------- ')
        



def split_rules2cyto(in_file, out_file):
    
    # print('split_rules ------------- ')
    split_rules(in_file, 'tmp_split.csv')
    
    # print('sort_csv ------------- ')
    sort_csv('tmp_split.csv', 'tmp_sort.csv')
    
    # print('remove_duplicates ------------- ')
    remove_duplicates('tmp_sort.csv', out_file)
    
    return True



def main_split_rules2cyto(subtype, algo):
    subtype_path = ''
    pnames_list = []
    
    if subtype == 'h3n2':
        subtype_path = h3n2_path
        pnames_list = h3n2_pnames_list
        
    elif subtype == 'ph1n1':
        subtype_path = ph1n1_path
        pnames_list = ph1n1_pnames_list
    
    elif subtype == 'eph1n1':
        subtype_path = eph1n1_path
        pnames_list = eph1n1_pnames_list
    
    elif subtype == 'b':
        subtype_path = b_path
        pnames_list = b_pnames_list
        
    file_path = subtype_path + 'spmf_file/' + algo + '/'
    outfile_path = file_path + 'results/'
    if not os.path.exists(outfile_path):
        os.makedirs(outfile_path)
    
    for p in pnames_list:
        infile_name = file_path + p + '_' + algo + '.dat'
        outfile_name = outfile_path + p + '_' + algo + '_cyto.csv'
        split_rules2cyto(infile_name, outfile_name)
    
    return True


# main_split_rules2cyto('ph1n1', 'TNS')
# main_split_rules2cyto('eph1n1', 'TNS')
# main_split_rules2cyto('h3n2', 'TNS')
# main_split_rules2cyto('b', 'TNS')

