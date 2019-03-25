import subprocess
import os


spmf_jar = './jar_SPMF/spmf.jar'

h3n2_path = './H3N2/'
ph1n1_path = './H1N1/pH1N1/'
eph1n1_path = './H1N1/epH1N1/'
b_path = './B/_b_all_lineages/'
# vicb_path = './B/b_vic/'
# yamb_path = './B/b_yam/'

h3n2_pnames_list_ori = ['HA', 'M1', 'M2', 'NA', 'NP', 'NS1', 'NS2', 'PA', 'PA-N155', 'PA-N182', 'PB1', 'PB1-F2', 'PB1-N40', 'PB2', 'HA_NA']
ph1n1_pnames_list_ori = ['HA', 'M1', 'M2', 'NA', 'NP', 'NS1', 'NS2', 'PA', 'PA-N155', 'PA-N182', 'PA-X', 'PB1', 'PB1-N40', 'PB2', 'HA_NA']
eph1n1_pnames_list_ori = ['HA', 'M1', 'M2', 'NA', 'NP', 'NS1', 'NS2', 'PA', 'PA-N155', 'PA-N182', 'PA-X', 'PB1', 'PB1-F2', 'PB1-N40', 'PB2', 'HA_NA']
b_pnames_list_ori = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'NB', 'M1', 'BM2', 'NS1', 'NS2', 'HA_NA']

def main_algo(subtype, algo, params):
    
    subtype_path = ''
    pnames_list = []
    
    if subtype == 'h3n2':
        subtype_path = h3n2_path
        pnames_list = h3n2_pnames_list_ori
        
    elif subtype == 'ph1n1':
        subtype_path = ph1n1_path
        pnames_list = ph1n1_pnames_list_ori
    
    elif subtype == 'eph1n1':
        subtype_path = eph1n1_path
        pnames_list = eph1n1_pnames_list_ori
    
    elif subtype == 'b':
        subtype_path = b_path
        pnames_list = b_pnames_list_ori
    
    outfile_spmf_dir = subtype_path + 'spmf/' + algo + '/'
    if not os.path.exists(outfile_spmf_dir):
        os.makedirs(outfile_spmf_dir)
            
    for p in pnames_list:
        print(p)
        in_file_path = subtype_path + 'spmf/'
        in_file = in_file_path + p + '_spmf.dat'
        out_file = outfile_spmf_dir + p + '_' + algo + '.dat'
        
        param_list = ['java', '-jar', spmf_jar, 'run', algo, in_file, out_file] + params
        subprocess.call(param_list)
    
    return True
