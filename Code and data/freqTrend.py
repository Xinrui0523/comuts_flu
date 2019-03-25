
# coding: utf-8

# In[47]:


from Bio import AlignIO
from Bio.Align import AlignInfo
from collections import Counter


# In[66]:


def get_frequency(file_path, start_year, end_year):
    
    out_freq_file = open(file_path + 'allfrequency.txt', 'w+')
    # out_freq_file.write('year, site_number, residue, #residue, frequency, #sequences')
    
    for i in range(start_year, end_year+1):
        seq_file_name = file_path + str(i)+'.fasta'
        alignment = AlignIO.read(seq_file_name, 'fasta')
        num_seqs = len(alignment)
        
        for j in range(alignment.get_alignment_length()):
            res_j = alignment[:, j]
            counts = Counter(res_j)
            freqs = {}
            
            for key, value in counts.items():
                f = value/num_seqs
                freqs[key] = f
                
                out_freq_file.write('%i\t%i\t%s\t%i\t%f\t%i\n' % (i, j+1, key, counts[key], freqs[key], num_seqs))            
    
    return True


# In[83]:


# h3n2_path = './H3N2/year_seqs/NS1_year/'
# get_frequency(h3n2_path, 1988, 2018)


# In[84]:


def main_get_frequency(subtype):
    h3n2_path = './H3N2/year_seqs/'
    ph1n1_path = './H1N1/pH1N1/year_seqs/'
    eph1n1_path = './H1N1/epH1N1/year_seqs/'
    b_path = './B/year_seqs/'
    
    h3n2_pnames_list = ['HA', 'M1', 'M2', 'NA', 'NP', 'HA_NA', 'NS2', 'PA', 'PA-N155', 
                        'PA-N182', 'PA-X', 'PB1', 'PB1-F2', 'PB1-N40', 'PB2']
    ph1n1_pnames_list = ['HA', 'M1', 'M2', 'NA', 'NP', 'NS1', 'NS2', 'PA', 'PA-N155', 
                         'PA-N182', 'PA-X', 'PB1', 'PB1-N40', 'PB2', 'HA_NA']
    eph1n1_pnames_list = ['HA', 'M1', 'M2', 'NA', 'NP', 'NS1', 'NS2', 'PA', 'PA-N155', 
                          'PA-N182', 'PA-X', 'PB1', 'PB1-F2', 'PB1-N40', 'PB2', 'HA_NA']
    b_pnames_list = ['PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'NB', 
                     'M1', 'BM2', 'NS1', 'NS2', 'HA_NA']
    
    h3n2_year_range = [1968, 2018]
    ph1n1_year_range = [2009, 2017]
    eph1n1_year_range = [1995, 2009]
    b_year_range = [1993, 2016]
    
    subtype_path = ''
    if subtype == 'h3n2':
        subtype_path = h3n2_path
        year_range = h3n2_year_range
        
        for p in h3n2_pnames_list:
            file_path = subtype_path + p + '_year/'
            get_frequency(file_path, year_range[0], year_range[1])
    
    elif subtype == 'ph1n1':
        subtype_path = ph1n1_path
        year_range = ph1n1_year_range
        for p in ph1n1_pnames_list:
            file_path = subtype_path + p + '_year/'
            get_frequency(file_path, year_range[0], year_range[1])
    
    elif subtype == 'eph1n1':
        subtype_path = eph1n1_path
        year_range = eph1n1_year_range
        for p in eph1n1_pnames_list:
            file_path = subtype_path + p + '_year/'
            get_frequency(file_path, year_range[0], year_range[1])
    
    elif subtype == 'b':
        subtype_path = b_path 
        year_range = b_year_range
        for p in b_pnames_list:
            file_path = subtype_path + p + '_year/'
            get_frequency(file_path, year_range[0], year_range[1])
            
    return True
            


# In[85]:


# main_get_frequency('h3n2')
# main_get_frequency('ph1n1')
# main_get_frequency('eph1n1')
# main_get_frequency('b')


# We first find out the sites under positive selection and obtain the sites co-evolving with the positive-selection sites, which would be predicted as the sites to be mutated. In [29], the positive selection site is defined as “a site that has been mutated between successive years and then remains fixed in the population for at least 1 year”. 
# 
# pos_sites: the sites occurring frequently in our rules (i.e. the sites co-mutated frequently with other sites, or with large in-degree) as the sites to be mutated.

# In[ ]:


def plot_frequency_across_years(freq_file, site_list):
    # if more than two years, the frequency is higher than low level, it's considered to be the candidates.
    
    # using all the sequences before the predicting year to train
    
    with open(freq_file, 'r') as f:
        for line in f:
            items = line.split('\t')
            year = items[0]
            site = items[1]
            res = items[2]
            count = items[3]
            freq = items[4]
            total = items[5]
            
            if total <= 100:
                low_level = 0.05
            else:
                low_level = 0.03
            
            for i in range(start_year + cut, pred_year):
                for j in range(pro_length):
                    if freq >= low_level:
                        for k in range(0, cut):
                            
            
        
        
    


# In[57]:


def get_pos_sites():
    
    pos_sites = []
    
    return pos_sites


# In[ ]:


def get_coevolve_sites(pos_sites):
    
    coevolve_sites = []
    
    return coevolve_sites


# In[ ]:


def pred_muts(pos_sites, coevolve_sites):
    
    pred_sites = []
    
    return pred_sites

