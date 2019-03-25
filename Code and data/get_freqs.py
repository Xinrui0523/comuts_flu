#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import csv
from collections import Counter


# In[2]:


def write_selected_freqs(input_file, input_site, output_file):
    # data_path = 'C:/Users/YU007/project/CRISPR/data/'
    
    df = pd.read_csv(input_file,delimiter='\t',header=None)
    
    year = df[:][0]
    site = df[:][1]
    residue = df[:][2]
    freq = df[:][4]
    site_output = []
    freq_output = []
    year_output = []
    residue_output = []
    
    for j in range(len(site)):
        for i in range(len(input_site)):
            if site[j] == input_site[i]:
                site_output.append(site[j])
                year_output.append(year[j])
                residue_output.append(residue[j])
                freq_output.append(freq[j])
    x_aixs_label0 = []
    
    for i in range(len(site_output)):
        x1 = str(site_output[i])
        x2 = str(residue_output[i])
        x_aixs_label0.append(x2+x1)
    
    x_aixs_label = list(set(x_aixs_label0))
    year_num = len(set(year_output))
    year_count0 = Counter(year_output)
    year_count=list(year_count0.items())
    result = np.zeros((year_num,len(x_aixs_label)+1))
    np.set_printoptions(suppress=True)
    
    k = 0
    for i in range(len(year_output)):
        s = str(site_output[i])
        r = str(residue_output[i])
        c = x_aixs_label.index(r+s)
        result[k][0] = int(year_output[i])
        result[k][c+1] = freq_output[i]
        if (i<len(year_output)-1) and (year_output[i]!=year_output[i+1]):       
            k += 1
    rstname = ['year']+x_aixs_label
    
    f = open(output_file, 'w')
    
    for i in range(year_num+1):
        if i == 0:
            for j in range(len(x_aixs_label)+1):
                if j == 0:
                    f.write(str(rstname[j]))
                else:
                    f.write(',' + str(rstname[j]))
            f.write('\n')
        else:
            for j in range(len(x_aixs_label)+1):
                if j == 0:
                    f.write(str(result[i-1][j]))
                else:
                    f.write(',' + str(result[i-1][j]))
            f.write('\n')               
    f.close()
    
    return True


# In[5]:


ph1n1_path = './H1N1/pH1N1/year_seqs/'
ph1n1_na_list = [41, 369]
# ph1n1_ha_list = [185, 451, 69, 143, 216, 260]
ph1n1_ha_list = [202, 468, 86, 160, 233, 277]

ph1n1_na_input = ph1n1_path + 'NA_year/allfrequency.txt'
ph1n1_na_output = ph1n1_path + 'NA_year/hana_na_freq.csv'

ph1n1_ha_input = ph1n1_path + 'HA_year/allfrequency.txt'
ph1n1_ha_output = ph1n1_path + 'HA_year/hana_ha_freq.csv'


# In[6]:


write_selected_freqs(ph1n1_ha_input, ph1n1_ha_list, ph1n1_ha_output)
# write_selected_freqs(ph1n1_na_input, ph1n1_na_list, ph1n1_na_output)


# In[7]:


h3n2_path = './H3N2/year_seqs/'
h3n2_ha_input = h3n2_path + 'HA_year/allfrequency.txt'
h3n2_ha_output = h3n2_path + 'HA_year/ha_135.csv'

write_selected_freqs(h3n2_ha_input, [151], h3n2_ha_output)


# In[8]:


h3n2_na_input = h3n2_path + 'NA_year/allfrequency.txt'
h3n2_na_output = h3n2_path + 'NA_year/na_400.csv'

write_selected_freqs(h3n2_na_input, [400], h3n2_na_output)


# In[3]:


ph1n1_path = './H1N1/pH1N1/year_seqs/'
ph1n1_pb2_list = [54, 66, 195, 293, 731, 344, 354, 299, 456]

input_file = ph1n1_path + 'PB2_year/allfrequency.txt'
output_file = ph1n1_path + 'PB2_year/selected_frequency.csv'
# write_selected_freqs(input_file, ph1n1_pb2_list, output_file)


# In[4]:


ph1n1_ns1_list = [2, 125]
input_file = ph1n1_path + 'NS1_year/allfrequency.txt'
output_file = ph1n1_path + 'NS1_year/selected_frequency.csv'
# write_selected_freqs(input_file, ph1n1_ns2_list, output_file)


# In[5]:


eph1n1_path = './H1N1/epH1N1/year_seqs/'
eph1n1_pb1_list = [177, 327, 372, 375, 383]
eph1n1_pb1N40_list = [138, 288, 333, 336, 344]

pb1_input = eph1n1_path + 'PB1_year/allfrequency.txt'
pb1_output = eph1n1_path + 'PB1_year/selected_frequency.csv'

pb1N40_input = eph1n1_path + 'PB1-N40_year/allfrequency.txt'
pb1N40_output = eph1n1_path + 'PB1-N40_year/selected_frequency.csv'

write_selected_freqs(pb1_input, eph1n1_pb1_list, pb1_output)
write_selected_freqs(pb1N40_input, eph1n1_pb1N40_list, pb1N40_output)


# In[11]:


# eph1n1_ha_list = [82, 94, 208, 266, 415]
eph1n1_ha_list = [99, 111, 225, 283, 432]

ha_input = eph1n1_path + 'HA_year/allfrequency.txt'
ha_output = eph1n1_path + 'HA_year/selected_frequency.csv'
write_selected_freqs(ha_input, eph1n1_ha_list, ha_output)


# In[12]:


h3n2_path = './H3N2/year_seqs/'
# h3n2_ha = [164, 174, 193, 201]
h3n2_ha = [180, 190, 209, 217]

ha_input = h3n2_path + 'HA_year/allfrequency.txt'
ha_output = h3n2_path + 'HA_year/selected_frequency.csv'
write_selected_freqs(ha_input, h3n2_ha, ha_output)


# In[13]:


b_path = './B/_b_all_lineages/year_seqs/'
# b_ha = [108, 150, 166, 230]
b_ha = [123, 165, 181, 245]
b_input = b_path + 'HA_year/allfrequency.txt'
b_output = b_path + 'HA_year/selected_frequency.csv'
write_selected_freqs(b_input, b_ha, b_output)


# In[ ]:




