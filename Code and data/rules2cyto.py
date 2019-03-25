#!/usr/bin/env python
# coding: utf-8

import os
import csv
import pandas as pd



def format_rules(infile_prefix):
    infile = infile_prefix + '.rules'
    outfile_prefix = infile_prefix + '_cy'
    outfile = outfile_prefix + '.csv'
    
    f = open(outfile, 'w')
    fieldnames = ['target', 'source', 'support']
    
    writer = csv.DictWriter(f, delimiter=',', lineterminator='\n', fieldnames= fieldnames)
    writer.writeheader()
    
    with open(infile, 'r') as f:
        lines = f.readlines()
        for l in lines:
            # (confidence, frequency) target <= sources (support)
            components = l.split(' ')
            
            c = components[0]
            
            i1 = c.index('(')
            i2 = c.index(',')
            confidence = c[i1+1:i2]
            
            target = components[1]
            
            s = components[-1]
            is1 = s.index('(')
            is2 = s.index(')')
            support = s[is1+1:is2]
            # print(support)
            
            rule_flag = components.index('<=')
            sources = components[rule_flag+1:-1]
            
            for s in sources:
                # target, source, confidence, support
                writer.writerow({'target': target, 'source': s, 'support': support})
                
    f.close()

    return outfile_prefix
    


def sort_csv(in_csv):
    # Sorted the target and source columns
    f = open(in_csv+'.csv', 'r')
    records = csv.DictReader(f, delimiter=',', lineterminator='\n')
    sortedlist = sorted(records, key=lambda row:(row['target'],row['source']), reverse=False)
    f.close()
    
    outfile_prefix = in_csv + 't'
    outfile = outfile_prefix + '.csv'
    with open(outfile, 'w') as f:
        fieldnames = ['target', 'source', 'support']
        writer = csv.DictWriter(f, delimiter=',', lineterminator='\n', fieldnames=fieldnames)
        writer.writeheader()
        for row in sortedlist:
            writer.writerow(row)
    
    
    os.remove(in_csv+'.csv')
    
    return outfile_prefix



def add_duplicates(in_csv):
    
    data = pd.read_csv(in_csv+'.csv')
    
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
    if len(target_list) > 0:
        output.to_csv(in_csv + 'o.csv', encoding='utf-8', index=False)
    
    os.remove(in_csv+'.csv')
    
    return True
    


def rules2cyto(infile_prefix):
    
    format_file = format_rules(infile_prefix)
    sorted_file = sort_csv(format_file)
    add_duplicates(sorted_file)
    
    return True
