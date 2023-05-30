import io
import os
import pandas as pd
from functools import partial

def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

    #transforms the INFO column into new columns
def info_to_columns(vcf):
    #divide each element in the original string by the separator ";"
    #and then divide each element by the separator "="
    new_df = pd.DataFrame([dict(kv.split("=") for kv in s.split(";")) for s in vcf['INFO']])
    return pd.concat([vcf.drop('INFO', axis=1), new_df], axis=1)


def apply_functions(df, func_dict):
    for key in func_dict:
        df = df.apply(partial(eval(key),func_dict[key]["n"]), axis=func_dict[key]["axis"])
    return df


#stores the vcf dataframe frequencies (AF) into a dictionary
#having as key frequency_list
#and having as key for total the number of cells n
def vcf_frequencies_to_dict(vcf,n):
    return {"frequency_list": [float(i) for i in vcf["AF"].tolist()],"total":n}
    #return {"frequency_list":vcf["AF"].tolist(),"total":n}
    


"""
read_vcf("C:\\Users\\jm\\OneDrive\\tfm bioinformatica\\vcf\\MCR01.A03.C01.HaplotypeCaller.joint.vcf")



#take only the n rows of the vcf file which fullfill the condition given by the function f
#if n=0 it takes all the rows
def take_rows(vcf, f, n=0):
    vcf.apply(f, axis=1)
    if n == 0:
        return vcf[vcf.apply(f, axis=1)]
    else:
        return vcf[vcf.apply(f, axis=1)].head(n)


function = lambda x: len(x["REF"]) == 1 and len(x["ALT"]) == 1
f3 = lambda n,df: df.head(n[0])

def apply_functions(df, func_dict):
    for key in func_dict:
        df = df[df.apply(func_dict[key], axis=1)]
    return pd.DataFrame(result)

def apply_functions(df, func_list):
    for f in func_list:
        df = df[df.apply(f, axis=1)]
    return df

def apply_functions(df, func_dict):
    for key in func_dict:
        df = df[df.apply(partial(eval(key),func_dict[key]["n"]), axis=func_dict[key]["axis"])]
    return df



apply_functions(data,{"lambda n,df: df.head(n[0])": {"n":[2],"axis":1}})


def take_columns(vcf, n, f):
    columns = []
    for i in range(n):
        columns.append([])
    for index, row in vcf.iterrows():
        if f(row):
            for i in range(n):
                columns[i].append(row[i])
    return columns
"""


