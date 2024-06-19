#! /usr/bin/env python

import requests, sys
import csv
import argparse
import gzip, shutil

parser = argparse.ArgumentParser(description="Add phenotype to outrider output", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("file", help="provide outrider input file to add phenotype")
parser.add_argument("-p", "--pval", type=float, help="Max p-value to search phenotype for")
parser.add_argument("-z", "--zscore", type=float, help="Min z-score to search phenotype for")
parser.add_argument("-m", "--meancorr", type=int, help="Min mean corrected to search phenotype for")

def get_phenotype(gene):
    '''Returns ensemble phenotype'''
    server = "https://rest.ensembl.org"
    ext = "/phenotype/gene/homo_sapiens/"
    output_pheno=""
    output_omim=""
    # try:
    r = requests.get(server+ext+gene, headers={ "Content-Type" : "application/json"})
    # except:
    if not r.ok:
        output_pheno = "error in request"#r.raise_for_status()
        # sys.exit()
    else:
        decoded = r.json()
        for i in decoded:
            if i['source']=='MIM morbid':
                output_omim += i['source']+':'+i['description']+";" 
            output_pheno += (i['source']+':'+i['description']+";")
        if output_omim: output_pheno = output_omim
        elif not output_pheno: output_pheno = "nothing found"
    return output_pheno.replace(",","-").replace('\"',"")

def add_pheno_to_row(fc,i,row,header_row,pval_bound,zscore_bound,meancorr_bound):
    '''get phenotype
    add as new column to row
    and return as newrow'''
    output_pheno="outside boundaries"
    if i==0:
        output_pheno="phenotype"
    #for feature count input
    elif fc: 
        output_pheno = get_phenotype(row[header_row.index('gene_name')])
    #for outrider result input
    elif abs(float(row[header_row.index('zScore')])) > zscore_bound and float(row[header_row.index('pValue')]) < pval_bound and float(row[header_row.index('meanCorrected')]) >= meancorr_bound:
        output_pheno = get_phenotype(row[header_row.index('gene')])
    
    if fc:
        row.insert(-2,output_pheno)
    else:
        row.append(output_pheno)
    return(row)

def add_phenotype_to_file(file_in, fc=0, pval_bound=0.01, zscore_bound=3, meancorr_bound=0):
    '''Add phenotype as a column to imported file
    for rows that match boundaries
    and write in new file'''
    ext = file_in.split(".",2)[-1]
    file_out = file_in.split(".")[0]+"_pheno_added."+ext
    with open(file_out, 'w') as newcsvfile:
        writer = csv.writer(newcsvfile)
        with open(file_in, newline='') as csvfile:
            reader = csvfile.readlines()
            header_row = reader[0][0:-1].split(",")
            for i,row in enumerate(reader):
                newrow = add_pheno_to_row(fc,i,row[0:-1].split(","),header_row,pval_bound,zscore_bound,meancorr_bound)
                writer.writerow(newrow)
    with open(file_out, 'rb') as f_out:
        with gzip.open(file_out+".gz", 'wb') as f_outgz:
            shutil.copyfileobj(f_out, f_outgz)
    print("done")

def main():
    args = parser.parse_args()
    file = args.file
    pval = 0.01
    zscore = 2.5
    meancorr = 0
    if args.pval is not None: pval = args.pval
    if args.zscore is not None: zscore = args.zscore
    if args.meancorr is not None: meancorr = args.meancorr
    
    add_phenotype_to_file(file, 0, pval, zscore, meancorr)

main()