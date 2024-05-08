#! /usr/bin/env python

import requests, sys
import csv

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

def add_pheno_to_row(fc,i,row,header_row,pval_bound,zscore_bound):
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
    elif abs(float(row[header_row.index('zScore')])) > zscore_bound and float(row[header_row.index('pValue')]) < pval_bound:
        output_pheno = get_phenotype(row[header_row.index('gene')])
    
    if fc:
        row.insert(-2,output_pheno)
    else:
        row.append(output_pheno)
    return(row)

def add_phenotype_to_file(file_in, fc=0, pval_bound=0.01, zscore_bound=3):
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
                newrow = add_pheno_to_row(fc,i,row[0:-1].split(","),header_row,pval_bound,zscore_bound)
                writer.writerow(newrow)

#file='/Users/lbrussel/rnaseq-emc/notebook/rnaseq-voila-main/downloads/umcu_rnaseq_fib_untreated_res_outrider_genes_counts_copy.tsv'
file='/Users/lbrussel/rnaseq-emc/notebook/rnaseq-voila-main/outrider/umcu_rnaseq_fib_untreated_res_outrider_genes_counts_nopheno.tsv'
add_phenotype_to_file(file, 0, 0.01, 2.5)
#print(get_phenotype("ENSG00000255495"))