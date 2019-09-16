#!/usr/bin/env python

import argparse
import pandas as pd
from ete3 import NCBITaxa


def _get_args():
    '''This function parses and return arguments passed in'''
    parser = argparse.ArgumentParser(
        prog='dada2taxo',
        description='Dada output to taxonomy')
    parser.add_argument('dada_csv', help="path to dada csv file")
    parser.add_argument(
        '-s',
        default=None,
        help='sample name'
    )

    args = parser.parse_args()

    infile = args.dada_csv
    sname = args.s
    return(infile, sname)


def get_taxid(specname):
    taxid_d = ncbi.get_name_translator([specname])

    if len(taxid_d.keys()) > 0:
        specname = taxid_d.keys()[0]
        taxid = taxid_d[specname]
        res = taxid
    else:
        res = 0

    # print(specname, res)
    return(res)


def dada_to_taxo(dada_df, sample_name):
    dada_df = dada_df.drop("Unnamed: 0", axis=1)
    dada_df['specname'] = dada_df.index
    # print(dada_df)
    
    taxdict = {}
    notNA = pd.Series(dada_df.index, index = dada_df.index).str.contains('NA').index
    # print(notNA)
    d = dada_df.drop(notNA, axis=0)
    # d = dada_df
    d['TAXID'] = dada_df['specname'].apply(get_taxid)
    # print(d)
    for i in pd.Series(d.index):
        if d['TAXID'][i] not in taxdict.keys():
            taxdict[d['TAXID'][i]] = d['Freq'][i]
        else:
            taxdict[d['TAXID'][i]]  += d['Freq'][i]
    res = pd.Series(taxdict).to_frame(name=sample_name).fillna(0)
    print(res)
    return(res)


def get_basename(inputFile):
    return(inputFile.split("/")[-1].split('.')[0])


if __name__ == "__main__":
    INFILE, SNAME = _get_args()
    outname = get_basename(INFILE)+".dadataxo.csv"

    ncbi = NCBITaxa()

    d = pd.read_csv(INFILE, index_col=1)

    r = dada_to_taxo(d, SNAME)

    r.to_csv(outname)