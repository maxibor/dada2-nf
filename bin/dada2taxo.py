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

    args = parser.parse_args()

    infile = args.dada_csv
    return(infile)


def get_otu(serie):
    if "x.Species" in serie.index:
        specname = " ".join(list(serie)[-2:])
        tax = ncbi.get_name_translator([specname])[specname][0]
    else:
        tax = ncbi.get_name_translator([serie[-1]])[serie[-1]][0]
    return(tax)


def dada_to_taxo(dada_df):
    nrow = dada_df.shape[0]
    ncol = dada_df.shape[1]
    res = {}
    for i in range(0, nrow):
        dadatax = dada_df.iloc[i, :][:-1].dropna()
        try:
            tax = get_otu(dadatax)
            if tax not in list(res.keys()):
                res[tax] = dada_df.iloc[i, ncol-1]
            else:
                res[tax] += dada_df.iloc[i, ncol-1]
        except (KeyError, IndexError):
            print("Key Error")
            pass

    resdf = pd.DataFrame(list(res.items()), columns=[
        'TAXID', "sample"], index=list(res.keys())).drop('TAXID', 1)
    return(resdf)


def get_basename(inputFile):
    return(inputFile.split("/")[-1].split('.')[0])


if __name__ == "__main__":
    INFILE = _get_args()
    outname = get_basename(INFILE)+".dadataxo.csv"

    ncbi = NCBITaxa()

    d = pd.read_csv(INFILE, index_col=0)

    r = dada_to_taxo(d)

    r.to_csv(outname)
