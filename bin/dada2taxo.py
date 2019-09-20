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
    parser.add_argument(
        '-r',
        default='Species',
        help='Rank to keep (Genus | Species). Default = Species'
    )

    args = parser.parse_args()

    infile = args.dada_csv
    sname = args.s
    rank = args.r
    return(infile, sname, rank)


def get_taxid(specname):
    taxid_d = ncbi.get_name_translator([specname])
    if len(taxid_d.keys()) > 0:
        specname = list(taxid_d.keys())[0]
        taxid = taxid_d[specname]
        res = taxid[0]
    else:
        res = 0

    return(res)


def dada_to_taxo(dada_df, sample_name, rank):
    d = pd.read_csv(dada_df, index_col=0)
    if rank == 'Species':
        d = d.dropna(axis = 0)
        d['name'] = d['Genus']+' '+d['Species']
    elif rank == 'Genus':
        d['name'] = d['Genus']
        d = d.dropna(subset=['Genus'], axis=0)
    d['TAXID'] = d['name'].apply(get_taxid)
    dct = {}
    for i in list(d.index):
        if d['TAXID'][i] not in list(dct.keys()):
            dct[d['TAXID'][i]] = d['abundance'][i]
        else:
            dct[d['TAXID'][i]] += d['abundance'][i]
    s = pd.Series(dct).to_frame(name=sample_name)
    return(s)


def get_basename(inputFile):
    return(inputFile.split("/")[-1].split('.')[0])


if __name__ == "__main__":
    INFILE, SNAME, RANK = _get_args()
    outname = get_basename(INFILE)+".dadataxo.csv"

    ncbi = NCBITaxa()

    r = dada_to_taxo(INFILE, SNAME, RANK)

    r.to_csv(outname, index=True)