import logging
from Bio import Entrez
import pandas as pd
logger = logging.getLogger('The Reticulator')
logger.setLevel(logging.INFO)

Entrez.email = "b.temperton@exeter.ac.uk"

BALTIMORE_TYPES = ['ssDNA viruses', 'dsDNA viruses, no RNA stage']


def main():
    records = {}
    df = pd.read_csv('/Users/bt273/Downloads/BAV.viruses.csv', sep=',')
    count = 1230
    total = len (df.index)
    try:
        for i in df.contig.values:
            records[i] = getTaxonomy(i)
            count +=1
            print "processed %i of %i (%s)" % (count, total, i)

    except Exception as e:
        tax_df = pd.DataFrame.from_dict(records, orient='index')
        tax_df.to_csv('/Users/bt273/Dropbox/TempertonLab/ViralEcologyToolkit/data/BAV.tax.%i.csv' % count,
                          header=True)
        raise e

    tax_df = pd.DataFrame.from_dict(records, orient='index')
    tax_df.to_csv('/Users/bt273/Dropbox/TempertonLab/ViralEcologyToolkit/data/BAV.tax.%i.csv' % count,
                  header=True)


def getTaxonomy(name):

    rtnValue = {}

    handle = Entrez.esearch(db="Taxonomy", term=name)
    record = Entrez.read(handle)
    handle.close()
    if record['Count'] == '0':
        print 'Found no details for %s' % name
        return rtnValue

    tax_id = record['IdList'][0]
    handle = Entrez.efetch(db="Taxonomy", id=tax_id)
    record = Entrez.read(handle)
    handle.close()

    lineage = record[0]['LineageEx']
    for l in lineage:
        try:
            if l['Rank'] == 'family':
                rtnValue['family'] = l['ScientificName']
            elif l['Rank'] =='genus':
                rtnValue['genus'] = l['ScientificName']
            elif l['Rank'] == 'order':
                rtnValue['order'] = l['ScientificName']
            elif l['Rank'] == 'no rank':
                if l['ScientificName'] in BALTIMORE_TYPES:
                    rtnValue['baltimore'] = l['ScientificName']
        except KeyError:
            pass

    return rtnValue

if __name__ == '__main__':
    main()