## Collection of miscellaneous functions.
##
##
import pandas as pd
import re
import requests
import time
from bs4 import BeautifulSoup
from mygene import MyGeneInfo



def requests_get(url, interval=1.0, n_trials=5, **kwargs):
    '''
    This function sends a GET request to a given URL.

    Arguments
    ----------
    url (required)
        The URL of the request.

    interval
        The interval time between trials.
        Default: 1.0

    n_trials
        The number of trials.
        Default: 5

    kwargs
        Parameters of requests.get().
    '''
    for _ in range(n_trials):
        try:
            response = requests.get(url, **kwargs)
            response.raise_for_status()
            return response
        except:
            time.sleep(interval)

    response.raise_for_status()



def get_pathways(organism, **kwargs):
    '''
    This function returns the KEGG pathways for a given organism.

    Arguments
    ----------
    organism (required)
        The KEGG organism (e.g. hsa, mmu).

    kwargs
        Parameters of requests_get().
    '''
    response = requests_get('https://rest.kegg.jp/list/pathway/' + organism, **kwargs)
    return re.findall(organism + r'[0-9]+(?=\t)', response.text)



def download_KGML(pathway, save_to='.', encoding='utf-8', **kwargs):
    '''
    This function downloads the KGML of a given pathway.

    Arguments
    ----------
    pathway (required)
        The KEGG pathway (e.g. hsa04150, mmu04130).

    save_to
        The directory where the KGML is saved.
        Default: "." (i.e. current directory)

    encoding
        Encoding.
        Default: "utf-8"

    kwargs
        Parameters of requests_get().
    '''
    response = requests_get(f'https://rest.kegg.jp/get/{pathway}/kgml', **kwargs)
    contents = response.text

    with open (save_to + f'/{pathway}.xml', 'w', encoding=encoding) as save_file:
        save_file.write(contents)

    return 1



class KEGGIDConverter:
    '''
    This class converts KEGG IDs to gene symbols and vice versa.

    Arguments
    ----------
    organism (required)
        The KEGG organism (e.g. hsa, mmu).

    kwargs
        Parameters of requests_get().
    '''
    def __init__(self, organism, **kwargs):
        self.organism = organism

        # The conversion of KEGG ID and Entrez ID.
        conversions = requests_get(f'https://rest.kegg.jp/conv/{organism}/ncbi-geneid', **kwargs)

        self.KEGGvsENTREZ_data = []
        for row in conversions.text.strip().split('\n'):
            ENTREZID, KEGGID = row.split('\t')
            self.KEGGvsENTREZ_data.append({'KEGGID':   KEGGID,
                                           'ENTREZID': re.search(r'[0-9]+', ENTREZID).group()})
        self.KEGGvsENTREZ_data = pd.DataFrame(self.KEGGvsENTREZ_data)

        # Taxonomy ID.
        html = requests_get('https://www.genome.jp/kegg-bin/show_organism?org=' + organism)
        soup = BeautifulSoup(html.text, 'html.parser')

        taxon = soup.find('a', href=re.compile(r'.+Taxonomy.+'))
        self.taxonomy_id = re.search(r'(?<=id=)[0-9]+$', taxon.get('href')).group()

        # For the conversion of Entrez ID and Symbol.
        self.mg = MyGeneInfo()


    def KEGGIDtoENTREZID(self, KEGGIDs):
        '''
        This method converts given KEGG IDs to Entrez IDs.

        Arguments
        ----------
        KEGGIDs (required)
            The KEGG IDs.
        '''
        result = pd.DataFrame(KEGGIDs, columns=['KEGGID'])
        result = pd.merge(result, self.KEGGvsENTREZ_data, on='KEGGID', how='left').fillna('')
        result = result.drop_duplicates().reset_index(drop=True)
        return result


    def ENTREZIDtoKEGGID(self, ENTREZIDs):
        '''
        This method converts given Entrez IDs to KEGG IDs.

        Arguments
        ----------
        ENTREZIDs (required)
            The Entrez IDs.
        '''
        result = pd.DataFrame(ENTREZIDs, columns=['ENTREZID'])
        result = pd.merge(result, self.KEGGvsENTREZ_data, on='ENTREZID', how='left').fillna('')
        result = result.drop_duplicates().reset_index(drop=True)
        return result


    def ENTREZIDtoSYMBOL(self, ENTREZIDs):
        '''
        This method converts given Entrez IDs to gene symbols.

        Arguments
        ----------
        ENTREZIDs (required)
            The Entrez IDs.
        '''
        result = self.mg.querymany(ENTREZIDs,
                                   scopes  = 'entrezgene',
                                   species = self.taxonomy_id,
                                   verbose = False,
                                   field   = 'symbol,entrezgene',
                                   as_dataframe = True)
        result = result.reset_index()[['symbol','query']].fillna('')
        result = result.drop_duplicates().reset_index(drop=True)
        result = result.astype(str)
        result = result.rename(columns={'symbol': 'SYMBOL',
                                        'query':  'ENTREZID'})
        result = result.groupby('ENTREZID')['SYMBOL'].apply(self._remove_blanks).reset_index()
        result = result.explode('SYMBOL').fillna('')
        return result


    def SYMBOLtoENTREZID(self, SYMBOLs):
        '''
        This method converts given gene symbols to Entrez IDs.

        Arguments
        ----------
        SYMBOLs [list] (required)
            The gene symbols.
        '''
        result = self.mg.querymany(SYMBOLs,
                                   scopes  = 'symbol',
                                   species = self.taxonomy_id,
                                   verbose = False,
                                   field   = 'symbol,entrezgene',
                                   as_dataframe = True)
        result = result.reset_index()[['entrezgene','query']].fillna('')
        result = result.drop_duplicates().reset_index(drop=True)
        result = result.astype(str)
        result = result.rename(columns={'query':      'SYMBOL',
                                        'entrezgene': 'ENTREZID'})
        result = result.groupby('SYMBOL')['ENTREZID'].apply(self._remove_blanks).reset_index()
        result = result.explode('ENTREZID').fillna('')
        return result


    def KEGGIDtoSYMBOL(self, KEGGIDs):
        '''
        This method converts given KEGG IDs to gene symbols.

        Arguments
        ----------
        KEGGIDs (required)
            The KEGG IDs.
        '''
        result_KtoE = self.KEGGIDtoENTREZID(KEGGIDs)
        result_EtoS = self.ENTREZIDtoSYMBOL(result_KtoE.query('ENTREZID!=""')['ENTREZID'].unique())

        result = pd.merge(result_KtoE, result_EtoS, on='ENTREZID', how='left').fillna('')
        result = result.drop_duplicates().reset_index(drop=True)
        return result


    def SYMBOLtoKEGGID(self, SYMBOLs):
        '''
        This method converts given gene symbols to KEGG IDs.

        Arguments
        ----------
        SYMBOLs (required)
            The gene symbols.
        '''
        result_StoE = self.SYMBOLtoENTREZID(SYMBOLs)
        result_EtoK = self.ENTREZIDtoKEGGID(result_StoE.query('ENTREZID!=""')['ENTREZID'].unique())

        result = pd.merge(result_StoE, result_EtoK, on='ENTREZID', how='left').fillna('')
        result = result.drop_duplicates().reset_index(drop=True)
        return result


    def _remove_blanks(self, x):
        return [ i for i in x if i != '' ]