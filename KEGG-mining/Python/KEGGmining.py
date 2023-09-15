## Perform the KEGG-mining.
##
##
import pandas as pd
from tqdm import tqdm

from KEGGmining_utils import get_pathways
from KGMLGraph_utils import KGMLGeneGraph


class KEGGmining:
    '''
    This class performs the KEGG-miming.

    Arguments
    ----------
    organism (required)
        The KEGG organism (e.g. hsa, mmu).

    kgml_dir (required)
        The directory where KGMLs are saved.

    LRdata
        Data of ligand-receptor pairs.
        2 columns ligand and receptor are required.
        Default: None

    kwargs
        Parameters of get_pathways().
    '''
    def __init__(self, organism, kgml_dir, LRdata=None, **kwargs):
        self.organism = organism
        self.kgml_dir = kgml_dir
        self.pathways = get_pathways(organism, **kwargs)

        if LRdata is not None:
            self.LRdata = LRdata.copy()
        else:
            self.LRdata = pd.DataFrame(columns=['ligand','receptor'])


    def __call__(self, KEGGIDs, disable=False):
        '''
        This method performs the KEGG-mining for given KEGG IDs.

        Arguments
        ----------
        KEGGIDs (required)
            The KEGG IDs.

        disable
            Do not show the progress bar.
            Deafult: False (i.e. Bar shown)
        '''
        results = []

        for pathway in tqdm(self.pathways, desc='Searching', disable=disable):
            G = KGMLGeneGraph(f'{self.kgml_dir}/{pathway}.xml', LRdata=self.LRdata)
            result = self._find_connected_genes(G, KEGGIDs, G.LRdata)
            results.append(result)

        if len(results) > 0:
            results = pd.concat(results)
            return results.drop_duplicates().reset_index(drop=True)
        else:
            return pd.DataFrame(columns=['pathway','gene.KEGGID','connected.KEGGID','direction'])


    def _find_connected_genes(self, G, KEGGIDs, LRdata):
        '''
        This method finds the genes connected to given KEGG IDs.
        '''
        result = []

        for KEGGID in KEGGIDs:
            if G.has_node(KEGGID):
                if G.degree(KEGGID) > 0:
                    # Dowmstream.
                    downstreams = self._find_dowmstreams(G, KEGGID, LRdata)

                    for genes in downstreams:
                        result += [ {'pathway':          G.kgml.KEGGID,
                                     'gene.KEGGID':      KEGGID,
                                     'connected.KEGGID': gene,
                                     'direction':        'downstream'} for gene in genes ]

                    # Upstream.
                    upstreams = self._find_upstreams(G, KEGGID, LRdata)

                    for genes in upstreams:
                        result += [ {'pathway':          G.kgml.KEGGID,
                                     'gene.KEGGID':      KEGGID,
                                     'connected.KEGGID': gene,
                                     'direction':        'upstream'} for gene in genes ]
                else:
                    result.append({'pathway':          G.kgml.KEGGID,
                                   'gene.KEGGID':      KEGGID,
                                   'connected.KEGGID': '',
                                   'direction':        ''})

        if len(result) > 0:
            return pd.DataFrame(result).drop_duplicates().reset_index(drop=True)
        else:
            return pd.DataFrame(columns=['pathway','gene.KEGGID','connected.KEGGID','direction'])


    def _find_dowmstreams(self, G, KEGGID, LRdata):
        '''
        This method returns genes connected downstream to a given KEGG ID.
        '''
        LRmapper = LRdata.groupby('ligand')['receptor_as_node'].apply(set).reset_index()
        LRmapper = { row.ligand: row.receptor_as_node for row in LRmapper.itertuples() }

        downstreams = {KEGGID}
        hangouts    = {KEGGID}

        while True:
            genes_add = set()

            for gene in downstreams:
                receptors_rm  = LRmapper[gene] if gene in LRmapper else set()
                genes_add    |= set(G.successors(gene)) - receptors_rm
            genes_add -= hangouts

            if len(genes_add) == 0:
                break
            else:
                yield genes_add

                downstreams  = genes_add
                hangouts    |= genes_add


    def _find_upstreams(self, G, KEGGID, LRdata):
        '''
        This method returns genes connected upstream to a given KEGG ID.
        '''
        RLmapper = LRdata.groupby('receptor_as_node')['ligand'].apply(set).reset_index()
        RLmapper = { row.receptor_as_node: row.ligand for row in RLmapper.itertuples() }

        upstreams = {KEGGID}
        hangouts  = {KEGGID}

        while True:
            genes_add = set()

            for gene in upstreams:
                ligands_rm  = RLmapper[gene] if gene in RLmapper else set()
                genes_add  |= set(G.predecessors(gene)) - ligands_rm
            genes_add -= hangouts

            if len(genes_add) == 0:
                break
            else:
                yield genes_add

                upstreams  = genes_add
                hangouts  |= genes_add