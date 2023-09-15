## Perform the KEGG-mining for F0-3.
##
##
organism = 'mmu'

# Download KGMLs.
#
#
import os
import time
from tqdm import tqdm
from KEGGmining_utils import get_pathways, download_KGML

kgml_dir = '../data/KGMLs'

if not os.path.exists(kgml_dir):
    os.mkdir(kgml_dir)

pathways = get_pathways(organism)

for pathway in tqdm(pathways, desc='Downloading'):
    download_KGML(pathway, save_to=kgml_dir)
    time.sleep(0.1)



# Prepare data of DEgenes.
#
#
import pandas as pd

data_files = {'F0-3 Adult': ('../data/DEseq_results/Adult_5wksND_9organs/'
                             'DESeqs_8wksPrpsALLKOF0-3deltaX_8wksWT.xlsx'),
              'F0-3 Young': ('../data/DEseq_results/Young_1wND_24organs/'
                             'DESeqs_4wksPrpsALLKO1wND_4wksWT.xlsx')}

thresholds = {'adj.p': 0.05, 'log2FC': 0.5}


DEgene_data = []

for condition, file in data_files.items():
    DEseq_data = pd.read_excel(file)
    prefix = DEseq_data.columns[1].split('_')[0]
    organs = sorted(set( col.split('_')[1] for col in DEseq_data.columns[1:] ))

    for organ in organs:
        columns = {'adj.p':  f'{prefix}_{organ}_padj',
                   'log2FC': f'{prefix}_{organ}_log2FoldChange'}

        # adjusted p-value
        add_DEgene_data = DEseq_data[DEseq_data[columns['adj.p']] < thresholds['adj.p']]

        # log2FoldChange
        add_DEgene_data = add_DEgene_data[(add_DEgene_data[columns['log2FC']] >  thresholds['log2FC']) |
                                          (add_DEgene_data[columns['log2FC']] < -thresholds['log2FC'])]

        if len(add_DEgene_data) > 0:
            add_DEgene_data = add_DEgene_data[['gene_symbols', columns['adj.p'], columns['log2FC']]]
            add_DEgene_data = add_DEgene_data.rename(columns={'gene_symbols':    'DEgene',
                                                              columns['adj.p']:  'adj.p',
                                                              columns['log2FC']: 'log2FC'})
            add_DEgene_data['condition'] = condition
            add_DEgene_data['organ']     = organ
        else:
            add_DEgene_data = pd.DataFrame([{'condition': condition,
                                             'organ':     organ,
                                             'DEgene':    '',
                                             'adj.p':     None,
                                             'log2FC':    None}])

        DEgene_data.append(add_DEgene_data)

DEgene_data = pd.concat(DEgene_data).drop_duplicates()


organs_a = DEgene_data.query('condition=="F0-3 Adult"')['organ'].unique()
organs_y = DEgene_data.query('condition=="F0-3 Young"')['organ'].unique()
used_organs = [ organ for organ in organs_a if organ in organs_y ]

DEgene_data = DEgene_data[DEgene_data['organ'].isin(used_organs)]


# Convert gene symbols to KEGG IDs.
#
#
from KEGGmining_utils import KEGGIDConverter

converter = KEGGIDConverter(organism)

conv_data = converter.SYMBOLtoKEGGID(DEgene_data.query('DEgene!=""')['DEgene'].unique())
conv_data = conv_data.drop(columns=['ENTREZID'])
conv_data = conv_data.rename(columns={'SYMBOL': 'DEgene',
                                      'KEGGID': 'DEgene.KEGGID'})

DEgene_data = pd.merge(DEgene_data, conv_data, on='DEgene', how='left')



# Prepare data of ligand-receptor pairs.
#
#
cellchat_file = '../data/CellChatDB_for_mmu.xlsx'
cellchat_data = pd.read_excel(cellchat_file).fillna('')

LRdata = []
for row in cellchat_data.itertuples():
    ligand1   = row.ligand1.strip()
    ligand2   = row.ligand2.strip()
    receptor1 = row.receptor1.strip()
    receptor2 = row.receptor2.strip()
    if ligand1 != '':
        if receptor1 != '':
            LRdata.append({'ligand':   ligand1,
                           'receptor': receptor1})
        if receptor2 != '':
            LRdata.append({'ligand':   ligand1,
                           'receptor': receptor2})
    if ligand2 != '':
        if receptor1 != '':
            LRdata.append({'ligand':   ligand2,
                           'receptor': receptor1})
        if receptor2 != '':
            LRdata.append({'ligand':   ligand2,
                           'receptor': receptor2})
LRdata = pd.DataFrame(LRdata)

conv_data = converter.SYMBOLtoKEGGID(LRdata['ligand'].unique())
conv_data = conv_data.drop(columns=['ENTREZID'])
conv_data = conv_data.rename(columns={'SYMBOL': 'ligand',
                                      'KEGGID': 'ligand.KEGGID'})

LRdata = pd.merge(LRdata, conv_data, on='ligand', how='left').fillna('')

conv_data = converter.SYMBOLtoKEGGID(LRdata['receptor'].unique())
conv_data = conv_data.drop(columns=['ENTREZID'])
conv_data = conv_data.rename(columns={'SYMBOL': 'receptor',
                                      'KEGGID': 'receptor.KEGGID'})

LRdata = pd.merge(LRdata, conv_data, on='receptor', how='left').fillna('')


import re
from itertools import product
from KGMLGraph_utils import KGMLReader

LRpathways = [ organism + number for number in ['04060','04080','04512'] ]

LRdata_add = []
named = re.compile(r' |\+')
gene  = re.compile(f'^{organism}:')

for pathway in LRpathways:
    kgml = KGMLReader(f'{kgml_dir}/{pathway}.xml')

    for entry1, entry2 in kgml.relation_data[['entry1','entry2']].drop_duplicates().values:
        inputs  = named.split(kgml.entry_data.loc[entry1,'name'])
        outputs = named.split(kgml.entry_data.loc[entry2,'name'])
        genes_in  = [ entry for entry in inputs  if gene.match(entry) ]
        genes_out = [ entry for entry in outputs if gene.match(entry) ]
        LRdata_add += [ {'ligand.KEGGID':   pair[0],
                         'receptor.KEGGID': pair[1]} for pair in product(genes_in, genes_out)]
LRdata_add = pd.DataFrame(LRdata_add)

conv_data = converter.KEGGIDtoSYMBOL(LRdata_add['ligand.KEGGID'].unique())
conv_data = conv_data.drop(columns=['ENTREZID'])
conv_data = conv_data.rename(columns={'SYMBOL': 'ligand',
                                      'KEGGID': 'ligand.KEGGID'})

LRdata_add = pd.merge(LRdata_add, conv_data, on='ligand.KEGGID', how='left').fillna('')

conv_data = converter.KEGGIDtoSYMBOL(LRdata['receptor.KEGGID'].unique())
conv_data = conv_data.drop(columns=['ENTREZID'])
conv_data = conv_data.rename(columns={'SYMBOL': 'receptor',
                                      'KEGGID': 'receptor.KEGGID'})

LRdata_add = pd.merge(LRdata_add, conv_data, on='receptor.KEGGID', how='left').fillna('')

LRdata = pd.concat([LRdata, LRdata_add]).drop_duplicates()

#LRdata.to_excel('../data/ligand-receptor_pairs.xlsx', index=False)



# KEGG-mining
#
#
from KEGGmining import KEGGmining

used_LRdata = LRdata[['ligand.KEGGID','receptor.KEGGID']].drop_duplicates()
used_LRdata = used_LRdata[(used_LRdata['ligand.KEGGID']  !='') &
                          (used_LRdata['receptor.KEGGID']!='')]
used_LRdata = used_LRdata.rename(columns={'ligand.KEGGID':   'ligand',
                                          'receptor.KEGGID': 'receptor'})

kgm = KEGGmining(organism, kgml_dir, LRdata=used_LRdata)

results = []
for condition in data_files.keys():
    DEgenes = DEgene_data[DEgene_data['condition']==condition]
    DEgenes = DEgenes[DEgenes['DEgene.KEGGID']!='']['DEgene.KEGGID'].unique()

    result = kgm(DEgenes)
    result['condition'] = condition

    results.append(result)
results = pd.concat(results)

#results.to_excel('../data/KEGG-mining_results.xlsx', index=False)



# Get the ligand-receptor pairs for each KEGG pathway.
#
#
from KGMLGraph_utils import KGMLGeneGraph

LRpair_data = []
for pathway in tqdm(kgm.pathways):
    G = KGMLGeneGraph(f'{kgml_dir}/{pathway}.xml', LRdata=kgm.LRdata)
    LRpair_data.append(G.LRdata)
LRpair_data = pd.concat(LRpair_data).drop_duplicates()

LRpair_data['pathway'] = LRpair_data['receptor_as_node'].apply(lambda x: x.split('__')[1])
LRpair_data['LRpair']  = LRpair_data['ligand'] + '__' + \
                         LRpair_data['receptor_as_node'].apply(lambda x: x.split('__')[0])

#LRpair_data.to_excel('../data/ligand-receptor_pairs_on_KEGG.xlsx', index=False)

ligands   = LRpair_data.query('ligand!=""')['ligand'].unique()
receptors = LRpair_data.query('receptor_as_node!=""')['receptor_as_node'].unique()



# Generate edges from DEgenes to ligands.
#
#
DLedge_data = []
for condition in data_files.keys():
    DEgenes = DEgene_data[DEgene_data['condition']==condition]
    DEgenes = DEgenes[DEgenes['DEgene.KEGGID']!='']['DEgene.KEGGID'].unique()
    result  = results[results['condition']==condition]

    DLedges = result[(result['gene.KEGGID'].isin(DEgenes)) &
                     (result['connected.KEGGID'].isin(ligands)) &
                     (result['direction']=='downstream')].drop(columns=['direction'])
    DLedges = DLedges.rename(columns={'pathway':          'pathway.DtoL',
                                      'gene.KEGGID':      'from.KEGGID',
                                      'connected.KEGGID': 'ligand.KEGGID'})

    ligands_in_DEgenes = [ gene for gene in ligands if gene in DEgenes ]

    add_DLedges = []
    for pathway in kgm.pathways:
        kgml = KGMLReader(f'{kgml_dir}/{pathway}.xml')
        components = kgml.all_components()
        add_DLedges += [ {'pathway.DtoL':  pathway,
                          'from.KEGGID':   gene,
                          'ligand.KEGGID': gene}
                         for gene in ligands_in_DEgenes if gene in components ]
    add_DLedges = pd.DataFrame(add_DLedges)

    DLedges = pd.concat([DLedges, add_DLedges]).drop_duplicates()

    use_DEgene_data = DEgene_data[DEgene_data['condition']==condition]
    use_DEgene_data = use_DEgene_data[['organ','DEgene','DEgene.KEGGID']].drop_duplicates()

    DLedges = pd.merge(DLedges, use_DEgene_data,
                       left_on  = 'from.KEGGID',
                       right_on = 'DEgene.KEGGID').drop(columns=['DEgene.KEGGID'])
    DLedges = DLedges.rename(columns={'organ':  'from.organ',
                                      'DEgene': 'from.DEgene'})

    DLedges = pd.merge(DLedges, LRdata[['ligand','ligand.KEGGID']].drop_duplicates(),
                       on = 'ligand.KEGGID').drop_duplicates()

    DLedges = DLedges.groupby(
                  ['from.organ','from.DEgene','from.KEGGID',
                   'ligand','ligand.KEGGID'])['pathway.DtoL'].apply(lambda x: ' '.join(sorted(x))).reset_index()

    DLedges['condition'] = condition
    DLedge_data.append(DLedges)

DLedge_data = pd.concat(DLedge_data)



# Generate edges from ligand-receptor pairs to DEgenes.
#
#
LRedge_data = []
for row in LRpair_data.itertuples():
    LRedge_data.append({'pathway':          row.pathway,
                        'ligand.KEGGID':    row.ligand,
                        '_receptor.KEGGID': row.LRpair,
                        'receptor.KEGGID':  row.receptor})
LRedge_data = pd.DataFrame(LRedge_data).drop_duplicates()

LRedge_data = pd.merge(LRedge_data,
                       LRdata[['ligand','ligand.KEGGID']].drop_duplicates(),
                       on = 'ligand.KEGGID')

LRedge_data = pd.merge(LRedge_data,
                       LRdata[['receptor','receptor.KEGGID']].drop_duplicates(),
                       on = 'receptor.KEGGID').drop(columns=['receptor.KEGGID'])

LRedge_data = LRedge_data.rename(columns={'_receptor.KEGGID': 'receptor.KEGGID'})


upstream_data = []
for condition in data_files.keys():
    DEgenes = DEgene_data[DEgene_data['condition']==condition]
    DEgenes = DEgenes[DEgenes['DEgene.KEGGID']!='']['DEgene.KEGGID'].unique()
    result  = results[results['condition']==condition]

    RDedge_data = result[(result['gene.KEGGID'].isin(DEgenes)) &
                         (result['connected.KEGGID'].isin(receptors)) &
                         (result['direction']=='upstream')].drop(columns=['direction'])
    RDedge_data = RDedge_data.rename(columns={'gene.KEGGID':      'DEgene.KEGGID',
                                              'connected.KEGGID': 'receptor_as_node'})

    RDedge_data = pd.merge(RDedge_data, LRpair_data, on=['pathway','receptor_as_node'])
    RDedge_data = RDedge_data[['pathway','LRpair','DEgene.KEGGID','receptor']].drop_duplicates()
    RDedge_data = RDedge_data.rename(columns={'LRpair':   '_receptor.KEGGID',
                                              'receptor': 'receptor.KEGGID'})

    add_RDedge_data = LRpair_data[LRpair_data['receptor'].isin(DEgenes)]
    add_RDedge_data = add_RDedge_data[['pathway','LRpair','receptor']].drop_duplicates()
    add_RDedge_data = add_RDedge_data.rename(columns={'LRpair':   '_receptor.KEGGID',
                                                      'receptor': 'receptor.KEGGID'})
    add_RDedge_data['DEgene.KEGGID'] = add_RDedge_data['receptor.KEGGID'].copy()

    RDedge_data = pd.concat([RDedge_data, add_RDedge_data]).drop_duplicates()

    use_DEgene_data = DEgene_data[DEgene_data['condition']==condition]
    use_DEgene_data = use_DEgene_data[['organ','DEgene','DEgene.KEGGID']].drop_duplicates()

    RDedge_data = pd.merge(RDedge_data, use_DEgene_data, on='DEgene.KEGGID')

    RDedge_data = pd.merge(RDedge_data,
                           LRdata[['receptor','receptor.KEGGID']].drop_duplicates(),
                           on = 'receptor.KEGGID').drop(columns=['receptor.KEGGID'])
    RDedge_data = RDedge_data.rename(columns={'_receptor.KEGGID': 'receptor.KEGGID'})


    upstream_edges = pd.merge(LRedge_data, RDedge_data,
                                  on = ['pathway','receptor','receptor.KEGGID']).drop_duplicates()
    upstream_edges['receptor.KEGGID'] = upstream_edges['receptor.KEGGID'].apply(lambda x: x.split('__')[1])
    upstream_edges = upstream_edges.drop_duplicates()

    upstream_edges = upstream_edges.groupby(
                             ['ligand','ligand.KEGGID',
                              'receptor','receptor.KEGGID',
                              'organ','DEgene','DEgene.KEGGID'])['pathway'].apply(lambda x: ' '.join(sorted(x))).reset_index()

    upstream_edges = upstream_edges.reindex(
                             columns = ['pathway','ligand','ligand.KEGGID',
                                        'receptor','receptor.KEGGID',
                                        'organ','DEgene','DEgene.KEGGID'])

    upstream_edges['condition'] = condition
    upstream_data.append(upstream_edges)

upstream_data = pd.concat(upstream_data)



# Generate edges between DEgenes via ligand-receptor pairs.
#
#
edge_data = []
for condition in data_files.keys():
    LRDedge_data = upstream_data[upstream_data['condition']==condition]
    LRDedge_data = LRDedge_data.rename(columns={'pathway':       'pathway.LRtoD',
                                                'organ':         'to.organ',
                                                'DEgene':        'to.DEgene',
                                                'DEgene.KEGGID': 'to.KEGGID'})

    DLedges = DLedge_data[DLedge_data['condition']==condition]

    edges = pd.merge(DLedges, LRDedge_data, on=['ligand','ligand.KEGGID'])
    edges = edges.reindex(columns=['from.organ','from.DEgene','from.KEGGID','pathway.DtoL',
                                   'ligand','ligand.KEGGID','receptor','receptor.KEGGID',
                                   'pathway.LRtoD','to.organ','to.DEgene','to.KEGGID'])
    edges = edges.drop_duplicates()

    edges['condition'] = condition
    edge_data.append(edges)

edge_data = pd.concat(edge_data)


# Save data.
#
#
DEgene_data = DEgene_data[['organ','DEgene','DEgene.KEGGID']]
DEgene_data = DEgene_data.drop_duplicates()
DEgene_data.to_excel('../data/DEgene_data.xlsx', index=False)

edge_data = edge_data.drop(columns=['condition'])
edge_data = edge_data.drop_duplicates()
edge_data.to_excel('../data/edge_data.xlsx', index=False)

upstream_data = upstream_data.drop(columns=['condition'])
upstream_data = upstream_data.drop_duplicates()
upstream_data.to_excel('../data/upstream_data.xlsx', index=False)

DLedge_data = DLedge_data.drop(columns=['condition'])
DLedge_data = DLedge_data.drop_duplicates()
DLedge_data.to_excel('../data/DEgene-ligand_data.xlsx', index=False)