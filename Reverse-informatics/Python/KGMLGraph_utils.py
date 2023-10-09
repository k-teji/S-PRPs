## Collection of functions for graph construction.
##
##
import networkx as nx
import pandas as pd
import re
import xml.etree.ElementTree as ET
from itertools import product



class KGMLReader:
    '''
    This class stores <entry>, <relation> and <reaction> of
    a given KGML as dataframe.

    Arguments
    ----------
    file (required)
        The filepath of the KGML.
    '''
    def __init__(self, file):
        tree = ET.parse(file)
        root = tree.getroot()

        self.title    = root.get('title')  # Pathway title.
        self.organism = root.get('org')    # KEGG organism.
        self.number   = root.get('number') # Numeric part.
        self.KEGGID   = self.organism + self.number

        # <entry>

        entry_data = []
        groups     = []
        for entry in root.findall('entry'):
            entry_id   = int(entry.get('id'))
            entry_type = entry.get('type')

            if entry_type != 'group':
                entry_components = []
                entry_name       = entry.get('name')
            else:
                entry_components = sorted(self._get_components(root, entry))

                component_names = []
                for component_id in entry_components:
                    component = root.find(f'.//entry[@id="{component_id}"]')
                    component_names.append(component.get('name'))

                    groups.append({'entry_id': component_id,
                                   'groups':   entry_id})

                entry_name = '+'.join(sorted(component_names))

            entry_reaction = entry.get('reaction')
            if entry_reaction is None:
                entry_reaction = ''

            entry_graphics = entry.find('graphics')
            entry_shape    = entry_graphics.get('type')
            if entry_shape != 'line':
                entry_coords = (int(entry_graphics.get('x')),
                                int(entry_graphics.get('y')))
            else:
                entry_coords = tuple( int(coord) for coord in entry_graphics.get('coords').split(',') )
                entry_coords = [ entry_coords[:2], entry_coords[2:] ]

            entry_data.append({'entry_id':   entry_id,
                               'type':       entry_type,
                               'name':       entry_name,
                               'components': entry_components,
                               'reaction':   entry_reaction,
                               'shape':      entry_shape,
                               'coords':     entry_coords})

        entry_columns = ['entry_id','type','name','components','groups','reaction','shape','coords']

        if len(entry_data) > 0:
            self.entry_data = pd.DataFrame(entry_data)

            if len(groups) > 0:
                group_data = pd.DataFrame(groups)
                group_data = group_data.groupby('entry_id')['groups'].apply(sorted).reset_index()
                self.entry_data = pd.merge(self.entry_data,
                                           group_data,
                                           on  = 'entry_id',
                                           how = 'left')
            else:
                self.entry_data['groups'] = ''

            self.entry_data['groups'] = self.entry_data['groups'].fillna('').apply(list)
            self.entry_data = self.entry_data.reindex(columns=entry_columns)
        else:
            self.entry_data = pd.DataFrame(columns=entry_columns)

        self.entry_data.set_index('entry_id', inplace=True)
        self.n_entries = len(self.entry_data.index.unique())

        # <relation>

        relation_data = []
        for rel_id, relation in enumerate(root.findall('relation')):
            rel_entry1 = int(relation.get('entry1'))
            rel_entry2 = int(relation.get('entry2'))
            rel_type   = relation.get('type')

            rel_subtypes = relation.findall('subtype')
            if len(rel_subtypes) > 0:
                for subtype in rel_subtypes:
                    sub_name  = subtype.get('name')
                    sub_value = subtype.get('value')
                    if sub_name in ['compound','hidden compound']:
                        sub_value = int(sub_value)

                    relation_data.append({'relation_id':   rel_id,
                                          'type':          rel_type,
                                          'entry1':        rel_entry1,
                                          'entry2':        rel_entry2,
                                          'subtype':       sub_name,
                                          'subtype_value': sub_value})
            else:
                relation_data.append({'relation_id':   rel_id,
                                      'type':          rel_type,
                                      'entry1':        rel_entry1,
                                      'entry2':        rel_entry2,
                                      'subtype':       '',
                                      'subtype_value': ''})

        relation_columns = ['relation_id','type','entry1','entry2','subtype','subtype_value']

        if len(relation_data) > 0:
            self.relation_data = pd.DataFrame(relation_data)
        else:
            self.relation_data = pd.DataFrame(columns=relation_columns)

        self.relation_data.set_index('relation_id', inplace=True)
        self.n_relations = len(self.relation_data.index.unique())

        # <reaction>

        reaction_data = []

        for reaction in root.findall('reaction'):
            react_id   = int(reaction.get('id'))
            react_name = reaction.get('name')
            react_type = reaction.get('type')
            react_substrates = [ int(subs.get('id')) for subs in reaction.findall('substrate') ]
            react_products   = [ int(prod.get('id')) for prod in reaction.findall('product')   ]

            reaction_data.append({'reaction_id': react_id,
                                  'name':        react_name,
                                  'type':        react_type,
                                  'substrates':  react_substrates,
                                  'products':    react_products})

        react_columns = ['reaction_id','name','type','substrates','products']

        if len(reaction_data) > 0:
            self.reaction_data = pd.DataFrame(reaction_data)
        else:
            self.reaction_data = pd.DataFrame(columns=react_columns)

        self.reaction_data.set_index('reaction_id', inplace=True)
        self.n_reactions = len(self.reaction_data.index.unique())


    def _get_components(self, root, entry):
        components = set()
        for component in entry.findall('component'):
            component_id = int(component.get('id'))
            component    = root.find(f'.//entry[@id="{component_id}"]')

            if component.get('type') != 'group':
                components.add(component_id)
            else:
                components |= self._get_components(root, component)

        return components


    def __repr__(self):
        return f'<KGML id="{self.KEGGID}">'

    def __getitem__(self, entry_id):
        return self.entry_data.loc[entry_id,:]


    def all_components(self):
        '''
        This method returns all components contained in the KGML.
        '''
        return sorted( component for name      in self.entry_data['name'].unique()
                                 for component in re.split(r' |\+', name) )


    def find_entry_containing(self, KEGGID):
        '''
        This method finds the <entry> containing a given KEGG ID.

        Arguments
        ----------
        KEGGID (required)
            The KEGG ID.
        '''
        return [ entry.Index for entry in self.entry_data.itertuples()
                             if KEGGID in re.split(' |\+', entry.name)]



class KGMLGraph(nx.DiGraph):
    '''
    This class is a directed graph constructed from a given KGML.

    The nodes are the attributes id of <entry>.
    The edges are based on <reaction> and <relation>.

    Arguments
    ----------
    file (required)
        The filepath of the KGML.
    '''
    def __init__(self, file):
        super().__init__()

        self.kgml = KGMLReader(file)

        self.add_nodes_from(self.kgml.entry_data.index)

        # entry ID --> group IDs containing it
        mapper = dict()
        for entry in self.kgml.entry_data.itertuples():
            if len(entry.groups) > 0:
                mapper[entry.Index] = entry.groups
            else:
                mapper[entry.Index] = [entry.Index]

        # Edges by <reaction>.
        for reaction in self.kgml.reaction_data.itertuples():
            enzymes    = mapper[reaction.Index]
            substrates = { group_id for subs     in reaction.substrates
                                    for group_id in mapper[subs] }
            products   = { group_id for prod     in reaction.products
                                    for group_id in mapper[prod] }

            self.add_edges_from( (subs, enz) for subs, enz in product(substrates, enzymes) )
            self.add_edges_from( (enz, prod) for enz, prod in product(enzymes, products)   )

            if reaction.type == 'reversible':
                self.add_edges_from( (enz, subs) for subs, enz in product(substrates, enzymes) )
                self.add_edges_from( (prod, enz) for enz, prod in product(enzymes, products)   )

        # Edges by <relation>.
        for entry1, entry2 in self.kgml.relation_data[['entry1','entry2']].drop_duplicates().values:
            groups1 = mapper[entry1]
            groups2 = mapper[entry2]

            rel_data = self.kgml.relation_data[(self.kgml.relation_data['entry1']==entry1) &
                                               (self.kgml.relation_data['entry2']==entry2)]

            compounds   = rel_data[rel_data['subtype'].isin(['compound','hidden compound'])]
            n_compounds = len(compounds)

            if n_compounds > 0:
                compounds = compounds['subtype_value'].unique()
                compounds = { group_id for cpd      in compounds
                                       for group_id in mapper[cpd] }

                if len(rel_data) == n_compounds:
                    edges = []

                    for grp, cpd in product(groups1, compounds):
                        if not (self.has_edge(grp, cpd) or self.has_edge(cpd, grp)):
                            edges += [(grp, cpd), (cpd, grp)]
                    for grp, cpd in product(groups2, compounds):
                        if not (self.has_edge(grp, cpd) or self.has_edge(cpd, grp)):
                            edges += [(grp, cpd), (cpd, grp)]

                    self.add_edges_from(edges)
                else:
                    n_bidirect = len(rel_data[rel_data['subtype']=='binding/association'])

                    if n_bidirect > 0:
                        for grp1, cpd, grp2 in product(groups1, compounds, groups2):
                            self.add_edges_from([(grp1, cpd), (cpd, grp2),
                                                 (grp2, cpd), (cpd, grp1)])
                    else:
                        for grp1, cpd, grp2 in product(groups1, compounds, groups2):
                            self.add_edges_from([(grp1, cpd), (cpd, grp2)])
            else:
                n_bidirect = len(rel_data[rel_data['subtype'].isin(['','binding/association'])])

                if n_bidirect > 0:
                    for grp1, grp2 in product(groups1, groups2):
                        self.add_edges_from([(grp1, grp2), (grp2, grp1)])
                else:
                    self.add_edges_from( (grp1, grp2) for grp1, grp2 in product(groups1, groups2) )

        # The indirect effects.
        indirect_effects = self.kgml.relation_data[self.kgml.relation_data['subtype']=='indirect effect']
        indirect_effects = indirect_effects[['entry1','entry2']].drop_duplicates()
        self.indirect_effects = [ (row.entry1, row.entry2) for row in indirect_effects.itertuples() ]


    @property
    def title(self):
        return self.kgml.title

    @property
    def organism(self):
        return self.kgml.organism

    @property
    def number(self):
        return self.kgml.number

    @property
    def KEGGID(self):
        return self.kgml.KEGGID

    @property
    def entry_data(self):
        return self.kgml.entry_data

    @property
    def n_entries(self):
        return self.kgml.n_entries

    @property
    def relation_data(self):
        return self.kgml.relation_data

    @property
    def n_relations(self):
        return self.kgml.n_relations

    @property
    def reaction_data(self):
        return self.kgml.reaction_data

    @property
    def n_reactions(self):
        return self.kgml.n_reactions

    def __repr__(self):
        return f'<KGMLGraph id="{self.kgml.KEGGID}">'

    def __getitem__(self, entry_id):
        return self.kgml.__getitem__(entry_id)



class KGMLGeneGraph(KGMLGraph):
    '''
    This class is a directed graph constructed from a given KGML.
    The nodes are the names of genes.

    Arguments
    ----------
    file (required)
        The filepath of the KGML.

    LRdata
        Data of ligand-receptor pairs.
        2 columns ligand and receptor are required.
        Default: None.
    '''
    def __init__(self, file, LRdata=None):
        super().__init__(file)

        self.LRdata = []

        if LRdata is None:
            mapper = { entry_id: self._extract_genes(entry_id) for entry_id in self.nodes() }
        else:
            mapper    = dict()
            receptors = set(LRdata['receptor'])

            for entry_id in self.nodes():
                entry_genes     = set(self._extract_genes(entry_id))
                entry_receptors = (entry_genes & receptors)

                if len(entry_receptors) == 0:
                    mapper[entry_id] = list(entry_genes)
                else:
                    entry_LRdata  = LRdata[LRdata['receptor'].isin(entry_receptors)]
                    entry_ligands = set(entry_LRdata['ligand'])
                    input_ligands = dict()
                    is_receptor   = False

                    for input_id in self.predecessors(entry_id):
                        input_genes = set(self._extract_genes(input_id))
                        if len(input_genes & entry_ligands) > 0:
                            input_ligands[input_id] = input_genes
                            is_receptor = True

                    if is_receptor:
                        mapper[entry_id] = [ f'{gene}__{self.kgml.KEGGID}__{entry_id}'
                                             for gene in entry_genes ]

                        for ligands in input_ligands.values():
                            self.LRdata += [ {'ligand':   ligand,
                                              'receptor': receptor.split('__')[0],
                                              'receptor_as_node': receptor}
                                             for ligand, receptor in product(ligands, mapper[entry_id]) ]
                    else:
                        mapper[entry_id] = list(entry_genes)

        if len(self.LRdata) > 0:
            self.LRdata = pd.DataFrame(self.LRdata).drop_duplicates().reset_index(drop=True)
        else:
            self.LRdata = pd.DataFrame(columns=['ligand','receptor','receptor_as_node'])

        # Re-construct a directed graph using genes.
        for entry_id, entry_genes in mapper.items():
            inputs   = list({ entry_in  for entry_in  in self.predecessors(entry_id) if entry_in  != entry_id })
            outputs  = list({ entry_out for entry_out in self.successors(entry_id)   if entry_out != entry_id })
            has_loop = self.has_edge(entry_id, entry_id)

            self.remove_node(entry_id)

            if len(entry_genes) == 0:
                self.add_edges_from( edge for edge in product(inputs, outputs) )
            else:
                self.add_edges_from( edge for edge in product(inputs, entry_genes)  )
                self.add_edges_from( edge for edge in product(entry_genes, outputs) )
                if has_loop:
                    self.add_edges_from( edge for edge in product(entry_genes, entry_genes) )

                self.add_nodes_from(entry_genes)

        # The indirect effects.
        indirect_effects = [ (mapper[entry_in], mapper[entry_out]) for entry_in, entry_out in self.indirect_effects ]
        indirect_effects = [ {'From': edge[0], 'To': edge[1]} for groups_in, groups_out in indirect_effects
                                                              for edge in product(groups_in, groups_out) ]
        if len(indirect_effects) == 0:
            self.indirect_effects = pd.DataFrame(columns=['From','To'])
        else:
            self.indirect_effects = pd.DataFrame(indirect_effects).drop_duplicates().reset_index(drop=True)


    def _extract_genes(self, entry_id):
        '''
        This method returns genes contained in <entry> of a given entry ID.
        '''
        entry_name = self.kgml.entry_data.loc[entry_id,'name']
        return re.findall(self.kgml.organism + r':[^ +]+', entry_name)


    def __repr__(self):
        return f'<KGMLGeneGraph id="{self.kgml.KEGGID}">'