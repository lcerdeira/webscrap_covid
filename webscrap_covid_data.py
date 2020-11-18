#!/usr/bin/env python3
import requests
import json
import os
import networkx
import pandas as pd
from ete3 import Tree
from collections import defaultdict
from Bio import Phylo
from networkx.readwrite.graphml import write_graphml

"""
    Script to extract the tree nodes from the NextStrain JSON output
    ("http://data.nextstrain.org/ncov.json") and Phylogenetic tree in Newick format with branch lengths
    in units of divergence (nextstrain_ncov_global_tree.nwk).

"""

def get_node_values(nodes, parent=None, parent_country=None):
    """
    Extract the tree nodes from the NextStrain JSON output
    ("http://data.nextstrain.org/ncov.json") into a list of
    dictionaries. Keeps country and parent node data.

    Notes: Now this raw json file is GZIP, treat before!
    """
    rows = []
    row = {}
    for node in nodes:
        row['name'] = node['name']
        row['parent'] = parent
        
        n = node['node_attrs']
        
        row['div'] = node['node_attrs']['div']
        row['date'] = node['node_attrs']['num_date']['value']
        row['date_lower'], row['date_upper'] = node['node_attrs']['num_date']['confidence']
        # row['country'] = n['country']['value']
        row['country'] = node['node_attrs']['country']['value']
        if 'confidence' in node['node_attrs']['country']:
            row['country_confidence'] = node['node_attrs']['country']['confidence'][row['country']]
        else:
            row['country_confidence'] = None
            
        row['country_entropy'] = node['node_attrs']['country'].get('entropy', None)
        
        if parent_country is not None:
            row['parent_country'] = parent_country['value']
            row['parent_country_confidence'] = parent_country['confidence'][row['parent_country']]
            row['parent_country_entropy'] = parent_country['entropy']
            
        rows.append(row)
        
        if 'children' in node:
            rows.extend(get_node_values(node['children'], 
                                        parent=row['name'], 
                                        parent_country=node['node_attrs']['country']))
    
    seen_node_names = []
    deduped_rows = []
    for r in rows:
        name = r['name']
        if name not in seen_node_names:
            deduped_rows.append(r)
            seen_node_names.append(name)

    return deduped_rows

# metadata_limited = pd.read_table(os.path.expanduser('nextstrain_ncov_global_metadata.tsv'))

# for strain in metadata_limited['Strain']:
#      if 'NODE_' in strain:
#          print(strain)

if not os.path.exists('ncov.json'):
    # url = "http://data.nextstrain.org/ncov.json"
    url = "https://nextstrain.org/charon/getDataset?prefix=/ncov/global" 

    response = requests.get(url, stream=True)
    data = response.json()

    with open('ncov.json', 'w') as fh:
        fh.write(json.dumps(data))
else:
    with open('ncov.json', 'r') as fh:
         data = json.loads(fh.read())

metadata = pd.DataFrame(get_node_values(data['tree']['children']))
# metadata.head(5)

tree = Tree(os.path.expanduser('nextstrain_ncov_global_tree.nwk'), format=1)
       
"""
Series has more than zero or one values.
"""
def first_or_none(series):
    if len(series) == 0:
        return None
    elif len(series) == 1:
        return series.values[0]
    else:
        raise ValueError("Series has more than zero or one values.")

international_events = []
local_pairs = []

for node in tree.traverse('preorder'):
    if node.up is None:
        continue

    parent_name = node.up.name
    country = metadata[metadata['name'] == node.name]['country']
    parent_country = metadata[metadata['name'] == node.name]['parent_country']
    date = metadata[metadata['name'] == node.name]['date']
    date_lower = metadata[metadata['name'] == node.name]['date_lower']
    date_upper = metadata[metadata['name'] == node.name]['date_upper']
    div = metadata[metadata['name'] == node.name]['div']
    country_entropy = metadata[metadata['name'] == node.name]['country_entropy']
    
    parent_country = first_or_none(parent_country)
    country = first_or_none(country)
    date = first_or_none(date)
    date_lower = first_or_none(date_lower)
    date_upper = first_or_none(date_upper)
    div = first_or_none(div)
    country_entropy = first_or_none(country_entropy)
    
    # we are at the root, no parent, skip to next node

    row = (parent_name,
           node.name,
           parent_country,
           country,
           str(date),
           date_lower,
           date_upper,
           div,
           country_entropy,
           node.dist,
          )
    
    if (country != parent_country
        #and country is not None
        ):
        international_events.append(row)
        
    if (country == parent_country
        #and country is not None
        ):
        local_pairs.append(row)

internat_df = pd.DataFrame(international_events, 
                           columns=['parent_strain', 'strain', 'parent_country', 'country', 'date', 'date_lower','date_upper','div','country_entropy', 'dist'])
internat_df.to_csv('international_events.tsv', sep='\t', index=False)

local_df = pd.DataFrame(local_pairs, 
                        columns=['parent_strain', 'strain', 'parent_country', 'country', 'date', 'date_lower','date_upper','div','country_entropy', 'dist'])

local_df.to_csv('local_pairs.tsv', sep='\t', index=False)

filter_internal_nodes = False

n_pop = len(list(tree.traverse()))

for e in international_events:
    src_name = e[0]
    dst_name = e[1]
    src_country = e[2]
    dst_country = e[3]
    
    if filter_internal_nodes and src_name.startswith('NODE_'):
        continue
    
    nodes = tree.search_nodes(name=src_name)
    
    if len(nodes) > 1:
        raise ValueError("Whaaa ?!?")
    
    if filter_internal_nodes:
        desc_count = 0
        for desc in nodes[0].iter_descendants():
            if not desc.name.startswith('NODE_'):
                desc_count += 1
    else:
        desc_count = len(list(nodes[0].iter_descendants()))
    
    proportion = desc_count / n_pop
    print(f'{src_name} ({src_country}) -> \t{dst_name} ({dst_country}): \t{desc_count} ({proportion * 100:.3} %)')

country_counts = defaultdict(int)
for node in tree.traverse():
    country =  first_or_none(metadata[metadata['name'] == node.name]['country'])
    country_counts[country] += 1
    
country_counts

filter_internal_nodes = False

n_pop = len(list(tree.traverse()))

for i, row in internat_df.iterrows():
    src_name = row['parent_strain']
    dst_name = row['strain']
    src_country = row['parent_country']
    dst_country = row['country']
    
    nodes = tree.search_nodes(name=src_name)
    
    if len(nodes) > 1:
        raise ValueError("Whaaa ?!?")
    
    desc_count = len(list(nodes[0].iter_descendants()))
    
    total_proportion = desc_count / n_pop
    internat_df.at[i, 'desc_count'] = desc_count
    internat_df.at[i, 'total_proportion'] = total_proportion
    
    n_same_country_desc = 0
    for n in nodes[0].iter_descendants():
        country =  first_or_none(metadata[metadata['name'] == n.name]['country'])
        if country == dst_country:
            n_same_country_desc += 1
        
    country_prop = n_same_country_desc / country_counts[dst_country]
    internat_df.at[i, 'country_proportion'] = country_prop

print(f'{src_name} ({src_country}) -> \t{dst_name} ({dst_country}): \t{desc_count} ({proportion * 100:.3} %)')

internat_df.to_csv('international_events.tsv', sep='\t', index=False)
# internat_df

#print international events (filtered)
internat_df[internat_df['country'] == 'Australia']

#print local events (filtered)
local_df[local_df['country'] == 'Australia']

#Alternative with Biopython: https://biopython.org/wiki/Phylo - read tree, #convert to networkx graph.

tree = next(Phylo.parse('nextstrain_ncov_global_tree.nwk', 'newick'))
G = Phylo.to_networkx(tree)
write_graphml(G, 'nextree_global.graphml')