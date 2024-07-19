#!/usr/bin/env python
import os
import _pickle as pickle
from collections import namedtuple

Node = namedtuple('Node', ['taxid', 'name', 'rank', 'parent', 'children', 'aliases'])
parse_line = lambda l: list(map(lambda s: s.strip(), l.split('|')))

class TaxonomyTree():

  def __init__(self, taxdump_dir, force_reload=False):
    """build taxonomy tree from nodes.dmp and names.dmp
    
    >>> tree = taxtree.TaxonomyTree('taxdmp')  
    """
    pickle_fn, names_fn, merged_fn, nodes_fn = map(lambda fn: os.path.join(taxdump_dir, fn),
        ('tree.pkl', 'names.dmp', 'merged.dmp', 'nodes.dmp'))

    if os.path.exists(pickle_fn) and not force_reload:
      with open(pickle_fn, 'rb') as f:
        self.nodes, self.aliases = pickle.load(f)
    else:
      names = {}
      with open(names_fn) as f:
        for l in f:
          taxid, name, _, name_class = parse_line(l)[:4]
          if name_class == 'scientific name':
            names[int(taxid)] = name
      
      # track all nodes whose taxids have been changed
      with open(merged_fn) as f:
        self.aliases = {old: new for old, new in 
                  map(lambda l: map(int, parse_line(l)[:2]), f)}
        merged = {}
        for old, new in self.aliases.items():
          if new not in merged: merged[new] = {old}
          else: merged[new].add(old)

      with open(nodes_fn) as f:
        self.nodes = {}
        unseen = {}
        for l in f:
          child_, parent_, rank = parse_line(l)[:3]
          child, parent = int(child_), int(parent_) if parent_.isdigit() else None

          self.nodes[child] = Node(
              taxid=child,
              name=names.get(child, None),
              rank=rank,
              parent=parent,
              children=set(), 
              aliases=merged.get(child, None))
          # node 1 has says its a child of itself
          if parent != child:
            if parent not in self.nodes:
              # some parent nodes are referenced before they come up as children
              if parent not in unseen:
                unseen[parent] = set((child,))
              else:
                unseen[parent].add(child)
            else:
              self.nodes[parent].children.add(child)
        # second pass creating links between unseen parents
        for parent in unseen:
          self.nodes[parent].children.update(unseen[parent])

      with open(pickle_fn, 'wb') as f:
        pickle.dump((self.nodes, self.aliases), f)

  def __getitem__(self, taxid):
    """implements `tree[taxid]` aka dict access, considering aliases
    
    >>> tree[11207]
    Node(taxid=1979162, name=None, rank='species', parent=39744, children={33730, 11208, 31608, 31609, 31610}, aliases={11200, 12837, 11207})
    """
    if taxid in self.nodes: return self.nodes[taxid]
    elif taxid in self.aliases: return self.__getitem__(self.aliases[taxid])
    else: raise KeyError('Could not find taxid: {}'.format(taxid))


  def __contains__(self, taxid):
    """implements `taxid in tree`, considering aliases"""
    return taxid in self.nodes or taxid in self.aliases
  
  def ascend(self, taxid):
    """walk up the tree and return list of nodes encountered"""
    return [self[taxid]] + (
        self.ascend(self[taxid].parent) if taxid in self and taxid is not self[taxid].parent 
        else [])
  def descend(self, taxid):
    """walk down the tree and return set of all child nodes"""
    all_children = {taxid}
    node = self[taxid]
    for child in node.children:
      all_children.update(self.descend(child))
    return all_children


if __name__ == '__main__':
  import time
  print('starting cold load...')
  t = time.time()
  tree = TaxonomyTree('taxdmp', force_reload=True)
  print('cold load in {} secs...'.format(time.time() - t))
  print('starting warm load...')
  t = time.time()
  tree = TaxonomyTree('taxdmp')
  print('warm load in {} secs...'.format(time.time() - t))
