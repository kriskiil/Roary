#!/usr/bin/env python3
# Copyright (C) <2015> EMBL-European Bioinformatics Institute

# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# Neither the institution name nor the name roary_plots
# can be used to endorse or promote products derived from
# this software without prior written permission.
# For written permission, please contact <marco@ebi.ac.uk>.

# Products derived from this software may not be called roary_plots
# nor may roary_plots appear in their names without prior written
# permission of the developers. You should have received a copy
# of the GNU General Public License along with this program.
# If not, see <http://www.gnu.org/licenses/>.

__author__ = "Marco Galardini"
__version__ = '0.1.0'

def get_options():
    import argparse

    # create the top-level parser
    description = "Create plots from roary outputs"
    parser = argparse.ArgumentParser(description = description,
                                     prog = 'roary_plots.py')

    #parser.add_argument('tree', action='store',
    #                    help='Newick Tree file', default='accessory_binary_genes.fa.newick')
    parser.add_argument('spreadsheet', action='store',
                        help='Roary gene presence/absence spreadsheet', default='gene_presence_absence.csv')

    parser.add_argument('--labels', action='store_true',
                        default=False,
                        help='Add node labels to the tree (up to 10 chars)')
    parser.add_argument('--format',
                        choices=('png',
                                 'tiff',
                                 'pdf',
                                 'svg'),
                        default='png',
                        help='Output format [Default: png]')
    parser.add_argument('-N', '--skipped-columns', action='store',
                        type=int,
                        default=14,
                        help='First N columns of Roary\'s output to exclude [Default: 14]')
    parser.add_argument('--tree-only', action='store_true', help="Only save a tree in phylip format")
    
    parser.add_argument('--version', action='version',
                         version='%(prog)s '+__version__)

    return parser.parse_args()


def ctree_to_phylotree(t,nodenames=None):
    """Creates a Bio.Phylo.Tree from a Bio.Cluster.Tree using single linkage clustering"""
    nodes = list()
    for n in t:
        node = Phylo.BaseTree.Clade()
        for child in (n.left, n.right):
            if child >= 0: # Leaf node
                sub_node = Phylo.BaseTree.Clade(branch_length=n.distance, name=nodenames[child])
                sub_node.cumlength = n.distance
            else: # Internal node
                sub_node = nodes[-child-1]
                sub_node.branch_length = n.distance-max(
                    [c.cumlength for c in sub_node.clades])
                sub_node.cumlength = n.distance
            node.clades.append(sub_node)
        nodes.append(node)
    return Phylo.BaseTree.Tree(root=nodes[-1])
    
def binary_distance_tree(m, nodenames=None):
    """Calculates a tree from single linkage clustering using Jaccard distances."""
    n = m.shape[1]
    if nodenames is None:
        nodenames = list(range(n))
    distmat = squareform(pdist(m.transpose(), 'jaccard'))
    
    print(len(distmat))
    ctree = Cluster.treecluster(data=None,distancematrix=distmat,method='s')
    return ctree_to_phylotree(ctree, nodenames)
    

if __name__ == "__main__":
    options = get_options()

    import matplotlib
    matplotlib.use('Agg')

    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_style('white')

    import os
    import sys
    import pandas as pd
    import numpy as np
    from Bio import Phylo
    from Bio import Cluster
    from scipy.spatial.distance import pdist, jaccard
    from scipy.spatial.distance import squareform


    # Load roary
    roary = pd.read_table(options.spreadsheet,
                         sep='\t',
                         #sep=',',
                         low_memory=False)
    # Set index (group name)
    roary.set_index('Gene', inplace=True)
    # Drop the other info columns
    roary.drop(list(roary.columns[:options.skipped_columns-1]), axis=1, inplace=True)

    # Transform it in a presence/absence matrix (1/0)
    roary.replace('.{2,100}', 1, regex=True, inplace=True)
    roary.replace(np.nan, 0, regex=True, inplace=True)
    roary = roary.astype(dtype="bool")

    # Sort the matrix by the sum of strains presence
    # idx = roary.sum(axis=1).sort_values(ascending=False).index
    # roary_sorted = roary.loc[idx]

    # # Pangenome frequency plot
    # plt.figure(figsize=(7, 5))

    # plt.hist(roary.sum(axis=1), roary.shape[1],
    #          histtype="stepfilled", alpha=.7)

    # plt.xlabel('No. of genomes')
    # plt.ylabel('No. of genes')

    # sns.despine(left=True,
    #             bottom=True)
    # plt.savefig('pangenome_frequency.%s'%options.format, dpi=300)
    # plt.clf()


    #t = Phylo.read(options.tree, 'newick')
    t = binary_distance_tree(roary, list(roary.keys()))
    with open("pangenome.distance.nwk","w") as fh:
        Phylo.NewickIO.write([t],fh)
    if options.tree_only == True:
        sys.exit()
    # Sort the matrix according to tip labels in the tree
    roary_sorted = roary[[x.name for x in t.get_terminals()]]
    
    troary = roary_sorted.transpose()
    t_transposed = binary_distance_tree(troary, list(troary.keys()))

    # Sort the matrix according to tip labels in the tree
    troary_sorted = troary[[x.name for x in t_transposed.get_terminals()]]

    roary_sorted = troary_sorted.transpose()

    # Max distance to create better plots
    mdist = max([t.distance(t.root, x) for x in t.get_terminals()])




    # Plot presence/absence matrix against the tree
    with sns.axes_style('whitegrid'):
        fig = plt.figure(figsize=(17, 10))

        ax1=plt.subplot2grid((1,40), (0, 10), colspan=30)
        a=ax1.matshow(roary_sorted.T, cmap=plt.cm.Blues,
                   vmin=0, vmax=1,
                   aspect='auto',
                   interpolation='none',
                    )
        ax1.set_yticks([])
        ax1.set_xticks([])
        ax1.axis('off')

        ax = fig.add_subplot(1,2,1)
        # matplotlib v1/2 workaround
        try:
            ax=plt.subplot2grid((1,40), (0, 0), colspan=10, facecolor='white')
        except AttributeError:
            ax=plt.subplot2grid((1,40), (0, 0), colspan=10, axisbg='white')

        fig.subplots_adjust(wspace=0, hspace=0)

        ax1.set_title('Roary matrix\n(%d gene clusters)'%roary.shape[0])

        if options.labels:
            fsize = 12 - 0.1*roary.shape[1]
            if fsize < 7:
                fsize = 7
            with plt.rc_context({'font.size': fsize}):
                Phylo.draw(t, axes=ax, 
                           show_confidence=False,
                           label_func=lambda x: str(x)[:10],
                           xticks=([],), yticks=([],),
                           ylabel=('',), xlabel=('',),
                           xlim=(-mdist*0.1,mdist+mdist*0.45-mdist*roary.shape[1]*0.001),
                           axis=('off',),
                           title=('Tree\n(%d strains)'%roary.shape[1],), 
                           do_show=False,
                          )
        else:
            Phylo.draw(t, axes=ax, 
                       show_confidence=False,
                       label_func=lambda x: None,
                       xticks=([],), yticks=([],),
                       ylabel=('',), xlabel=('',),
                       xlim=(-mdist*0.1,mdist+mdist*0.1),
                       axis=('off',),
                       title=('Tree\n(%d strains)'%roary.shape[1],),
                       do_show=False,
                      )
        plt.savefig('pangenome_matrix.%s'%options.format, dpi=300)
        plt.clf()

    # Plot the pangenome pie chart
    plt.figure(figsize=(10, 10))

    core     = roary[(roary.sum(axis=1) >= roary.shape[1]*0.99) & (roary.sum(axis=1) <= roary.shape[1]     )].shape[0]
    softcore = roary[(roary.sum(axis=1) >= roary.shape[1]*0.95) & (roary.sum(axis=1) <  roary.shape[1]*0.99)].shape[0]
    shell    = roary[(roary.sum(axis=1) >= roary.shape[1]*0.15) & (roary.sum(axis=1) <  roary.shape[1]*0.95)].shape[0]
    cloud    = roary[roary.sum(axis=1)  < roary.shape[1]*0.15].shape[0]

    total = roary.shape[0]
    
    def my_autopct(pct):
        val=int(round(pct*total/100.0))
        return '{v:d}'.format(v=val)

    a=plt.pie([core, softcore, shell, cloud],
          labels=['core\n(%d <= strains <= %d)'%(roary.shape[1]*.99,roary.shape[1]),
                  'soft-core\n(%d <= strains < %d)'%(roary.shape[1]*.95,roary.shape[1]*.99),
                  'shell\n(%d <= strains < %d)'%(roary.shape[1]*.15,roary.shape[1]*.95),
                  'cloud\n(strains < %d)'%(roary.shape[1]*.15)],
          explode=[0.1, 0.05, 0.02, 0], radius=0.9,
          colors=[(0, 0, 1, float(x)/total) for x in (core, softcore, shell, cloud)],
          autopct=my_autopct)
    plt.savefig('pangenome_pie.%s'%options.format, dpi=300)
    plt.clf()
