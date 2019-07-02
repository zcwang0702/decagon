from collections import defaultdict
import networkx as nx


# Returns dictionary from combination ID to pair of stitch IDs,
# dictionary from combination ID to list of polypharmacy side effects, 
# and dictionary from side effects to their names.
# bio-decagon-combo.csv: STITCH 1; STITCH 2; Polypharmacy Side Effect ID; Side Effect Name
def load_combo_se(fname):
    """
    :param fname: bio-decagon-combo.csv, Polypharmacy side effects triples
    :return: combo2stitch, drug pair
             combo2se, side effect list for certain drug
             e2name, side effect list
    """
    combo2stitch = {}
    combo2se = defaultdict(set)
    se2name = {}
    fin = open(fname)
    print('***'*10, 'Reading: %s' % fname, '***'*10)
    fin.readline()
    for line in fin:
        stitch_id1, stitch_id2, se, se_name = line.strip().split(',')
        combo = stitch_id1 + '_' + stitch_id2
        combo2stitch[combo] = [stitch_id1, stitch_id2]
        combo2se[combo].add(se)
        se2name[se] = se_name
    fin.close()
    n_interactions = sum([len(v) for v in combo2se.values()])

    all_drug = []
    for drug_pair in combo2stitch.values():
        all_drug += drug_pair

    print('Drug combinations number: %d Unique side effects number: %d' % (len(combo2stitch), len(se2name)))
    print('Drug-drug interactions number: %d' % n_interactions)
    print('Unique drug number: %d' % len(set(all_drug)))
    return combo2stitch, combo2se, se2name


# Returns networkx graph of the PPI network
# and a dictionary that maps each gene ID to a number
# bio-decagon-ppi.csv: Protein-protein interaction network
def load_ppi(fname):
    fin = open(fname)
    print('***'*10, 'Reading: %s' % fname, '***'*10)
    fin.readline()
    edges = []
    for line in fin:
        gene_id1, gene_id2 = line.strip().split(',')  # two-protein pair
        edges += [[gene_id1, gene_id2]]
    nodes = set([u for e in edges for u in e])
    print('PPI Edges: %d' % len(edges))
    print('PPI Nodes: %d' % len(nodes))
    net = nx.Graph()
    net.add_edges_from(edges)
    net.remove_nodes_from(nx.isolates(net))
    net.remove_edges_from(net.selfloop_edges())
    node2idx = {node: i for i, node in enumerate(net.nodes())}  # protein id : index
    return net, node2idx


# Returns dictionary from Stitch ID to list of individual side effects,
# and dictionary from side effects to their names
# bio-decagon-mono.csv: Side effects of individual drugs
def load_mono_se(fname):
    stitch2se = defaultdict(set)
    se2name = {}
    fin = open(fname)
    print('***'*10, 'Reading: %s' % fname, '***'*10)
    fin.readline()
    for line in fin:
        contents = line.strip().split(',')
        stitch_id, se, = contents[:2]
        se_name = ','.join(contents[2:])
        stitch2se[stitch_id].add(se)
        se2name[se] = se_name
    n_interactions = sum([len(v) for v in stitch2se.values()])
    print('Individual drug side effect total number: %d' % n_interactions)
    print('Unique drug number: %d' % len(stitch2se))

    all_se = []
    for se_dict in stitch2se.values():
        all_se += list(se_dict)

    print('Unique side effect number : %d' % len(set(all_se)))
    return stitch2se, se2name


# Returns dictionary from Stitch ID to list of drug targets
# bio-decagon-targets.csv	Drug-target protein associations, used for model training
# bio-decagon-targets-all.csv	Drug-target protein associations culled from several curated databases, full dataset
# STITCH: drug ID; Gene: protein ID
def load_targets(fname):
    stitch2proteins = defaultdict(set)
    fin = open(fname)
    print('***'*10, 'Reading: %s' % fname, '***'*10)
    fin.readline()
    for line in fin:
        stitch_id, gene = line.strip().split(',')
        stitch2proteins[stitch_id].add(gene)
    n_interactions = sum([len(v) for v in stitch2proteins.values()])
    print('Drug-protein interaction: %d' % n_interactions)
    return stitch2proteins


# Returns dictionary from side effect to disease class of that side effect,
# and dictionary from side effects to their names.
# bio-decagon-effectcategories.csv: Side effect categories
# Side Effect ID	Side Effect Name	Disease Class

def load_categories(fname):
    se2name = {}
    se2class = {}
    fin = open(fname)
    print('***'*10, 'Reading: %s' % fname, '***'*10)
    fin.readline()
    for line in fin:
        se, se_name, se_class = line.strip().split(',')
        se2name[se] = se_name
        se2class[se] = se_class

    print('unique side effect class number: %d' % len(set(se2class.values())))
    print('unique side effect name: %d' % len(set(se2name.values())))
    return se2class, se2name
