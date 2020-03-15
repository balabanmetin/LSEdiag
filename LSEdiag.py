import treeswift as tw
import subprocess
import tempfile
import re
from io import StringIO
from optparse import OptionParser
from shutil import copyfile

def write_phylip_dist(obs_dist):
    output = StringIO()
    output.write(str(len(obs_dist))+"\n")
    keys = obs_dist.keys()
    for i in keys:
        output.write(str(i))
        for j in keys:
            output.write("\t" + str(obs_dist[i][j]))
        output.write("\n")
    return output.getvalue()
            


if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-d", "--distances", dest="dist_fp",
                      help="path to the table of observed distances", metavar="FILE")
    parser.add_option("-t", "--tree", dest="tree_fp",
                      help="path to the newick tree (optional). If not provided, a tree is inferred from given distance matrix.", metavar="FILE")
    parser.add_option("-r", "--treeout", dest="treeout_fp",
                      help="path to the directory where the phylogeny should be output. If an input tree provided, it will be the same tree. If not, it will be the tree inferred by FastME", metavar="FILE")


    (options, args) = parser.parse_args()
    dist_fp = options.dist_fp
    tree_fp = options.tree_fp
    treeout_fp = options.treeout_fp


    tbl = open(dist_fp)
    tags = list(re.split("\s+", tbl.readline().rstrip()))[1:]
    obs_dist=dict()

    for line in tbl.readlines():
            dists = list(re.split("\s+", line.strip()))
            query_name = dists[0]
            obs_dist[query_name] = dict(zip(tags, map(float, dists[1:])))
    
    tbl.close()

    
    if not tree_fp:
        tree_fp = tempfile.NamedTemporaryFile(delete=True, mode='w+t').name
        dist_phy = tempfile.NamedTemporaryFile(delete=True, mode='w+t')
        nldef = tempfile.NamedTemporaryFile(delete=True, mode='w+t')

        dist_phy.write(write_phylip_dist(obs_dist))
        dist_phy.flush()

        s = ["fastme", "-i", dist_phy.name, "-o", tree_fp]
        subprocess.call(s, stdout = nldef, stderr = nldef)
    
    if treeout_fp:
        copyfile(tree_fp, treeout_fp + "/tree.nwk")
    treestr = open(tree_fp).readline().strip()
    tree = tw.read_tree(treestr, "newick")
    pdc = tree.distance_matrix(leaf_labels=True)
     
    try:
        errs = dict()
        glob = 0
        glob_fm = 0
        errs_fm = dict()
        for l1 in tree.labels(leaves=True, internal=False):
            tot = 0
            tot_fm = 0
            for l2 in tree.labels(leaves=True, internal=False):
                if not l1 == l2:
                    cont =  (pdc[l1][l2] - obs_dist[l1][l2])**2
                    if cont > 0 and obs_dist[l1][l2] > 0:
                        tot += cont
                        glob += cont
                        tot_fm += cont/(pdc[l1][l2])**2
                        glob_fm += cont/(pdc[l1][l2])**2
            errs[l1] = tot
            errs_fm[l1] = tot_fm
        print("%s\t%s\t%s\t%s\t%s"% ("Organism", "FM", "FM-%", "LSE", "LSE-%"))
        for k,v in sorted(errs_fm.items(), key=lambda kv: -kv[1]/glob_fm):
            print("%s\t%.4f\t%.4f\t%.4f\t%.4f" % (k,v,v/glob_fm,errs[k],errs[k]/glob))
    except:
        import pdb,traceback; traceback.print_exc(); pdb.set_trace()
