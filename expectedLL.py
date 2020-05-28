from dendropy import Tree, DnaCharacterMatrix
import numpy
import math
import dendropy
import numpy.linalg as la




def give_index(c):
    if c == "A":
        return 0
    elif c == "C":
        return 1
    elif c == "G":
        return 2
    elif c == "T":
        return 3
# =======================================================================================================
def p_matrix(br_length,model):
    p = numpy.zeros((4, 4))
    if model == 'JC69':
        for i in range(4):
            for j in range(4):
                p[i][j] = 0.25 - 0.25 * math.exp(-4 * br_length / 3)
            p[i][i] = 0.25 + 0.75 * math.exp(-4 * br_length / 3)
    elif model == 'GTR':
        mu = 0
        freq = numpy.zeros((4, 4))
        q = numpy.zeros((4, 4))
        sqrtPi = numpy.zeros((4, 4))
        sqrtPiInv = numpy.zeros((4, 4))
        exchang = numpy.zeros((4, 4))
        s = numpy.zeros((4, 4))
        fun = numpy.zeros(1)
        a, b, c, d, e = rates
        f = 1

        freq = numpy.diag(pi)
        sqrtPi = numpy.diag(numpy.sqrt(pi))
        sqrtPiInv = numpy.diag(1.0 / numpy.sqrt(pi))
        mu = 1 / (2 * ((a * pi[0] * pi[1]) + (b * pi[0] * pi[2]) + (c * pi[0] * pi[3]) + (d * pi[1] * pi[2]) + (
                e * pi[1] * pi[3]) + (pi[2] * pi[3])))
        exchang[0][1] = exchang[1][0] = a
        exchang[0][2] = exchang[2][0] = b
        exchang[0][3] = exchang[3][0] = c
        exchang[1][2] = exchang[2][1] = d
        exchang[1][3] = exchang[3][1] = e
        exchang[2][3] = exchang[3][2] = f

        q = numpy.multiply(numpy.dot(exchang, freq), mu)

        for i in range(4):
            q[i][i] = -sum(q[i][0:4])

        s = numpy.dot(sqrtPi, numpy.dot(q, sqrtPiInv))

        eigval, eigvec = la.eig(s)
        eigvec_inv = la.inv(eigvec)

        left = numpy.dot(sqrtPi, eigvec)
        right = numpy.dot(eigvec_inv, sqrtPiInv)

        p = numpy.dot(left, numpy.dot(numpy.diag(numpy.exp(eigval * br_length)), right))

    return p
# =======================================================================================================
def expectedLL(tree,dna,model):

    partial = numpy.zeros((2 * tips,4))
    expll = numpy.zeros((2 * tips, 4))

    for node in tree.postorder_node_iter():
        node.index = 0
        node.annotations.add_bound_attribute("index")

    s = tips + 1
    for id,node in enumerate(tree.postorder_node_iter()):
        if not node.is_leaf():
            node.index = s
            s += 1
        else:
            for idx, name in enumerate(dna):
                if idx + 1 == int(node.taxon.label):
                    node.index = idx+1
                    break
    pos = 0
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            i = give_index(str(dna[pos]))
            pos += 1
            partial[node.index][i] = 1
            expll[node.index][i] = 1
        else:
            # print(node.index)
            children = node.child_nodes()
            partial[node.index] = numpy.dot(p_matrix(children[0].edge_length,model), partial[children[0].index])
            expll[node.index] = partial[node.index] * children[0].edge_length
            for i in range(1, len(children)):
                partial[node.index] *= numpy.dot(p_matrix(children[i].edge_length,model), partial[children[i].index])
                expll[node.index] *= (numpy.dot(p_matrix(children[i].edge_length, model), partial[children[i].index]) * children[i].edge_length)

    expected_ll = numpy.zeros(tips - 1)
    p_index = 0
    for par  in range(tips+1,tree.seed_node.index+1,1):
        temp = 0
        temp = numpy.sum(expll[par]) * 0.25
        expected_ll[p_index] = round(numpy.log(temp), 7)
        p_index += 1


    return  expected_ll
#=======================================================================================================================
def fillvector(tree,alignment,model):
    expected_vector_a = numpy.zeros((tips - 1 ,alignment_len ))
    expected_vector_b = numpy.zeros((tips - 1, alignment_len))
    for l in range(alignment_len):
        col = ""
        for t in range(tips):
            col += str(alignment[t][l])
        expected_vector_b[:, l] = expectedLL(tree, col, model)
    return expected_vector_b
#=======================================================================================================================
pi = [0.2184,0.2606,0.3265,0.1946]
rates = [2.0431,0.0821,0,0.067,0]
f = 1
tree = Tree.get_from_path('/home/nehleh/0_Research/PhD/Data/simulationdata/recombination/faketree/RAxML_bestTree.RaxMLtree', 'newick')
alignment = dendropy.DnaCharacterMatrix.get(file=open("/home/nehleh/0_Research/PhD/Data/simulationdata/recombination/faketree/recombined.fasta"), schema="fasta")
tips = len(alignment)
alignment_len = alignment.sequence_size


result = fillvector(tree,alignment,'GTR')
for i in range(result.shape[0]):
    avg = numpy.sum(result[i]) /alignment_len
    print(" expected likelihood node {}  is: {}".format(i+1+tips, avg))


