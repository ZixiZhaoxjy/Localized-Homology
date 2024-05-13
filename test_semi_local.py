# a demo for computing localized homology

from blowup_complex import BlowupComplex
from complex import DelaunayComplex
from simplicial_cover import HyperCube
import numpy as np
from matplotlib import pyplot as plt
from semi_blowup_complex import SemiBlowupComplex
import util
import persistence as persis
import cechmate as cm
import phat
import gudhi
from complex import HyperCube



if __name__ == '__main__':

    import numpy as np
    import gudhi as gd


    def compute_rips_filtration_components_with_lists(point_cloud, max_edge_length=1.0):


        rips_complex = gd.RipsComplex(points=point_cloud, max_edge_length=max_edge_length)

        simplex_tree = rips_complex.create_simplex_tree(max_dimension=3)

        filtration_dict = {}
        edges = []
        two_simplices = []

        for simplex, filtration_value in simplex_tree.get_filtration():
            simplex_list = list(simplex)
            filtration_dict[tuple(simplex)] = round(filtration_value, 2)
            if len(simplex) == 2:
                edges.append(simplex_list)
            elif len(simplex) == 3:  # two-simplex
                print(simplex_list)
                two_simplices.append(simplex_list)

        return edges, two_simplices, filtration_dict


    # 示例
    testlist = [[ 0,  2], [ 2,  0],[ 2,  2]]#[[0,0],[0,2],[2,0],[2,2]]
    testdata = np.array(testlist)
    point_cloud = [tuple(row) for row in testdata]

    # 计算并获取结果
    edges, two_simplices, dist = compute_rips_filtration_components_with_lists(point_cloud, max_edge_length=4)

    # print(dist,edges)
    a = HyperCube(testdata , edges,
                  two_simplices,  # edge_list
                  1,  # N
                  1,  # M, the numbers of covers is N*M = 2*1
                  0.5,  # overlapping degree(0-2)
                  True
                  )
    print(dist)
    print(a.cover())

    semi_blowup = SemiBlowupComplex(a.cover(), dist)
    semi_blowup.compute_persistence(verbose=True, show_diag = False)
    dgms = semi_blowup.dgms
    print(semi_blowup.dgms)

    fig, ax = plt.subplots(ncols=1, figsize=(10, 10))
    ax.scatter(testdata[:, 0], testdata[:, 1])
    plt.show()
    #plt.savefig('figures/test_semi/four_points.png')
