import numpy as np
from kraken2ref.taxonomytree import TaxonomyTree
from kraken2ref.poll import Poll

def generate_random_data():
    bulk = np.random.randint(100, 1000, 90)
    lower_outliers = np.random.randint(50, 100, 7)
    upper_outliers = np.random.randint(1000, 10_000, 3)
    data = list(np.concatenate((bulk, lower_outliers, upper_outliers)))
    return data


random_freqs1 = [369, 891, 463, 281, 258, 534, 615, 855, 699, 936, 106, 498, 451, 744, 752, 341, 343, 613, 674, 523, 731, 464, 675, 215, 480, 798, 897, 516, 881, 610, 828, 793, 615, 234, 811, 334, 712, 505, 221, 149, 784, 466, 592, 343, 256, 628, 878, 716, 856, 138, 830, 137, 314, 150, 374, 662, 146, 561, 509, 124, 154, 894, 983, 886, 289, 238, 350, 696, 295, 969, 910, 752, 240, 128, 251, 781, 367, 445, 934, 848, 813, 880, 435, 929, 187, 293, 689, 338, 623, 197, 75, 89, 96, 93, 61, 74, 86, 8362, 7050, 5008]
random_freqs2 = [611, 208, 328, 159, 373, 820, 190, 597, 924, 147, 483, 197, 647, 923, 205, 244, 164, 465, 789, 552, 898, 501, 689, 990, 894, 746, 622, 294, 575, 400, 335, 709, 486, 859, 797, 254, 643, 653, 354, 711, 275, 279, 466, 400, 174, 663, 716, 376, 836, 993, 999, 912, 175, 860, 685, 901, 115, 759, 302, 746, 512, 483, 835, 373, 105, 898, 242, 996, 440, 282, 926, 489, 262, 292, 225, 938, 425, 840, 884, 104, 897, 495, 921, 197, 336, 769, 635, 513, 630, 713, 57, 86, 84, 69, 57, 85, 70, 1248, 8913, 4846]
random_freqs3 = [195, 566, 743, 599, 680, 348, 270, 386, 255, 728, 799, 998, 892, 107, 491, 131, 105, 982, 726, 536, 846, 765, 679, 233, 723, 150, 368, 691, 293, 810, 962, 833, 473, 481, 662, 825, 359, 644, 975, 435, 311, 504, 384, 195, 727, 544, 591, 273, 208, 268, 384, 223, 581, 370, 969, 637, 725, 585, 473, 829, 275, 135, 537, 992, 224, 441, 740, 933, 647, 493, 959, 693, 525, 448, 663, 280, 161, 987, 798, 658, 610, 899, 941, 577, 599, 782, 932, 680, 419, 349, 81, 53, 63, 95, 88, 93, 96, 8163, 4976, 9625]

def create_nodes_and_data_dict(freqs_list):
    random_data = [(0, 0, 8997), (0, 0, 8998), (0, 0, 8999)] + [(element, element, 9000+i) for i, element in enumerate(freqs_list)]
    stem = [(11, "S"), (12, "S1"), (13, "S2")]
    leaves = []
    i = 14
    while len(leaves) < 100:
        leaves.append((i, "S3"))
        i += 1
    nodes = stem + leaves
    data_dict = dict(zip(nodes, random_data))
    tree = TaxonomyTree(nodes=nodes)
    return tree, data_dict

def test_polling_with_skew():
    tree1, data_dict1 = create_nodes_and_data_dict(random_freqs1)
    poll_t100 = Poll(taxonomy_tree=tree1, data_dict=data_dict1, threshold=100)
    poll_t100.poll_leaves(method="skew")
    result_t100 = poll_t100.filt_leaves
    # print(poll_t100.class_method)
    assert result_t100 == [(111, 'S3')]

    poll_t900 = Poll(taxonomy_tree=tree1, data_dict=data_dict1, threshold=900)
    poll_t900.poll_leaves(method="skew")
    result_t900 = poll_t900.filt_leaves
    # print(poll_t900.class_method)
    assert sorted(result_t900) == sorted([(97, 'S3'), (92, 'S3'), (23, 'S3'), (83, 'S3'), (76, 'S3'), (113, 'S3'), (112, 'S3'), (111, 'S3')])

# test_polling_with_skew()

def test_polling_with_tiles():
    tree1, data_dict1 = create_nodes_and_data_dict(random_freqs1)
    poll_t100 = Poll(taxonomy_tree=tree1, data_dict=data_dict1, threshold=100)
    poll_t100.poll_leaves(method="tiles")
    result_t100 = poll_t100.filt_leaves
    assert result_t100 == [(111, 'S3'), (112, 'S3'), (113, 'S3')]

    poll_t900 = Poll(taxonomy_tree=tree1, data_dict=data_dict1, threshold=900)
    poll_t900.poll_leaves(method="tiles")
    result_t900 = poll_t900.filt_leaves
    assert result_t900 == []

# test_polling_with_tiles()

def test_polling_with_kmeans():
    tree1, data_dict1 = create_nodes_and_data_dict(random_freqs1)
    poll_t100 = Poll(taxonomy_tree=tree1, data_dict=data_dict1, threshold=100)
    poll_t100.poll_leaves(method="kmeans")
    result_t100 = poll_t100.filt_leaves
    assert result_t100 == [(111, 'S3'), (112, 'S3'), (113, 'S3')]

    poll_t900 = Poll(taxonomy_tree=tree1, data_dict=data_dict1, threshold=900)
    poll_t900.poll_leaves(method="kmeans")
    result_t900 = poll_t900.filt_leaves
    assert result_t900 == [(111, 'S3'), (112, 'S3'), (113, 'S3')]

# test_polling_with_kmeans()
