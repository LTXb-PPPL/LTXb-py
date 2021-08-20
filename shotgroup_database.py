import numpy as np


def get_shots(key='011020'):
	if key == '010820':  # 08Jan20 first beam into plasma!
		return np.array([101585, 101588, 101594])
	if key == '011020':  # 10Jan20 perveance scan
		shots = 101600 + np.array([60, 61, 62, 63, 64, 65, 66, 67, 68, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79])
		return shots
	if key == '011320':  # 13Jan20 thompson data collecting
		return 101000 + np.array(
			[693, 694, 696, 699, 701, 704, 707, 708, 710, 711, 712, 713, 714, 716, 717, 718, 719, 720, 721])
	if key == 'cold_notable':
		shots = 101000 + np.array(
			[611, 621, 624, 625, 694, 707, 708, 710, 711, 718, 719, 699, 701, 704, 711, 712, 713, 720, 721, 731, 732,
			 742, 749, 750, 754, 733, 736, 746, 747, 748, 757, 776, 778, 785, 787, 793, 763, 781, 782])
		return shots
	if key == 'hot_notable':
		shots = 101000 + np.array(
			[818, 819, 826, 828, 845, 846, 820, 821, 829, 830, 844, 869, 870, 871, 872, 873, 874, 896, 897, 902, 903,
			 898, 899, 904, 911, 913])
		return shots
	if key == 'pauls_paper':
		shots = 100900 + np.array([81, 85, 88, 89])
		return shots

# 14Jan20 more thompson data collecting
nbi_shots = 101700 + np.array([30, 46, 47, 48, 49, 50, 51, 54, 57])

# good 40A shots to average from 13&14Jan
nbi_shots = 101700 + np.array([7, 8, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 30, 31, 32, 33, 36, 37, 42])
