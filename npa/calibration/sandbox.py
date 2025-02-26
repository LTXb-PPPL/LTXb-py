import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.use('TkAgg')  # allows plotting in debug mode
clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']


student_tuples = [
    ('john', 'A', 15),
    ('jane', 'B', 12),
    ('dave', 'B', 10),
]
print(sorted(student_tuples, key=lambda butt: butt[2]))

if __name__ == '__main__':
	pass
