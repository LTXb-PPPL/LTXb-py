import matplotlib.pyplot as plt

fig = plt.figure()
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222, projection='polar')
ax3 = fig.add_subplot(224, projection='polar')
ax4 = fig.add_subplot(223)
plt.tight_layout()
plt.show()
# eq_fn = 'Z:/transp/t106536/106536R02_05.eqdsk'
# read_eqdsk3(eq_fn, plot=True)

if __name__ == '__main__':
	pass