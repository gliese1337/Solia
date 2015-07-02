from itertools import combinations
import matplotlib.pyplot as plt
from mpmath import cbrt
import sys
import OrbitData

times = []
masses = None
positions = []
for time, bodies, _ in OrbitData.read(sys.stdin):
	times.append(time)
	positions.append(map(lambda b: b[1:3], bodies))
	masses = map(lambda b: b[0], bodies)

scale = min(masses)
masses = [12*cbrt(m/scale)**2 for m in masses]

print "Read all data"

def update_plot(i, scat, pos, mass):
	scat.set_offsets(pos[i])
	scat.set_sizes(mass)
	return scat

minx = min(map(lambda p: min(map(lambda c: c[0], p)), positions))
maxx = max(map(lambda p: max(map(lambda c: c[0], p)), positions))
miny = min(map(lambda p: min(map(lambda c: c[1], p)), positions))
maxy = max(map(lambda p: max(map(lambda c: c[1], p)), positions))

xbuf = (maxx - minx)/10
maxx += xbuf
minx -= xbuf

ybuf = (maxy - miny)/10
maxy += ybuf
miny -= ybuf 

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim(minx,maxx)
ax.set_ylim(miny,maxy)

scat = plt.scatter([],[])

import matplotlib.animation as animation
anim = animation.FuncAnimation(fig, update_plot, fargs = (scat, positions, masses),
		frames = len(positions)-1, interval = 1)

plt.show()