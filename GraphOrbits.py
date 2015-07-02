from itertools import combinations
import matplotlib.pyplot as plt
import argparse
import numpy
import math
import sys
import OrbitData

parser = argparse.ArgumentParser()
parser.add_argument('--smooth', dest='smooth', type=int, default=0)
parser.add_argument('--tscale', dest='tscale', type=float, default=1)
parser.add_argument('--start', dest='start', type=int, default=0)
parser.add_argument('--stop', dest='stop', type=int, default=0)
parser.add_argument('--f', dest='fname', type=str, default="")

args = parser.parse_args()

smooth = args.smooth
start = args.start
stop = args.stop
tscale = max(args.tscale,1)

xs = []
distances = []
for time, _, dlist in OrbitData.read(sys.stdin):
	xs.append(time)
	distances.append(dlist)

print "Read all data"

distances = numpy.transpose(distances)

print "Transposed distance arrays"

def movingaverage(interval, window_size):
	window = numpy.ones(int(window_size))/float(window_size)
	return numpy.convolve(interval, window, 'same')

if smooth > 1:
	avrgdt = sum(b-a for a,b in zip(xs,xs[1:]))/(len(xs)-1)
	window = int(round(smooth/(2*avrgdt)))
	xs = xs[window:-window]
	distances = [movingaverage(dl, window)[window:-window] for dl in distances]

starti = 0
while xs[starti] < start: starti += 1

if stop <= start:
	xs = xs[starti:]
	distances = [dl[starti:] for dl in distances]
else:
	stopi = len(xs) - 1
	while xs[stopi] > stop: stopi -= 1
	xs = xs[starti:stopi]
	distances = [dl[starti:stopi] for dl in distances]

if tscale > 1:
	xs = [t/tscale for t in xs]

rows = len(distances)
for i, dlist in enumerate(distances):
	plt.subplot(rows, 1, i+1)
	plt.plot(xs, dlist)
	ymin, ymax = min(dlist), max(dlist)
	buffer = (ymax - ymin)/10.0
	plt.ylim(ymin-buffer, ymax+buffer)
	print "Finished Plot ", (i+1)

if args.fname == "":
	plt.show()
else:
	plt.savefig(args.fname, bbox_inches='tight')