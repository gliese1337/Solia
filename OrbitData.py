from itertools import combinations
from mpmath import hypot

def read(input):
	count, time = 0, 0
	bodies = None
	state = "num"
	for line in input:
		line = line.strip()
		if line == "": continue
		if state == "num":
			count, time = map(float, line.split(' '))
			bodies = []
			state = "body"
		elif state == "body":
			mass, x, y, vx, vy = map(float, line.split(' '))
			bodies.append((mass, x, y, vx, vy))
			count -= 1
			if count == 0:
				state = "dist"
		elif state == "dist":
			distances = map(float, line.split(' '))
			state = "num"
			yield time, bodies, distances

def write(output, time, bodies):
	output.write("%d %f\n" % (len(bodies), time))
	for b in bodies:
		output.write("%f %f %f %f %f\n" % (b.mass, b.x, b.y, b.vx, b.vy))
	distances = []
	for p, b in combinations(bodies, 2):
		distances.append(str(hypot(p.x-b.x, p.y-b.y)))
	output.write(" ".join(distances)+'\n')