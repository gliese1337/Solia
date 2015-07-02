from mpmath import *
import OrbitData
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--step', dest='step', type=int, default=1)
parser.add_argument('--out', dest='out', type=int, default=1)
parser.add_argument('--end', dest='end', type=int)

args = parser.parse_args()

G = mpf(6.67384e-11)

class Body():
	def __init__(self,mass,x,y,vx,vy):
		self.mass = mpf(mass)
		self.x = mpf(x)
		self.y = mpf(y)
		self.vx = mpf(vx)
		self.vy = mpf(vy)
		self.ax = 0
		self.ay = 0

def step(bodies, dt, first):
	vdt = dt * (G/2.0 if first else G)
	distances = []
	for i, p in enumerate(bodies):
		px, py = p.x, p.y
		for b in bodies[i+1:]:
			dx = px-b.x;
			dy = py-b.y
			h = hypot(dx, dy)
			h3 = h**3

			#G is accounted for in vdt
			p.ax -= b.mass * dx / h3
			p.ay -= b.mass * dy / h3

			b.ax += p.mass * dx / h3
			b.ay += p.mass * dy / h3

		#Update velocity & position
		p.vx += p.ax*vdt
		p.vy += p.ay*vdt
		p.x += p.vx*dt
		p.y += p.vy*dt
		p.ax = 0
		p.ay = 0

reader = OrbitData.read(sys.stdin)
time, b_info, _ = reader.next()

minmass = min(map(lambda t: t[0], b_info))
G = G*minmass

bodies = [Body(data[0]/minmass, *data[1:]) for data in b_info]

timestep = args.step
output_interval = args.out
end = args.end

step(bodies, timestep, True)
for s in xrange(time, end, timestep):
	step(bodies, timestep, False)
	if s % output_interval == 0:
		OrbitData.write(sys.stdout, s, bodies)