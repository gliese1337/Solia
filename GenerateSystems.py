from mpmath import *
import OrbitData

class Body():
	def __init__(self,mass,x,y,vx,vy):
		self.mass = mass
		self.x=x
		self.y=y
		self.vx=vx
		self.vy=vy

		
G = mpf(6.67384e-11)
R = mpf(2591111127.56519)
M = mpf(1.59128e29)

def calcInitial(mass, dist, angle):
	x = cos(angle)*dist
	y = sin(angle)*dist
	speed = sqrt(G*M/dist)
	vy = cos(angle)*speed
	vx = -sin(angle)*speed
	return Body(mass, x, y, vx, vy)

Sun = Body(M,0,0,0,0)

for i in xrange(1,10):
	Solia = calcInitial(5.972e24,R+i*5e7/4,radians(0))
	Moon = calcInitial(5.972e24/4,R-i*5e7,radians(300))
	with open("Solia%d.in"%(i,), "w") as f:
		OrbitData.write(f, 0, [Sun, Solia, Moon])
	print "nbody.exe -o 2500 -t 1e8 < Solia%d.in > Solia%d.out" % (i,i)
	print "python GraphOrbits.py --smooth 2543000 --f Solia%d.png < Solia%d.out" % (i,i)