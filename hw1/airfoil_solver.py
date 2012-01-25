from numpy import *
import matplotlib.pyplot as pylab
import sys

from flow import Flow

def thicknessY(tau, xc):
    # tau = thickness/chord
    # xc = x/c
    return 5*tau*(0.2969*sqrt(xc) - 0.1260*(xc) - 0.3537*(xc)**2 + 0.2843*(xc)**3 - 0.1015*(xc)**4)

def dthicknessY(tau, xc):
    # tau = thickness/chord
    # xc = x/c
    return 5*tau*(0.2969*0.5/sqrt(xc) - 0.1260 - 2*0.3537*(xc) + 3*0.2843*(xc)**2 - 4*0.1015*(xc)**3)

def sourceU(strength, x, y):
    r2 = x*x + y*y
    return strength/(2*pi) * array([x/r2, y/r2])

class Panel():
    def __init__(self, position0, position1):
        self._position0 = array(position0)
        self._position1 = array(position1)

    def centerPos(self):
        return 0.5*(self._position0 + self._position1)

    def normal(self):
        dxy = self._position0 - self._position1
        normal = array([-dxy[1], dxy[0]])
        return normal / sqrt(dot(normal,normal))

    def length(self):
        delta = self._position0 - self._position1
        return dot(delta, delta)


class Airfoil():
    def __init__(self, tau=0.12, chord=1.0, nPanels=6):
        self.tau = tau
        self.chord = chord

        self.nPanels = nPanels
        # make number of panels even by adding one if odd
        if self.nPanels % 2 == 1:
            self.nPanels += 1
            print "setting number of panels to " + str(self.nPanels) + " to make them even"
    
        self.nChordPoints = self.nPanels/2 + 1
        xcs = linspace(0,1,self.nChordPoints)
        xcsTop      = xcs[1:-1]
        xcsBottom   = xcsTop[::-1]

        def equalArea():
            frontBunchingParam = 2.0

            sys.stdout.write("equalizing panel areas... ")
            sys.stdout.flush()
            xcs = array(xcsTop)
            badSteps = 0
            
            for k in range(300):
                xs = hstack((0.0, array(xcs), 1.0))
                ys = array([0.0]+[self.halfThickness(chord*xc) for xc in xcs]+[0.0])
                dxs = xs[1:] - xs[:-1]
                dys = ys[1:] - ys[:-1]
                deltas = sqrt( dxs*dxs + dys*dys )
                dydxs = [dthicknessY(tau,xc) for xc in xcs]

                mat = zeros((len(xcs)+1, len(xcs)))
                for row in range(len(xcs)+1):
                    x10 = dxs[row]
                    y10 = dys[row]
                    d10 = deltas[row]

                    if row != 0:
                        dy0dx = dydxs[row-1]
                        mat[row, row-1] = -(x10 + y10*dy0dx)/d10
                    if row != len(xcs):
                        dy1dx = dydxs[row]
                        mat[row,   row] =  (x10 + y10*dy1dx)/d10
                z = zeros((len(xcs), 1))
                e = eye(len(xcs))
                diff = hstack((z,e)) - hstack(((1 + frontBunchingParam/len(xcs))*e,z))
                mat = dot(diff, mat)
                rs  = dot(diff, deltas)
                step = -linalg.solve(mat, rs)

                if not isnan(sum(step)):
                    xcs += step
                    while any([xc <= 0 or xc >= 1 for xc in xcs]):
                        badSteps += 1
                        alpha = 0.5
                        xcs += (alpha - 1)*step
                        step = alpha*step
                    if sum(abs(step)) < 1e-12:
                        print "finished after " +str(k+1)+" iterations (" + str(badSteps) + " bad steps)"
                        return list(xcs)

            print "didn't converge after " +str(k+1)+" iterations (" + str(badSteps) + " bad steps)"
            return list(xcs)

        xcsTop = equalArea()[::-1]
        xcsBottom = xcsTop[::-1]

        positions = [[chord, 0]] + \
                    [[chord*xc,  self.halfThickness(chord*xc)] for xc in xcsTop] + \
                    [[0,0]] + \
                    [[chord*xc, -self.halfThickness(chord*xc)] for xc in xcsBottom]

#        positions = [[0,0]] + \
#                    [[chord*xc,  self.halfThickness(chord*xc)] for xc in xcsTop] + \
#                    [[chord, 0]] + \
#                    [[chord*xc, -self.halfThickness(chord*xc)] for xc in xcsBottom]

        self.panels = []
        for k in range(len(positions)-1):
            self.panels.append(Panel(positions[k], positions[k+1]))
        self.panels.append(Panel(positions[len(positions)-1], positions[0]))

    def halfThickness(self, x):
        if x < 0 or x > self.chord:
            return 0.0

        xc = x/self.chord
        return self.chord*thicknessY(self.tau, xc)

    def plot(self):
        xs = [p._position0[0] for p in self.panels + [self.panels[0]]]
        ys = [p._position0[1] for p in self.panels + [self.panels[0]]]
        pylab.plot(xs,ys,'-')

    def plotNormals(self):
        xcs = [p.centerPos()[0] for p in self.panels + [self.panels[0]]]
        ycs = [p.centerPos()[1] for p in self.panels + [self.panels[0]]]
        ucs = [p.normal()[0] for p in self.panels + [self.panels[0]]]
        vcs = [p.normal()[1] for p in self.panels + [self.panels[0]]]
        pylab.quiver(xcs,ycs,ucs,vcs,scale=15)

normUinf = 10.0
alpha = 0.0#*pi/180
rho = 1.2
uinf = array( [cos(alpha)*normUinf, sin(alpha)*normUinf] )

foil = Airfoil(tau=0.12, chord=1.0, nPanels=200)
flow = Flow(uinf, foil.panels)

flow.solveGammas(foil.panels)

# calculate forces by "integrating" pressure
forces = flow.force(foil.panels, rho)
forcesWind = dot(matrix([[cos(alpha), sin(alpha)],[-sin(alpha), cos(alpha)]]), forces)
print ""
print "forces calculated from pressure:"
print "body frame force:              " + str(forces)
print "wind frame force (drag, lift): " + str(forcesWind)
print "CL: " + str(forcesWind[0,1]/(0.5*rho*normUinf*normUinf))
print "CD: " + str(forcesWind[0,0]/(0.5*rho*normUinf*normUinf))

# calculate forces from "vorticity" pressure
print ""
print "forces calculated from vorticity:"
vorticity = sum([v._gamma for v in flow.vortices])
print "vorticity: " +str(vorticity)
vortForce = rho*cross(array( [cos(alpha)*normUinf, sin(alpha)*normUinf, 0] ), array([ 0, 0, vorticity]))
vortForceWind = dot(matrix([[cos(alpha), sin(alpha),0],[-sin(alpha), cos(alpha),0],[0,0,1]]), vortForce)
print "body frame force:                " + str(vortForce)
print "wind frame force (drag, lift,_): " + str(vortForceWind)
print "CL: " + str(vortForceWind[0,1]/(0.5*rho*normUinf*normUinf))
print "CD: " + str(vortForceWind[0,0]/(0.5*rho*normUinf*normUinf))

#flow.plotPrimitives()
foil.plot()
#foil.plotNormals()
#flow.plot(lambda x: foil.halfThickness(x), xRange=linspace(-0.05,0.05,51), yRange=linspace(-0.05,0.05,51))
#flow.plot(lambda x: foil.halfThickness(x), xRange=linspace(-0.2,1.2,31), yRange=linspace(-0.3,0.3,31))
#flow.plot(lambda x: foil.halfThickness(x))
#flow.plotForces(foil.panels, rho)
#flow.plotCps(foil.panels)
flow.plotStreamlines()
#flow.plotSurfaceVelocities(foil.panels)

pylab.axis('equal')
pylab.show()
