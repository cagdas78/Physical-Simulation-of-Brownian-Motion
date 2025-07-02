#Brownian Motion written by Cagdas Allahverdi (from Toros University)
import turtle
from turtle import *
import math
import random
import time
import matplotlib.pyplot as plt

turtle.title("Brownian Motion (BM) for Everyone by Cagdas Allahverdi")
print("Brownian Motion (BM): The Beginning")

#Dimensions of Rectangle
width=450
height=450
expand=1.1
setup(width*expand,height*expand,0,0)

colormode(255)
ht()
speed(0)

#Number of Particles (Includes the Brownian particle)
particleNumber=14

#Draw the Rectangle (In other words, the Box having four Walls)
def boxDrawing():
    t1=turtle.Turtle()
    t1.ht()
    t1.speed(0)
    t1.penup()
    t1.goto(0, -height/2)
    t1.pendown()
    t1.seth(0)
    t1.forward(width/2)
    t1.left(90)
    t1.forward(height)
    t1.left(90)
    t1.forward(width)
    t1.left(90)
    t1.forward(height)
    t1.left(90)
    t1.forward(width/2)

#Draw the Particle
def particleDrawing(x,y,r,g,b,radius):
    penup()
    goto(x,y)
    pendown()
    #fillcolor((r,g,b)) #Paint the Particle
    #begin_fill()
    circle(radius)
    #end_fill()

#Move the Particle and Bounce off the Walls
def particleMovement(x,y,theta,step,radius,index,stop,wl):

    x=x+step*math.cos(theta) #Move the particle along x-axis
    y=y+step*math.sin(theta) #Move the particle along y-axis

    if (x>=width/2-radius): #Particle bouncing off the right wall (Reflection from the right wall)
        wl=1 #Collision with the wall
        if index==0: #Is it Brownian particle?
            stop=1
        if (theta==0) or (theta==2*math.pi):
            theta=math.pi
        elif (theta<math.pi/2):
            theta=math.pi-theta
        else:
            theta=3*math.pi-theta
    elif (x<=-width/2+radius): #Particle bouncing off the left wall
        wl=1
        if index==0:
            stop=1
        if (theta==math.pi):
            theta=0
        elif (theta<math.pi):
            theta=math.pi-theta
        else:
            theta=3*math.pi-theta
    if (y>=height/2-2*radius): #Particle bouncing off the top wall
        wl=1
        if index==0:
            stop=1
        if (theta==math.pi/2):
            theta=3*math.pi/2
        elif (theta<math.pi/2):
            theta=2*math.pi-theta
        else:
            theta=2*math.pi-theta
    elif (y<=-height/2): #Particle bouncing off the bottom wall
        wl=1
        if index==0:
            stop=1
        if (theta==3*math.pi/2):
            theta=math.pi/2
        elif (theta<1.50*math.pi):
            theta=2*math.pi-theta
        else:
            theta=2*math.pi-theta

    if (x>width/2-radius): #Particle bouncing off the right wall (Reflection from the right wall)
        x=x-2*(x-width/2+radius)
    if (x<-width/2+radius): #Particle bouncing off the left wall
        x=x+2*(-width/2-x+radius)
    if (y>height/2-2*radius): #Particle bouncing off the top wall
        y=y-2*(y-height/2+2*radius)
    if (y<-height/2): #Particle bouncing off the bottom wall
        y=y+2*(-height/2-y)

    return x, y, theta, step, radius, index, stop, wl

#X'-Component of Velocity
def velocityXPrimeBefore(v,theta,alpha):
    velocityXPB=v*(math.cos(theta)*math.cos(alpha)+math.sin(theta)*math.sin(alpha)) #Tangential velocity; vx'=v*cos(theta-alpha)
    return velocityXPB

#Y'-Component of Velocity
def velocityYPrimeBefore(v,theta,alpha):
    velocityYPB=v*(math.sin(theta)*math.cos(alpha)-math.cos(theta)*math.sin(alpha)) #Normal velocity; vy'=v*sin(theta-alpha)
    return velocityYPB

#Calculate Scattering Angles of Two Particles (After the collision of two particles)
def scatteringAngle(v1T,v1N,v2T,v2N,b1,b2,alpha):
    if (v1T>0) and (v1N>0):
        angle1=b1
    elif (v1T>0) and (v1N<0):
        angle1=2*math.pi-b1
    elif (v1T<0) and (v1N>0):
        angle1=math.pi-b1
    elif (v1T<0) and (v1N<0):
        angle1=math.pi+b1
    elif (v1T==0) and (v1N>0):
        angle1=math.pi/2
    elif (v1T==0) and (v1N<0):
        angle1=3*math.pi/2
    elif (v1T>0) and (v1N==0):
        angle1=0
    else:
        angle1=math.pi

    if (v2T>0) and (v2N>0):
        angle2=b2
    elif (v2T>0) and (v2N<0):
        angle2=2*math.pi-b2
    elif (v2T<0) and (v2N>0):
        angle2=math.pi-b2
    elif (v2T<0) and (v2N<0):
        angle2=math.pi+b2
    elif (v2T==0) and (v2N>0):
        angle2=math.pi/2
    elif (v2T==0) and (v2N<0):
        angle2=3*math.pi/2
    elif (v2T>0) and (v2N==0):
        angle2=0
    else:
        angle2=math.pi

    angle1=angle1+alpha
    angle2=angle2+alpha

    if (angle1<0):
        angle1=2*math.pi+angle1
    if (angle2<0):
        angle2=2*math.pi+angle2
    if (angle1>2*math.pi):
        angle1=angle1-2*math.pi
    if (angle2>2*math.pi):
        angle2=angle2-2*math.pi

    return angle1, angle2

xPosition=[]
yPosition=[]
theta=[]
velocityBefore=[]
mass=[]
radius=[]

particlePath=[]
particleStep=[]

xPositionBrownianParticleCollision=[0]
yPositionBrownianParticleCollision=[0]
tTimeBrownianParticleCollision=[0]

#Generate Particles
for n in range(particleNumber):
    col=[random.uniform(0, 2*math.pi)]
    theta.append(col)
    col=[random.uniform(0.8, 0.8)]
    velocityBefore.append(col)
    col=random.uniform(1,1)
    mass.append(col)
    col=random.uniform(0.3, 0.3)
    radius.append(col)
    col=[0]
    particleStep.append(col)

rM=max(radius)

#Describe Bigger Particle (Brownian Particle)
mass[0]=pow(10,3) #ITS MASS
radius[0]=181 #ITS RADIUS
velocityBefore[0]=[0.0001] #ITS VELOCITY
theta[0]=[0]

#Brownian Particle Tracing
xPositionTracker=[0]
yPositionTracker=[-1*radius[0]]
velocityBeforeTracker=[0.0001] # Give velocityBefore[0]'s value
distanceToBrownian=[]

collisionCounter=0 #The number of all collisions
particleCollisionCounter=0 #Particle-particle collisions (Including Brownian particle-particle collisions)
brownianParticleCollisionCounter=0 #Brownian particle-particle collisions
brownianParticleTotalDistance=0 #Total distance taken by Brownian particle

t2 = turtle.Turtle()
t2.ht()
#Track the Brownian Particle (Draw its Trajectory)
def particleTracker(x1,y1,x2,y2):
    t2.penup()
    t2.goto(x1,y1+radius[0]) #Brownian particle initial-x and initial-y position
    t2.pendown()
    t2.goto(x2,y2+radius[0]) #Brownian particle final-x and final-y position

counter=0

#Place Particles in The Box (Randomly)
while True:
    col1=[random.uniform(-width/2+4.00*rM, width/2-4.00*rM)]
    col2=[random.uniform(-height/2+4.00*rM, height/2-4.00*rM)]
    if math.sqrt(col1[-1]**2+col2[-1]**2) > (radius[0]+4.00*rM):
        xPosition.append(col1)
        yPosition.append(col2)
        counter=counter+1
        if counter==particleNumber:
            break

xPosition[0]=[0]
yPosition[0]=[-radius[0]]
energy=[]

i=0

#Draw the Box
boxDrawing()

rank=0
stop=0
finish=0

start=time.time()
tau=0

while True:

    for j in range(particleNumber):

        if j==0:
            r, g, b = 255, 127, 127 #Color of the Brownian particle
        else:
            r, g, b = 86, 86, 86 #Color of the other particles

        distanceToBrownian.clear()

        for p in range(1,particleNumber):
            distanceToBrownian.append(math.sqrt((xPosition[0][-1]-xPosition[p][-1])**2+(yPosition[0][-1]+radius[0]-yPosition[p][-1]-radius[p])** 2)-radius[0]-radius[p])

        timeInterval = 1.0 #0.9 for watching # Time Interval (ps)

        #if j==0: # Put "#" to see all particles on the screen
        particleDrawing(xPosition[j][i], yPosition[j][i], r, g, b, radius[j]) #Show the particles
        if j==0 and i!=0:
            particleTracker(xPosition[j][i-1], yPosition[j][i-1], xPosition[j][i], yPosition[j][i])
            xPositionTracker.append(xPosition[j][i])
            yPositionTracker.append(yPosition[j][i])
            velocityBeforeTracker.append(velocityBefore[j][i])
            brownianParticleTotalDistance=brownianParticleTotalDistance+math.sqrt((yPosition[j][i]-yPosition[j][i-1])**2+(xPosition[j][i]-xPosition[j][i-1])**2)

        step=velocityBefore[j][i]*timeInterval #Displacement of the particle
        particleStep[j].append(step)
        wall=0 #No collision with the wall
        newLocation=particleMovement(xPosition[j][i], yPosition[j][i], theta[j][i], step, radius[j], j, stop, wall)
        if newLocation[7]==1: #To Check Particle-Wall Collision.
            collisionCounter=collisionCounter+1
            particlePath.append(sum(particleStep[j]))
            particleMeanFreePath = sum(particlePath) / len(particlePath)
            particleStep[j] = [0]
        xPosition[j].append(newLocation[0])
        yPosition[j].append(newLocation[1])
        theta[j].append(newLocation[2])
        rank=newLocation[5]
        finish=newLocation[6]
        if (rank==0) and (finish==1): #If the Brownian particle hits the wall
            break #Stop movement
        velocityBefore[j].append(velocityBefore[j][i])
        energy.append(0.5*mass[j]*velocityBefore[j][i]**2) #Kinetic energy of particle

        #Check for Collision; if any, calculate next parameters
        for k in range(j + 1, particleNumber):
            minDistance=radius[j]+radius[k]
            if (xPosition[j][-1]-xPosition[k][-1])**2+(yPosition[j][-1]+radius[j]-yPosition[k][-1]-radius[k])**2 <= minDistance ** 2:
                collisionCounter=collisionCounter+1 #Particle-Particle collision occurred.
                particleCollisionCounter=particleCollisionCounter+1 #Number of particle-particle collisions (Including Brownian particle-particle collision)
                if j==0: #If True, it is Collision of the Brownian Particle
                    brownianParticleCollisionCounter=brownianParticleCollisionCounter+1 #Number of Brownian particle-particle collisions
                    xPositionBrownianParticleCollision.append(xPosition[j][-1])
                    yPositionBrownianParticleCollision.append(yPosition[j][-1]+radius[0])
                    tau=tau+len(particleStep[j])*timeInterval
                    tTimeBrownianParticleCollision.append(tau)
                m1=(yPosition[j][-1]+radius[j]-yPosition[k][-1]-radius[k])/(xPosition[j][-1]-xPosition[k][-1]) #The slope of the line passing through the centers of the colliding particles. This line is called Y'.
                m2=-1/m1 #The slope of the line tangent to the point of contact of the colliding particles. This line is called X'.
                alpha=math.atan(m2) #Angle between X and X' axes. As a note, X-Y is the stationary Cartesian coordinate system.
                vjxp=velocityXPrimeBefore(velocityBefore[j][i],theta[j][i],alpha)
                vkxp=velocityXPrimeBefore(velocityBefore[k][i],theta[k][i],alpha)
                vjyp=velocityYPrimeBefore(velocityBefore[j][i],theta[j][i],alpha)
                vkyp=velocityYPrimeBefore(velocityBefore[k][i],theta[k][i],alpha)
                vjTangentialBefore=vjxp #Tangential velocity of particle j before the collision at X'-Y' coordinate system
                vkTangentialBefore=vkxp #Tangential velocity of particle k before the collision at X'-Y' coordinate system
                vjNormalBefore=vjyp #Normal velocity of particle j before the collision at X'-Y' coordinate system
                vkNormalBefore=vkyp #Normal velocity of particle k before the collision at X'-Y' coordinate system

                vjTangentialAfter = vjTangentialBefore #Tangential velocity of particle j after the collision at X'-Y' coordinate system
                vkTangentialAfter = vkTangentialBefore #Tangential velocity of particle k after the collision at X'-Y' coordinate system
                vjNormalAfter = ((mass[j]-mass[k])*vjNormalBefore+2*mass[k]*vkNormalBefore)/(mass[j]+mass[k]) #Normal velocity of particle j after the collision at X'-Y' coordinate system
                vkNormalAfter = ((mass[k]-mass[j])*vkNormalBefore+2*mass[j]*vjNormalBefore)/(mass[j]+mass[k]) #Normal velocity of particle k after the collision at X'-Y' coordinate system
                velocityBefore[j][-1]=math.sqrt(vjTangentialAfter**2+vjNormalAfter**2) #Velocity of particle j after the collision
                velocityBefore[k][-1]=math.sqrt(vkTangentialAfter**2+vkNormalAfter**2) #Velocity of particle k after the collision

                betaj = math.atan(abs(vjNormalAfter / vjTangentialAfter)) #Beta for particle j after the collision (Beta is the scattering angle of the particle at X'-Y' coordinate system. It is between 0-90 degrees)
                betak = math.atan(abs(vkNormalAfter / vkTangentialAfter)) #Beta for particle k after the collision (Beta is the scattering angle of the particle at X'-Y' coordinate system. It is between 0-90 degrees)

                newTheta=scatteringAngle(vjTangentialAfter,vjNormalAfter,vkTangentialAfter,vkNormalAfter,betaj,betak,alpha) #Scattering angles of two collided particles
                theta[j][-1]=newTheta[0] #Scattering angle of particle j with respect to X-Y coordinate system
                theta[k][-1]=newTheta[1] #Scattering angle of particle k with respect to X-Y coordinate system

                if j!=0:
                    particlePath.append(sum(particleStep[j]))
                if k!=0:
                    particlePath.append(sum(particleStep[k]))
                particleMeanFreePath=sum(particlePath)/len(particlePath)
                #print("Mean Free Path= " + str(particleMeanFreePath))
                particleStep[j]=[]
                particleStep[k]=[]
    i=i+1

    if (rank == 0) and (finish == 1):
        plt.hist(energy) #Plot energy distribution
        plt.title('Energy Distribution')
        plt.xlabel('Energy')
        plt.ylabel('Particle Number')
        #plt.pause(0.05)
        #plt.clf()
        #print("Total Energy= "+str(sum(energy)))

    time.sleep(0.05) #1/24
    turtle.tracer(1, 0)
    clear()
    energy.clear()

    if i==2:
        for j in range(particleNumber):
            xPosition[j].pop(0)
            yPosition[j].pop(0)
            theta[j].pop(0)
            velocityBefore[j].pop(0)
        i=i-1

    workingTime=(time.time()-start)/60.0

    if (rank==0) and (finish==1): #Other choice: if workingTime>=10.0:
        filename = "BrownianMotion16112024a14102u"
        f = open(filename + ".txt", "w")
        f.write("Title=" + filename + "\n")
        f.write("Running time in minutes(min)=" + str(workingTime) + "\n")
        f.write("Number of all particles=" + str(particleNumber) + "\n")
        f.write("Number of all (Brownian particle-particle, particle-particle, Brownian particle-wall(Once) and particle-wall) collisions=" + str(collisionCounter) + "\n")
        f.write("Number of particle-particle (excluding Brownian particle) collisions=" + str(particleCollisionCounter-brownianParticleCollisionCounter) + "\n")
        f.write("Number of Brownian particle-particle collisions=" + str(brownianParticleCollisionCounter) + "\n")
        f.write("Number of particle-wall (excluding Brownian particle-wall) collisions=" + str(collisionCounter-particleCollisionCounter-1) + "\n")
        f.write("Number of Brownian particle-wall collisions=" + str(1) + "\n") #Once for Stop
        f.write("Total distance taken by Brownian particle=" + str(brownianParticleTotalDistance)+ " px" + "\n")
        f.write("Mean free path of particles (except Brownian particle)=" + str(particleMeanFreePath) + " px" + "\n")
        f.write("Radius of small particles=" + str(radius[1]) + " px" + "\n")
        f.write("Mass of small particles=" + str(mass[1]) + " kg" + "\n")
        f.write("Radius of Brownian (bigger) particle=" + str(radius[0]) + " px" + "\n")
        f.write("Mass of Brownian particle=" + str(mass[0]) + " kg"+"\n")
        f.write("Final velocity of Brownian particle=" + str(velocityBefore[0][-1]) + "\n")
        f.write("Final energy of Brownian particle=" + str(0.5*mass[0]*velocityBefore[0][-1]**2) + "\n")
        f.write("Time interval=" + str(timeInterval) + "\n")
        f.write("Brownian Motion (BM):" + "\n")
        f.write("x coordinate" + "\t" + "y coordinate" + "\t" + "velocity" + "\n")
        for p in range(len(xPositionTracker)):
            f.write(str(xPositionTracker[p]) + "\t" + str(yPositionTracker[p]+radius[0]) + "\t" + str(velocityBeforeTracker[p]) + "\n")
        f.close()
        screen = turtle.getscreen()
        canvas = screen.getcanvas()
        canvas.postscript(file=filename + ".eps", width=width*expand, height=height*expand)

        filename = "BrownianMotion16112024a14102v"
        f = open(filename + ".txt", "w")
        f.write("Title=" + filename + "\n")
        f.write("Collision Trajectory of Brownian Particle" + "\n")
        f.write("Collision Time/Tau (a.u.)" + "\t" + "x-Position of Collision (px)" + "\t" + "y-Position of Collision (px)" + "\n")
        for p in range(len(tTimeBrownianParticleCollision)):
            f.write(str(tTimeBrownianParticleCollision[p]) + "\t" + str(xPositionBrownianParticleCollision[p]) + "\t" + str(yPositionBrownianParticleCollision[p]) + "\n")
        f.close()

        print("Brownian Motion (BM): The End")
        break
done()