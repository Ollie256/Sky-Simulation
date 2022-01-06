import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

g=0.76
r0 = 6365.362
HR = 8
HM = 1.2
ratmo = 6465.362
Isun = 1300

betaR = [33e-3, 13e-3, 5e-3]
betaM = 2.1e-3

rayOrigin = [0, 0, r0]

def sunDirf(sunAzimuth, sunElevation):
    sunDir = [np.cos(sunElevation*np.pi/180)*np.cos(sunAzimuth*np.pi/180), np.cos(sunElevation*np.pi/180)*np.sin(sunAzimuth*np.pi/180), np.sin(sunElevation*np.pi/180)]
    return sunDir

def rayDirf(rayAzimuth, rayElevation):
    rayDir = [np.cos(rayElevation*np.pi/180)*np.cos(rayAzimuth*np.pi/180), np.cos(rayElevation*np.pi/180)*np.sin(rayAzimuth*np.pi/180), np.sin(rayElevation*np.pi/180)]   
    return rayDir

#rayDir = rayDirf(0, 20)
#sunDir = sunDirf(0, 20)

numOpticalDepth = 8
numView = 16

intensityRed = []
intensityGreen = []
intensityBlue = []

I = []
H = []
S = []

rgb = []
hexrgb = []



def phaseR(mu):
    return 3*(1+mu**2)/(16*np.pi)

def phaseM(mu):
    return 3*(1-g**2)*(1+mu**2)/((8*np.pi * (2+g**2) * (1 + g**2 -  2*g*mu)**1.5))

def opticalDepth(rayOrigin, rayDir, rayLength):
    samplePoint = rayOrigin
    ds = rayLength / (numOpticalDepth-1)
    opticalDepthR = 0
    opticalDepthM = 0
    for i in range(numOpticalDepth):
        samplePoint = rayOrigin + (rayDir * ds * i)
        h = np.linalg.norm(samplePoint) - r0
        if h < 0:
            break
        opticalDepthR += np.exp(-h/HR) * ds
        opticalDepthM += np.exp(-h/HM) * ds
        
    return opticalDepthR, opticalDepthM

def rayLength(sphereRadius, rayOrigin, rayDir):
    b = 2*np.dot(rayOrigin, rayDir)
    c = np.dot(rayOrigin, rayOrigin) - sphereRadius**2
    d = b**2  - 4*c
    if d >= 0:
        t = (-b + np.sqrt(d))/2
    else:
        t = -b/2
    return t

def Hf(R, G, B):
    if G >= B:
        H = np.arccos((R - G + R - B)/(2 * np.sqrt((R - G)**2 + (R - B)*(G - B))))
        H = H*180/np.pi
    else:
        H = 360 - np.arccos((R - G + R - B)/(2 * np.sqrt((R - G)**2 + (R - B)*(G - B))))*180/np.pi

def sampleLight(rayOrigin, rayDir, viewRayLength):
    samplePoint = np.array(rayOrigin)
    ds = viewRayLength / (numView-1)
    sumR = [0, 0, 0]
    sumM = [0, 0, 0]
    for i in range(numView):
        
        samplePoint = rayOrigin + (np.array(rayDir) * ds * i)
        
        mu = np.dot(rayDir, sunDir)
        
        lightLength = rayLength(ratmo, samplePoint, np.array(sunDir))
        
        lightOpticalDepthR = opticalDepth(samplePoint, np.array(sunDir), lightLength)[0]
        viewOpticalDepthR = opticalDepth(samplePoint, -np.array(rayDir), ds*i)[0]
        
        lightOpticalDepthM = opticalDepth(samplePoint, np.array(sunDir), lightLength)[1]
        viewOpticalDepthM = opticalDepth(samplePoint, -np.array(rayDir), ds*i)[1]
        #print(lightOpticalDepthR + lightOpticalDepthM, i)
        tau = np.array(betaR)*(lightOpticalDepthR + viewOpticalDepthR) + 1.1*betaM*(lightOpticalDepthM + viewOpticalDepthM)
        
        T = np.array([np.exp(-tau[0]), np.exp(-tau[1]), np.exp(-tau[2])])
        
        
        sumR += np.array(T * np.exp(-(np.linalg.norm(samplePoint) - r0)/ HR) * ds)
        sumM += np.array(T * np.exp(-(np.linalg.norm(samplePoint) - r0)/ HM) * ds)
        
        print((((sumR * betaR * phaseR(mu)) + (sumM * betaM * phaseM(mu)))[0] + ((sumR * betaR * phaseR(mu)) + (sumM * betaM * phaseM(mu)))[1] + ((sumR * betaR * phaseR(mu)) + (sumM * betaM * phaseM(mu)))[2])*Isun, mu, rayDir[2]*180/np.pi)
        
    return ((sumR * betaR * phaseR(mu)) + (sumM * betaM * phaseM(mu)))#*Isun


#fig = plt.figure()
#print(sampleLight(rayOrigin, rayDir, rayLength(ratmo, rayOrigin, rayDir)))
#print(rayLength(ratmo, rayOrigin, rayDir))   
#ims = []
for k in range(70, 71, 1):
    intensityTotal = []
    
    azimuthPlot = []
    elevationPlot = []
    azimuthElevationPlot = [] 
    for i in range(68, 73, 1):
        for j in range(0, 1, 2):
            rayDir = rayDirf(j, i)
            sunDir = sunDirf(0, k)
            
            
            
            blue = sampleLight(rayOrigin, rayDir, rayLength(ratmo, rayOrigin, rayDir))[0] #/ Isun
            green = sampleLight(rayOrigin, rayDir, rayLength(ratmo, rayOrigin, rayDir))[1]# / Isun
            red = sampleLight(rayOrigin, rayDir, rayLength(ratmo, rayOrigin, rayDir))[2] #/ Isun
            
            normalise = red+green+blue
            #print(normalise, i)
            #print(red, green, blue)
            rgb.append([0, 0, 1])
            
            #print(red/normalise, green/normalise, blue/normalise)
            #print((red+green+blue)*normalise)
            #print([red/255, green/255, blue/255])
            #hexrgb.append('#%02x%02x%02x' % (int(red), int(green), int(blue)))
            #total = (blue + green + red)*normalise
            #I.append((blue + green + red)/3)
            #H.append(Hf(red, green, blue))
            #S.append(1 - 3 * min(red, green, blue)/(red + green + blue))
            #intensityBlue.append(blue)
            #intensityGreen.append(green)
            #intensityRed.append(red)
            #intensityTotal.append(total)
            azimuthPlot.append(j)
            elevationPlot.append(i)
            #azimuthElevationPlot.append((j, i))
        #print(red, green, blue, i)
    #normalise = max(rgb)
    rgb = rgb/np.amax(rgb)
    
    #rgbN = rgb/normalise
    #norm = mpl.colors.Normalize(vmin=0, vmax=max(intensityTotal))
    plt.scatter(azimuthPlot, elevationPlot, c=rgb)#, norm=norm,)
    #plt.colorbar()
    #plt.scatter(0, k, c=1100, s=50)
    plt.show()
    plt.savefig('Sunset2_' + str(k) + '.png', dpi=500)        
        

#cmap = mpl.cm.coolwarm
#norm = mpl.colors.Normalize(vmin=0, vmax=Isun)

#print(max(intensityTotal))  
#plt.plot(intensityRed, elevationPlot, color='red')   
#plt.plot(intensityGreen, elevationPlot, color='green')  
#plt.plot(intensityBlue, elevationPlot, color='blue')  
#plt.scatter(azimuthPlot, elevationPlot, c=intensityTotal)#, cmap='ocean')   
#plt.savefig('Sky.png', dpi=500) #plt.colorbar()
"""
ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True,
                                repeat_delay=500)

writer = PillowWriter(fps=20)
ani.save("demo2.gif", writer=writer)

plt.show()
"""

print(rayLength(ratmo, [0, 0, r0], rayDirf(0, -1)))