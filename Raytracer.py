""""
Reid Trisler 
Ray Tracing

This program renders objects (Spheres and a Plane) using a recursive ray tracing method. The ray tracing method is enhanced by
implementing 4x oversampling as well as shadows with the use of shadow feeler rays
"""
from tkinter import *
import cmath
import math
import copy
import sys

# Tkinter Canvas 
CanvasHeight = 800
CanvasWidth = 800
Ia = 0.3 # intensity of the ambient light in the scene
Ip = .7 # intensity of the point light source in the scene
# distance 
d = 500
centerofprojection = [0,0,-d]



# Sphere Object Class - instantiated with centerPoint, radius, localRGB, Kd, Ks, specIndex, localWeight, reflectWeight, refractWeight
# Contains methods that calculate where a ray will intersect with the sphere.
class Sphere:
    centerPoint = [0,0,0]
    radius = 0
    localRGB = [0,0,0]
    Kd = 0
    Ks = 0
    specIndex = 0
    localweight = 0
    reflectWeight = 0
    refractWeight = 0

    def __init__(self, centerPoint, radius, localRGB, Kd, Ks, specIndex, localWeight, reflectWeight, refractWeight):
        self.centerPoint = centerPoint
        self.radius = radius
        self.localRGB = localRGB
        self.Kd = Kd
        self.Ks = Ks
        self.specIndex = specIndex
        self.localWeight = localWeight
        self.reflectWeight = reflectWeight
        self.refract = refractWeight

    # The intesect method returns an array containing the intersection point of the ray and the sphere
    def intersect(self, startPoint, ray, currentT):
        # The ray-sphere intersection equation: 𝒂 ∙ 𝒕^𝟐 + 𝒃 ∙ 𝒕 + 𝒄 = 0
        # a = i^2 + j^2 + k^2 (i, j, and k represents the traced ray so we add each component to a)
        a = 0
        for i in range(3):
            a += ray[i] ** 2
        # 𝒃 = 𝟐 ∙ 𝒊 ∙ (𝑿𝟏 − 𝒍) + 𝟐 ∙ 𝒋 ∙ (𝒀𝟏 − 𝒎) + 𝟐 ∙ 𝒌 ∙ (𝒁𝟏 − 𝒏)
        b = 0
        for i in range(3):
            b += 2 * ray[i] * (startPoint[i] - self.centerPoint[i])
        # 𝒄 = 𝒍^𝟐 + 𝒎^𝟐 + 𝒏^𝟐 + 𝑿𝟏^𝟐 + 𝒀𝟏^𝟐 + 𝒁^𝟏𝟐 + 𝟐 ∙ (−𝒍 ∙ 𝑿𝟏 − 𝒎 ∙ 𝒀𝟏 − 𝒏 ∙ 𝒁𝟏) − 𝒓𝟐
        c = 0
        for i in range(3):
            c += self.centerPoint[i] ** 2
            c += startPoint[i] ** 2
        c2 = 0
        for i in range(3):
            c2 += -self.centerPoint[i] * startPoint[i]
        c2 *= 2
        c += c2 - self.radius ** 2
        # compute the discriminant
        discriminant = b**2 - 4 * a * c
        # If the discriminant (𝑏2 − 4 ∙ 𝑎 ∙ 𝑐) < 0 then there are no real roots (no intersection).
        if discriminant < 0:
            return sys.maxsize, 0
        #If the discriminant (𝑏2 − 4 ∙ 𝑎 ∙ 𝑐) = 0 then there is one real root (one intersection).
        #If the discriminant (𝑏2 − 4 ∙ 𝑎 ∙ 𝑐) > 0 then there are two real roots (two intersections, select nearest).
        else:
            t = (-b - math.sqrt(discriminant)) / (2 * a)
            # check for a closer object 
            if t > currentT or t < 0.001:
                return sys.maxsize, 0
            # get the intersection points
            intersection = []
            for i in range(3):
                intersection.append(startPoint[i] + (ray[i] * t))
            # check the z value 
            if intersection[2] < 0:
                return sys.maxsize, 0
            # otherwise return the calulated t and intersection
            return t, intersection

    def getlocalColor(self, intersection):
        return self.localRGB

# Checkerboard Class - Instantiated with a surface normal, acnchor point, Kd, Ks, specIndex, localWeight, reflectWeight, refractWeight
# Contains methods to get the intersection of the ray with the plane as well as getting the color of the alternating red and white board pattern   
class Checkerboard:
    normal = [0,0,0]
    anchor = [0,0,0]
    Kd = 0
    Ks = 0
    specIndex = 0
    localWeight = 0
    reflectWeight = 0
    refractWeight = 0

    def __init__(self, normal, anchor, Kd, Ks, specIndex, localWeight, reflectWeight):
        self.normal = normalize(normal)
        self.anchor = anchor
        self.Kd = Kd
        self.Ks = Ks
        self.specIndex = specIndex
        self.localWeight = localWeight
        self.reflectWeight = reflectWeight

    # the intersect with a sphere takes in a start point and a ray and returns the intersection point(if any) of the ray and the Checkerboard Plane
    def intersect(self, startPoint, ray, currentT):
        #𝑡 = −(𝐴 ∙ 𝑋1 + 𝐵 ∙ 𝑌1 + 𝐶 ∙ 𝑍1 − 𝐷) / 𝐴 ∙ 𝑖 + 𝐵 ∙ 𝑗 + 𝐶 ∙ 𝑘
        # so first we can calulate the denominator and make sure its not zero
        denominator = 0
        for i in range(3):
            denominator += self.normal[i] * ray[i]

        if denominator == 0: # no intersection
            t = sys.maxsize
            intersection = 0
            return t, intersection

        # Next compute the plane offset D
        D = 0
        for i in range(3):
            D += self.normal[i] * self.anchor[i]
        # numerator for equation
        numerator = 0
        for i in range(3):
            numerator += self.normal[i] * startPoint[i]
        # subtract D from numerator and multiply by -1
        numerator -= D
        numerator *= -1
        # calculate the value of t
        t = numerator / denominator
        # check if there is another object close by comparing the current t val with the current lowest val
        if t > currentT:
            t = sys.maxsize
            intersection = 0
            return t, intersection
        if t > 0.001:
            # Compute the point of intersection
            intersection = []
            for i in range(3):
                intersection.append(startPoint[i] + ray[i] * t)

            if (intersection[2] < 0) or (intersection[2] > 2000):
                # no visible intersection
                t = sys.maxsize 
                intersection = 0
                return t, intersection
            
            return t, intersection
        # t <= .001
        else:
            
            t = sys.maxsize
            intersection = 0
            return t, intersection
    
    # function that determins the color of the board based on the intersection point
    def getlocalColor(self, intersection):
        if intersection[0] >= 0.0:
            color_flag = True
        else:
            color_flag = False

        if abs(intersection[0]) % 400.0 > 200.0:
            color_flag = not color_flag
        
        if abs(intersection[2]) % 400.0 > 200.0:
            color_flag = not color_flag

        if color_flag:
            return [1,0,0]
        else:
            return [1,1,1]

# Driver for traceRay method. Implements 4X Oversampling. Loops through every pixel and calculates the color of each
# pixel by using the traceRay method. Then the average of the 4 pixel colors surrounding the current pixel is taken
# and that color is what is painted in the scene at that pixel
def renderImage():
    # Create array to store all the colors of the pixels
    pixelColor = [[0 for i in range(CanvasHeight + 1)] for j in range(CanvasWidth + 1)]
    top = round(CanvasHeight/2)
    bottom = round(-CanvasHeight/2)
    left = round(-CanvasWidth/2)
    right = round(CanvasWidth/2)
    # get the colors for each of the pixels and place into our pixelColor array
    for y in range(top,bottom-1, -1):
        for x in range(left, right+1):
            pixelColor[right + x][top - y] = traceRay(centerofprojection, computeUnitVector(centerofprojection, [x,y,0]), 4, 1)
    for y in range(CanvasHeight):
        for x in range(CanvasWidth):
            averageofPixels = getAverageColor(pixelColor[x][y], pixelColor[x+1][y], pixelColor[x][y+1], pixelColor[x+1][y+1])
            w.create_line(x, y, x+1, y, fill = RGBColorHexCode(averageofPixels))

""" traceRay method with adaptive depth control that returns the color for a pixel. for reflected rays only
    When tracing a ray, if the ray doesn’t hit anything, the color of the sky dome is returned. When a ray
    intersects one or more objects, the closest intersection, that is not behind the ray’s starting point, is
    determined and the algorithm computes: (1) a local color for the intersected object and  (2) the reflection vector and traces the reflected ray recursively, until some
    maximum depth is reached, returning a reflected color. The local and reflected are combined to get the final color
"""
def traceRay(startPoint, ray, depth, totalReflectWeight):
    # Stop Tracing 
    if depth == 0:
        return [0,0,0]
    
    # intersect the ray with all objects to determine nearestObject (if any)
    tMin = sys.maxsize # initialize t to a very large number
    for object in scene:
        # get the intersection point and the distance to that point
        t, intersection = object.intersect(startPoint, ray, tMin)
        # Determine if another object is closer by comparing the distances
        if t < tMin:
            tMin = t
            currentintersection = intersection
            nearestObject = object
    # return skyColor if no intersection
    if tMin == sys.maxsize: 
        return [0.53, 0.81, 0.92]

    # If we have a sphere we need to calculate the surface normal in order to apply the illumination model
    if isinstance(nearestObject, Sphere):
        nearestObjectNorm = computeUnitVector(nearestObject.centerPoint, currentintersection)
        color = nearestObject.localRGB

    # Assign the norm and color of the checkerboard        
    elif isinstance(nearestObject, Checkerboard):
        nearestObjectNorm = normalize(nearestObject.normal)
        color = nearestObject.getlocalColor(currentintersection)
        
    # Getting the ambient diffuse, and specular components of the intensity
    ADS = getphongIntensity(nearestObjectNorm, nearestObject.Kd, nearestObject.Ks, nearestObject.specIndex, currentintersection)
    # add the components together and multiply by 2
    intensity = ADS[0] + ADS[1] + ADS[2]
    intensity *= 2
    #check if object is in shadow, if so reduce intensity val
    if inShadow(nearestObject, currentintersection):
        # lower intensity value
        intensity *= 0.25
    # get local color
    localColor = [color[0] * intensity, color[1]*intensity, color[2]*intensity]
    # get weight of color
    localWeight = nearestObject.localWeight
    # get weight of reflection and multiply it by the total reflection weight
    reflectWeight = nearestObject.reflectWeight
    totalReflectWeight *= reflectWeight
    # check to see if the total weight is less than 10% and if so return [0,0,0], other wise traceRay using the intersection as the starting point
    if totalReflectWeight > 0.1:
        reflectColor = traceRay(currentintersection, reflectRay(nearestObjectNorm, ray), depth - 1, totalReflectWeight)
    else:
        reflectColor = [0,0,0]
    
    # combine the local and reflected colors together using their respective weights
    returnColor = []
    for i in range(3):
        returnColor.append(localColor[i] * localWeight + reflectColor[i] * reflectWeight)
    return returnColor   

# Returns the ambient diffuse and specular components for the provided 
def getphongIntensity(normal, Kd, Ks, specIndex, intersection):
        # ambient diffuse component of illumination model
        ambient = Ia * Kd
        L = computeUnitVector(intersection, light)
        NdotL = normal[0] * L[0] + normal[1] * L[1] + normal[2] * L[2]
        if NdotL < 0:
            NdotL = 0
        diffuse = Ip * Kd * NdotL
        # reflection vector
        R = reflectlight(normal, L) # return vector is normalized in "reflect"
        RdotV = R[0] * V[0] + R[1] * V[1] + R[2] * V[2]
        if RdotV < 0:
            RdotV = 0
        specular = Ip * Ks * RdotV ** specIndex
        return [ambient, diffuse, specular]

# Computes the reflection vector given a surface normal and a lighting vector. Used in the phong illumination model
def reflectlight(N, L):
    R = []
    # normalize surface normal 
    N = normalize(N)
    # normalize light
    L = computeUnitVector(N, light)

    twoCosPhi = 2 * (N[0]*L[0] + N[1]*L[1] + N[2]*L[2])

    #if phi == 90
    if twoCosPhi == 0:
        for i in range(3):
            R.append(-L[i])
    # if phi < 90
    elif twoCosPhi > 0:
        for i in range(3):
            R.append(N[i] - (L[i] / twoCosPhi))
    # phi > 90
    else:
        for i in range(3):
            R.append(-N[i] + (L[i] / twoCosPhi))
    return normalize(R)


# Computes the reflection vector of a traced Ray using the surface normal and a traced ray
def reflectRay(N, T):
    # The modified equation for reflecting a traced ray, 𝑇, for incident angles 𝝓 < 𝟗𝟎˚ is:𝑅 = 𝑁 + 𝑇 ∙1 / 2 ∙ cos(𝜙)
    # Or, when 𝑁 and 𝑇 are unit vectors:𝑹 = 𝑵 +𝑻𝟐 ∙ (𝑵𝑿 ∙ −𝑻𝑿 + 𝑵𝒀 ∙ −𝑻𝒀 + 𝑵𝒁 ∙ −𝑻𝒁)
    R = []
    # normalize surface normal
    N = normalize(N)
    T = normalize(T)

    # 2 ∙ cos(𝜙)
    twoCosPhi = 2 * (N[0]*-T[0] + N[1]*-T[1] + N[2]*-T[2])

    #When 𝝓 = 𝟗𝟎˚ the traced ray grazes the surface and thus,𝑹 = 𝑻
    if twoCosPhi == 0:
        for i in range(3):
            R.append(T[i])
    # if phi < 90
    elif twoCosPhi > 0:
        for i in range(3):
            R.append(N[i] + (T[i] / twoCosPhi))
    # phi > 90
    else:
        for i in range(3):
            R.append(-N[i] - (T[i] / twoCosPhi))
    return normalize(R)

# Takes in an intensity and returns a hex color code
def colorHexCode(intensity):
    hexString = str(hex(min(round(255 * intensity), 255)))
    if hexString[0] == "-":
        trimmedHexString = "00"
    else:
        trimmedHexString = hexString[2:]
        if len(trimmedHexString) == 1:
            trimmedHexString = "0" + trimmedHexString
    return trimmedHexString

# Takes in intesnity in gets hex code of each component and concats them together
def RGBColorHexCode(intensity):
    red = colorHexCode(intensity[0])
    green = colorHexCode(intensity[1])
    blue = colorHexCode(intensity[2])
    colorString = "#" + red + green + blue
    return colorString

# Takes in a start and end point and returns the unit vector
def computeUnitVector(start, end):
    return normalize([end[0]-start[0], end[1]-start[1], end[2]-start[2]])

# normalizes vectors
def normalize(vector):
    sumofSquares = 0
    for i in range(len(vector)):
        sumofSquares += vector[i] ** 2
    magnitude = math.sqrt(sumofSquares)
    vect = []
    for i in range(len(vector)):
        vect.append(vector[i]/ magnitude)
    return vect

# Determines if a point is in shadow by tracin a ray back to the light source
def inShadow(nearestObject, currentIntersection):
    shadowRay = computeUnitVector(currentIntersection, light)
    # for each object that is not the nearest object 
    for object in scene:
        if object != nearestObject:
            # check if there is an intersection
            t, intersection = nearestObject.intersect(currentIntersection, shadowRay, sys.maxsize )
            if t != sys.maxsize:
                return True
    return False

# This function takes in 4 pixels and get the average color between them in order to improve the traced image
def getAverageColor(pixel1, pixel2, pixel3, pixel4):
    average = []
    for i in range(3):
        average.append((pixel1[i] + pixel2[i] + pixel3[i] + pixel4[i]) / 4)
    return average







#########################################################################################################################
# View vector
V = [0, 0, -1] 
V = normalize(V)
# light source
light = [500, 500, 0]

# Instantiating spheres in the scene given: centerPoint, radius, localRGBcolor, Kd, Ks, specindex, localweigyt, weight for reflections
blueSphere = Sphere([0,0,900],350,[0.5,0.5,1],0.5,0.5,16,0.2,0.7,0)
redSphere = Sphere([400,-100,300], 200, [1,0.5,0.5], 0.5, 0.5, 2, 0.5, 0.5, 0)
greenSphere = Sphere([-350,-200,300], 175, [0.5,1,0.5], 0.5, 0.5, 16, 0.5, 0.5, 0)
magentaSphere = Sphere([0,-200,150], 100, [1.0,0,1.0], 0.5, 0.5, 8, 0.6, 0.4, 0)


# plane constant
PlaneConstant = -400
# Define the checkerboard plane given: surface normal, an anchor point, Kd, Ks, specIndex, weight local, weight for reflections
board = Checkerboard([0,1,0], [0, PlaneConstant,0], 0.6, 0.4, 8, 0.8, 0.25)

# The scene consists of all objects
scene = [blueSphere, redSphere, greenSphere, magentaSphere,board]

# Tkinter Canvas
root = Tk()
outerframe = Frame(root)
outerframe.pack()
w = Canvas(outerframe, width=CanvasWidth, height=CanvasHeight)
renderImage() # Call our render image function
w.pack()
root.mainloop()