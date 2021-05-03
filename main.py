import math
import random
from Input_file import Mc_variables
mc_v = Mc_variables()
"""Finte model of Monte carlo simulation of light transport in tissue"""

#initialization of photons
for i in range(mc_v.photons):
  #collimated launch
  x, y ,z, u, v, w, th = 0, 0, 0, 0, 0, 1, 0 
  wt=1.0-mc_v.rs
  while wt>0.0: 
    a=random.random()
    d=-math.log(a)/(mc_v.ua+ mc_v.us)
    x+=d*u
    y+=d*v
    z+=d*w
    R = math.sqrt((x * x) + (y * y))
    #radial index
    if R >= mc_v.Rmax: 
      ir = 19
    else:
      ir = (int) (R / mc_v.dr)
    #lower boundery check
    if (z < 0.0): 
      th = mc_v.pi - th
      z = abs(z)
      w = -w
      if (abs(th) <= mc_v.thc):
        c = math.sin(th)
        th1 = math.asin(mc_v.n * c)
        b = math.cos(th)
        e = math.cos(th1)
        temp = (mc_v.n * b - 1.0 * e) / (mc_v.n * b + 1.0 * e)
        temp1 = (1.0 * b - mc_v.n * e) / (1.0 * b + mc_v.n * e)
        #fresnel reflection
        rf = (c * c) * ((temp * temp) + (temp1 * temp1)) / 2.0 
        mc_v.rb += ((1.0 - rf) * wt)
        wt = rf * wt
        mc_v.heatrb[ir] += (1.0 - rf) * wt
      #upper boundery check
      elif (z > mc_v.Hmax): 
        z = 2 * mc_v.Hmax - z
        w = -w

      #check for Total Internal Reflection
      if (abs(th) <= mc_v.thc): 
        c = math.sin(th)
        th1 = math.asin(mc_v.n * c)
        b = math.cos(th)
        e = math.cos(th1)
        temp = (mc_v.n * b - 1.0 * e) / (mc_v.n * b + 1.0 * e)
        temp1 = (1.0 * b - mc_v.n * e) / (1.0 * b + mc_v.n * e)
        rf = (c * c) * ((temp * temp) + (temp1 * temp1)) / 2.0
        mc_v.rt += ((1.0 - rf) * wt)
        wt = rf * wt
        mc_v.heatrt[ir] += (1.0 - rf) * wt

    #height index
    if (z >= mc_v.Hmax): 
      iz = 99
    else:
      iz = (int) (z / mc_v.dz) 

    #Drop
    mc_v.heata[iz][ir] += (1.0 - mc_v.alb) * wt
    wt = wt * mc_v.alb

    #roulette method
    if (wt < 0.001):
      if (random.random() > 0.1):
        wt = 0.0
      else:
        wt = wt * 10.0;

    #Spin
    phi = 2.0 * mc_v.pi * random.random()
    th = math.acos(2.0 * random.random() - 1.0)
    u = math.sin(th) * math.cos(phi)
    v = math.sin(th) * math.sin(phi)
    w = math.cos(th);

#Total absorption
for i in range(100): 
  for j in range(20):
    mc_v.absv += mc_v.heata[i][j]

#output
print("Heat Absorbed : {}".format(mc_v.absv/mc_v.photons))
print("Specular Reflection = {}".format(mc_v.rs))
print("Backscattered Reflection = {}".format(mc_v.rb/mc_v.photons))
print("Transmitted Reflection = {}".format(mc_v.rt / mc_v.photons))
print("Total Reflection = {}".format(mc_v.rs + ((mc_v.rb + mc_v.rt) / mc_v.photons)))
