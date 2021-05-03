import math
import random
from Input_file import Mc_variables

mc_v = Mc_variables()
"""Semi infinte model of Monte carlo simulation of light transport in tissue"""
heata = 0.0
#Initialization of photons
for i in range(mc_v.photons):
    #Collimated launch
  x, y, z, u, v, w, th = 0, 0, 0, 0, 0, 1, 0
  wt = 1.0 - mc_v.rs
  while wt>0.0:
    #print(wt)
    a = random.random()
      #Step size
    d = -math.log(a) / (mc_v.ua + mc_v.us)
    x += d * u
    y += d * v
    z += d * w
      #boundary check
    if z < 0:
      th = mc_v.pi - th
      z = -z
      w = -w
      if abs(th) <= mc_v.thc:
        c = math.sin(th)
        th1 = math.asin(mc_v.n * c)
        b = math.cos(th)
        e = math.cos(th1)
        temp = (mc_v.n * b - 1.0 * e) / (mc_v.n * b + 1.0 * e)
        temp1 = (1.0 * b - mc_v.n * e) / (1.0 * b + mc_v.n * e)
          #Fresnel reflection
        rf = (c * c) * ((temp * temp) + (temp1 * temp1)) / 2.0
        mc_v.rb += ((1.0 - rf) * wt)
        wt = rf * wt

        #Drop
    heata += (1.0 - mc_v.alb) * wt
    wt = wt * mc_v.alb
        #Roulette method
    if wt < 0.001:
      if random.random() > 0.1:
        wt = 0.0
      else:
        wt = wt * 10.0

        #Spin
    phi = 2.0 * mc_v.pi * random.random()
    th = math.acos(2.0 * random.random() - 1.0)
    u = math.sin(th) * math.cos(phi)
    v = math.sin(th) * math.sin(phi)
    w = math.cos(th)

#Output
print("Heat Absorbed : {}".format(heata / mc_v.photons))
print("Specular Reflection : {}".format(mc_v.rs))
print("Backscattered Reflection  : {}".format(mc_v.rb / mc_v.photons))
print("Total Reflection : {}".format((mc_v.rs + (mc_v.rb)/ mc_v.photons)))