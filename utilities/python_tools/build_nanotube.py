#
# source activate param
#
import numpy as np
from ase.build import nanotube
import diameter

m = int(input("m: "))
n = int(input("n: "))
length = int(input("length: "))


diameter.cnt_radius(m,n)
cnt = nanotube(m,n,length=length)
cnt.write("nanotube-" + str(m) + "-" + str(n) + "-" + str(length) + ".xyz")

print("Lenght {:.8f}".format(np.max(cnt.positions[:,2])))
