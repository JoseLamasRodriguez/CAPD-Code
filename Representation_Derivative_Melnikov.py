import matplotlib.pyplot as plt
import math
import shapely.geometry as sg
import shapely.ops as so
import numpy as np

def strtolst(line):
    return line.strip('\n][').split(', ')

file_data = open("data_derivative_Melnikov.txt","r+")
file_angles = open("angles_derivative_Melnikov.txt","r+")

angle_left = []
angle_right = []
data_left = []
data_right=[]

#Translation into types

for line in file_angles:
  angle_left.append(float(strtolst(line)[0].split(',')[0]))
  angle_right.append(float(strtolst(line)[0].split(',')[1]))

for line in file_data:
  data_left.append(float(strtolst(line)[0].split(',')[0]))
  data_right.append(float(strtolst(line)[0].split(',')[1]))


file_data.close()
file_angles.close()


#Transform into numpy arrays
angle_left = np.array(angle_left)
angle_right = np.array(angle_right)
data_left = np.array(data_left)
data_right = np.array(data_right)

rect_left = []
rect_right=[]

#Plot 
fig, axs = plt.subplots()
plt.grid(which='major', axis='y', zorder=-1.0)
plt.grid(which='major', axis='x', zorder=-1.0)

#Discontinuity at position 736-737 due to the singularity
for i in range(736):
  # plt.plot([angle_left[i],angle_right[i]],[data_right[i],data_right[i]],'b')
  # plt.plot([angle_left[i],angle_right[i]],[data_left[i],data_left[i]],'b')
  pts = sg.Polygon([(angle_left[i],data_left[i]),(angle_left[i],data_right[i]),(angle_right[i],data_right[i]),(angle_right[i],data_left[i])])
  rect_left.append(pts)

for i in range(737,angle_left.size):
  # plt.plot([angle_left[i],angle_right[i]],[data_right[i],data_right[i]],'b')
  # plt.plot([angle_left[i],angle_right[i]],[data_left[i],data_left[i]],'b')
  pts = sg.Polygon([(angle_left[i],data_left[i]),(angle_left[i],data_right[i]),(angle_right[i],data_right[i]),(angle_right[i],data_left[i])])
  rect_right.append(pts)
# pts = sg.Polygon([(angle_left[0],data_left[0]),(angle_left[0],data_right[0]),(angle_right[0],data_right[0]),(angle_right[0],data_left[0])])
# rect.append(pts)

new_shape_left = so.unary_union(rect_left)
new_shape_right = so.unary_union(rect_right)
#exterior coordinates split into two arrays, xs and ys
# which is how matplotlib will need for plotting
xs, ys = new_shape_left.exterior.xy
xxs,yys = new_shape_right.exterior.xy

axs.fill(xs, ys, alpha=0.5, fc='r', ec='none', label= 'CAP result')
axs.fill(xxs,yys, alpha=0.5,fc='r', ec='None')

plt.xlabel(r'$\theta$')
plt.ylabel(r'$\frac{d}{d\theta}\;M_+(\theta)$')
plt.vlines(math.sqrt(2)/3, -7, 1.1, 'black', '--', label='Singularity')
#plt.vlines(2*np.pi - 4./3, 0.370392,0.891352, 'b', '-', label = r'$\theta = -\frac 4 3$')

plt.hlines(0,0, 2*np.pi, 'm','--')
plt.xlim(0,2*np.pi)
plt.ylim(-7,1.1)
plt.plot()
plt.legend(loc='lower right')
plt.show()
# plt.ylim(-17,17)
# plt.show()

