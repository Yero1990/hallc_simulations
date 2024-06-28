import matplotlib.pyplot as plt


# SHMS Octagonal Collimator Size (Each of the octagonal points is a multiple of 1 or 1/2 of these values)

shms_hsize = 8.5
shms_vsize = 12.5

hms_hsize = 4.575
hms_vsize = 11.646

  
coord_shms = [[  shms_hsize,     shms_vsize/2.],
         [  shms_hsize/2.,  shms_vsize   ],
         [ -shms_hsize/2.,  shms_vsize   ],
         [ -shms_hsize,     shms_vsize/2.],
         [ -shms_hsize,    -shms_vsize/2.],
         [ -shms_hsize/2., -shms_vsize   ],
         [  shms_hsize/2., -shms_vsize   ],
         [  shms_hsize,    -shms_vsize/2.],
         [  shms_hsize,     shms_vsize/2.]]



coord_hms = [[  hms_hsize,     hms_vsize/2.],
         [  hms_hsize/2.,  hms_vsize   ],
         [ -hms_hsize/2.,  hms_vsize   ],
         [ -hms_hsize,     hms_vsize/2.],
         [ -hms_hsize,    -hms_vsize/2.],
         [ -hms_hsize/2., -hms_vsize   ],
         [  hms_hsize/2., -hms_vsize   ],
         [  hms_hsize,    -hms_vsize/2.],
         [  hms_hsize,     hms_vsize/2.]]


coord.append(coord[0])
xs, ys = zip(*coord)
plt.figure()
plt.plot(xs,ys)
plt.show()
