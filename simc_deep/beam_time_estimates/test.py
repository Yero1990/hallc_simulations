import numpy as np
import matplotlib.pyplot as plt

hms_hsize = 4.575  #cm
hms_vsize = 11.646 #cm

#define HMS octagon coordinates (9 points)
hcoord = [
    [ hms_hsize,    hms_vsize/2.],
    [ hms_hsize/2., hms_vsize   ],
    [-hms_hsize/2., hms_vsize   ],
    [-hms_hsize,     hms_vsize/2.],
    [ -hms_hsize,    -hms_vsize/2.],
    [-hms_hsize/2., -hms_vsize],
    [hms_hsize/2., -hms_vsize],
    [hms_hsize,    -hms_vsize/2.],
    [hms_hsize,     hms_vsize/2.]
]

shms_hsize = 17. #cm
shms_vsize = 25. #cm
#define HMS octagon coordinates (9 points)
scoord = [
    [ shms_hsize,     shms_vsize/2.],
    [ shms_hsize/2.,  shms_vsize   ],
    [-shms_hsize/2.,  shms_vsize   ],
    [-shms_hsize,     shms_vsize/2.],
    [ -shms_hsize,   -shms_vsize/2.],
    [-shms_hsize/2., -shms_vsize   ],
    [ shms_hsize/2., -shms_vsize   ],
    [ shms_hsize,    -shms_vsize/2.],
    [ shms_hsize,     shms_vsize/2.]
]


#coord.append(coord[0]) #repeat the first point to create a 'closed loop'

hxs, hys = zip(*hcoord) #create lists of x and y values
exs, eys = zip(*hcoord) #create lists of x and y values

plt.figure()

plt.plot(hxs,hys, color='r')

plt.show() # if you need...
