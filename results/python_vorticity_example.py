import numpy
import matplotlib.pyplot as plt
import cmocean

data = numpy.load('vortex-fields-2023-05-03-B1.npz')

# set up grid
X,Y = numpy.meshgrid(data['grid_x'], data['grid_y'])
Ntimesteps = data['Omega'].shape[2]

## SET UP FIGURE
plt.ion() # interactive plots
fig = plt.figure()
ax1 = fig.add_subplot(111)


## SET UP ANIMATION LOOP
def animate(idx):
    ax1.clear()
    cax1 = ax1.pcolormesh(X,Y,data['Omega'][:,:,idx], 
                          cmap=cmocean.cm.balance,
                          vmin=-150, vmax=150)
    return cax1

## DRAW THE FIRST FRAME AND SET UP COLORBAR
fig.colorbar(animate(0), ax=ax1)

## ANIMATE THE PLOT
from matplotlib.animation import FuncAnimation
ani = FuncAnimation(fig, animate, 
                    frames=Ntimesteps, 
                    interval=1000/10, # ten frames per second
                    repeat=True)
plt.show()


