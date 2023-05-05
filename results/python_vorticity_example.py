import numpy
import matplotlib.pyplot as plt
import cmocean

dataname = 'vortex-fields-2023-05-05-vU.npz'
data = numpy.load(dataname)

# set up grid
X,Y = numpy.meshgrid(
    data['grid_x'], 
    data['grid_y'],
    indexing='ij' # THIS IS IMPORTANT -- 
                  # MATLAB/JULIA and PYTHON USE DIFFERENT  CONVENTION
                  # MATLAB columns = x, rows = y
                  # NUMPY (without this parameter) columns = y, rows = x
    )
Ntimesteps = data['Omega'].shape[2]

## SET UP FIGURE
plt.ion() # interactive plots
fig = plt.figure()
ax1 = fig.add_subplot(111)

## SET UP ANIMATION LOOP
def animate_vorticity(idx):
    ax1.clear()
    cax1 = ax1.pcolormesh(X,Y,data['Omega'][:,:,idx], 
                          cmap=cmocean.cm.balance,
                          vmin=-150, vmax=150)
    Vx = data['Vx'][:,:,idx]
    Vy = data['Vy'][:,:,idx]
    
    qax1 = ax1.quiver(X,Y,
                      Vx,
                      Vy, scale=150, alpha=0.25)
    ax1.set_title('idx = %d -- t = %f' % (idx, data['t'][idx]) )
    return cax1, qax1


## DRAW THE FIRST FRAME AND SET UP COLORBAR
#fig.colorbar(animate_velocity(0), ax=ax1)

## ANIMATE THE PLOT
from matplotlib.animation import FuncAnimation
ani = FuncAnimation(fig, animate_vorticity, 
                    frames=Ntimesteps, 
                    interval=1000/10, # ten frames per second
                    repeat=True)

try:
    animationname = dataname + 'gif'
    print('Saving to '+ animationname)
    ani.save(animationname,fps=10)
except:
    print('Problems while trying to save animation')


plt.show()

