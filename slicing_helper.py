from pymatgen import Structure, Composition, Lattice
from pymatgen.core.sites import Site, PeriodicSite
from pymatgen.io.vasp.outputs import Elfcar, VolumetricData, Chgcar
from matplotlib import pyplot as plt
import plotly.graph_objects as go


class UnitCellData:
    
    def __init__(self, lattice=None, items=[]):
        self.lattice = lattice
        self.items = items
        
    def cartesian(self):
        #assert self.lattice is not None, 'Must have a lattice'
        new_points = [item[:] for item in self.items]
        for i in new_points:
            i[0:3] = self.lattice.get_cartesian_coords(i[0:3])
        return UnitCellData(items=new_points)
        
    @property
    def x(self):
        return [p[0] for p in self.items]
    
    @property
    def y(self):
        return [p[1] for p in self.items]
    
    @property
    def z(self):
        return [p[2] for p in self.items]
    
    @property
    def vals(self):
        return [p[3] for p in self.items]
    
    @property
    def elf(self):
        return [1 / (1 + p[3] ** 2) for p in self.items]
    
    def slice_x(self, x0):
        minx = min(self.x,key=lambda p: abs(p - x0))
        return UnitCellData(lattice=self.lattice, items=[i for i in self.items if i[0] == minx])
    
    def slice_y(self, y0):
        miny = min(self.y,key=lambda p: abs(p - y0))
        return UnitCellData(lattice=self.lattice, items=[i for i in self.items if i[1] == miny])
    
    def slice_z(self, z0):
        minz = min(self.z,key=lambda p: abs(p - z0))
        return UnitCellData(lattice=self.lattice, items=[i for i in self.items if i[2] == minz])
    
    def plane(self, rel):
        return UnitCellData(lattice=self.lattice, items=[i for i in self.items if rel(i[0], i[1], i[2]) == 0])
    
    def find_point(self, axis, r0, x1, y1):
        if axis == 'x':
            minx = min(self.x,key=lambda p: abs(p - r0))
            miny = min(self.y,key=lambda p: abs(p - x1))
            minz = min(self.z,key=lambda p: abs(p - y1))
            for pt in self.items:
                if pt[0] == minx and pt[1] == miny and pt[2] == minz:
                    return pt[3]
                
def generate_3d(formula, data_dict):
    vd = data_dict[formula]
    pts = list()
    x = 0
    for xaxis in vd.as_dict()['data']['total']:
        y = 0
        for yaxis in xaxis:
            z = 0
            for zaxis in yaxis:
                pts.append([x / (vd.dim[0] - 1), y / (vd.dim[1] - 1), z / (vd.dim[2] - 1), zaxis])
                z += 1
            y += 1
        x += 1
    return UnitCellData(lattice=vd.structure.lattice, items=pts)

def take_slice(d3, which, r0):
    if which == 'x':
        return d3.slice_x(r0)
    elif which == 'y':
        return d3.slice_y(r0)
    elif which == 'z':
        return d3.slice_z(r0)
    else:
        assert 0 == 1
        
def plot_slice(ucb, axis='x', r0=1, f=None):
    ucb_100_f = take_slice(ucb, axis, r0)
    ucb_100 = ucb_100_f.cartesian()
    if axis == 'x':
        X = ucb_100.y
        Y = ucb_100.z
    elif axis == 'y':
        X = ucb_100.x
        Y = ucb_100.z
    else:
        X = ucb_100.x
        Y = ucb_100.y
    Z = ucb_100.vals

    fig = go.Figure(data=go.Scatter(x=X,
                                    y=Y,
                                    mode='markers',
                                    marker_color=Z,
                                    text=Z)) # hover text goes here
    
    #print(len(X), len(Y))

    fig.update_layout(yaxis=dict(scaleanchor="x", scaleratio=1), 
                      title_text=f,
                      width=700, 
                      height=700, plot_bgcolor='rgb(0, 0, 0)',
                      xaxis_showgrid=False, 
                      yaxis_showgrid=False,
                      xaxis_zeroline=False, 
                      yaxis_zeroline=False)
    fig.show()
    
plt.rcParams['figure.figsize'] = [12, 8]
def plot_3d(ucb1, axis='x', r0=1, f=None, savefig=False):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ucb = ucb1.cartesian()
    x = ucb.x
    y = ucb.y
    z = ucb.z
    c = ucb.vals

    img = ax.scatter(x, y, z, c=c, cmap="Accent")
    #fig.colorbar(img)
    plt.title(f, fontsize=18)
    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.tick_params(axis='y', which='both', bottom=False, top=False, labelbottom=False)
    plt.tick_params(axis='z', which='both', bottom=False, top=False, labelbottom=False)
    if savefig:
        pass
        #plt.savefig('images/bader_basins/{}.png'.format(f), dpi=300)
    if not savefig:
        plt.show()
    
def contour_map(ucb, axis='x', r0=1):
    ucb_100_f = take_slice(ucb, axis, r0)
    ucb_100 = ucb_100_f.cartesian()
    Z = np.reshape(ucb_100.vals, (len(set(ucb_100_f.y)), len(set(ucb_100_f.z))))
    fig = go.Figure(data =
    go.Contour(
        z=np.transpose(Z)))
    fig.update_layout(yaxis=dict(scaleanchor="x", scaleratio=1), 
                      plot_bgcolor='rgb(255,255,255)',
                      width=700,
                      height=700)
    fig.show()