
from zone import Zone
import numpy as np
import numpy.random as rnd

class ConeGene:
    def __init__(self, lcone_gene, d2_gene, d1, n_points=5, alpha_max=45):    
        self.normalizable = False
        self.name = 'Cone'
        self.lcone_gene = lcone_gene
        self.d2_gene = d2_gene
        self.d1 = d1
        self.n_points = n_points
        self.alpha_max = alpha_max
    
    def get_rnd_value(self):
        Lcone = self.lcone_gene.get_rnd_value()
        d2 = self.d2_gene.get_rnd_value()
        points_xd = self.get_points(0, Lcone, self.d1, d2)
        return {
            'Lcone': Lcone,
            'd2': d2,
            'points_xd': points_xd
        }
        
    def get_points(self, x1, x2, d1, d2):
        z = Zone(x1, x2, d1/2, d2/2, alpha_max=self.alpha_max)
        zs = [z]
        res = []
        for i in range(self.n_points):
            areas = np.array([z.area for z in zs])
            areas /= np.sum(areas)
            i = rnd.choice(len(zs), p=areas)
            z = zs.pop(i)
            x, r = z.get_rnd_x_r()
            res.append((x, r*2))
            z1, z2 = z.get_split_zones(x, r)
            zs.append(z1)
            zs.append(z2)
        return np.array(sorted(res))
    
    def cross(self, gene1, gene2):
        t = np.random.random()
        Lcone = (1 - t) * gene1['Lcone'] + t * gene2['Lcone']
        d2 =    (1 - t) * gene1['d2']    + t * gene2['d2']
        points_xd = (1 - t) * gene1['points_xd']    + t * gene2['points_xd']
        return {
            'Lcone': Lcone,
            'd2': d2,
            'points_xd': points_xd
        }
        
    def mutate(self, mu):
        what_mutate = np.random.choice(['Lcone', 'd2', 'points_xd'], p=[0.25,0.1,0.65])
        if what_mutate == 'Lcone':
            Lcone = self.lcone_gene.mutate(mu['Lcone'])
            d2 = mu['d2']
            t = Lcone / mu['Lcone']
            xs = mu['points_xd'][:,0] * t
            ds = mu['points_xd'][:,1]
            points_xd = np.stack((xs, ds)).T
        elif what_mutate == 'd2':
            Lcone = mu['Lcone']
            d2 = self.d2_gene.mutate(mu['d2'])
            d1 = self.d1
            t = (d1 - d2) / (d1 - mu['d2'])
            xs = mu['points_xd'][:,0] 
            ds = d1 - (d1 - mu['points_xd'][:,1]) * t
            points_xd = np.stack((xs, ds)).T
        else:
            Lcone = mu['Lcone']
            d2 = mu['d2']
            d1 = self.d1
            i_mutate = np.random.choice(self.n_points)
            if i_mutate == 0:
                x1 = 0
                r1 = d1/2
            else:
                x1, r1  = mu['points_xd'][i_mutate-1]
                r1 /= 2
            
            if i_mutate == self.n_points - 1:
                x2 = Lcone
                r2 = d2/2
            else:
                x2, r2  = mu['points_xd'][i_mutate+1]
                r2 /= 2
            
            z = Zone(x1, x2, r1, r2, self.alpha_max)
            x2, r2 = z.get_rnd_x_r()
            points_xd = np.array(mu['points_xd'])
            points_xd[i_mutate] = (x2, r2*2)
        return {
            'Lcone': Lcone,
            'd2': d2,
            'points_xd': Zone(0, Lcone, self.d1/2, d2/2, self.alpha_max).fix_points(points_xd)
        }
    
    def plot(self, ax, gene, **kwargs):
        xs = [0, *(gene['points_xd'][:,0]), gene['Lcone']]
        ys = [self.d1/2, *(gene['points_xd'][:,1]/2), gene['d2']/2]
        ax.plot(xs, ys, **kwargs)
            
            
                
