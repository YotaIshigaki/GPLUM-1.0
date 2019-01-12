from classdef import *

PI = m.acos(-1.)
L = 14959787070000
M = 1.9884e33

def readSnap(filename, bHeader=True):
    """
    Read data of particles from a snapshot file

    Arguments
    ----------
    filename : character string. snapshot filename
    bHeader :  boolian. flag as to whether header exists in snapshot (default True)
    
    Returns
    ----------
    header : Header. header in snapshot
    pp :     dictionary {int : Particle}. particle system. key indicats particle id
    """
    
    i = 0
    pp = {}
    header = (0.,0,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.)

    with open(filename) as f:
        for line in f:
            part = [p.strip() for p in line.split("\t")]

            if i==0 and bHeader :
                header = Header(float(part[0]),\
                                int(part[1]),\
                                float(part[2]),\
                                float(part[3]),\
                                float(part[4]),\
                                float(part[5]),\
                                float(part[6]),\
                                float(part[7]),\
                                float(part[8]),\
                                float(part[9]),\
                                float(part[10]),\
                                float(part[11]) )
                
            else :
                idx = int(part[0])
                ptcl = Particle(float(part[1]),\
                                float(part[2]),\
                                float(part[3]),\
                                float(part[4]),\
                                float(part[5]),\
                                float(part[6]),\
                                float(part[7]) )
                ptcl.setNeighbor(int(part[8]),\
                                 int(part[9]) )
                pp[idx] = ptcl;

            i += 1

    return [header, pp]

def readSnapForParticles(filename, ids, bHeader=True):
    """
    Read data of paticular particles from a snapshot file

    Arguments
    ----------
    filename : character string. snapshot filename
    ids :      list of int. list of particle id of particles to be read
    bHeader :  boolian. flag as to whether header exists in snapshot
    
    Returns
    ----------
    header : Header. header in snapshot
    pp :     dictionary {int : Particle}. particle system. key indicats particle id
    """

    i = 0
    pp = {}
    
    with open(filename) as f:
        for line in f:
            part = [p.strip() for p in line.split("\t")]

            if i==0 and bHeader:
                header = Header(float(part[0]),\
                                int(part[1]),\
                                float(part[2]),\
                                float(part[3]),\
                                float(part[4]),\
                                float(part[5]),\
                                float(part[6]),\
                                float(part[7]),\
                                float(part[8]),\
                                float(part[9]),\
                                float(part[10]),\
                                float(part[11]) )
            else :
                idx = int(part[0])
                if idx in ids :
                    ptcl = Particle(float(part[1]),\
                                    float(part[2]),\
                                    float(part[3]),\
                                    float(part[4]),\
                                    float(part[5]),\
                                    float(part[6]),\
                                    float(part[7]) )
                    ptcl.setNeighbor(int(part[8]),\
                                     int(part[9]) )
                    pp[idx] = ptcl
                    
            i += 1

    return [header, pp]


def clacOrbitalElements(pp, m_sun=1.) :
    """
    Calculate orbital elements of paticles

    Arguments
    ----------
    pp :     dictionary {int : Particle}. particle system. key indicats particle id
    m_sun :  float. mass of a central star (default 1.0)
    
    Returns
    ---------
    m_max:   float. largest particle mass
    m_ave:   float. mean particle mass
    ecc_rms: float. root mean square of eccentricities
    inc_rms: float. root mean square of inclinations
    """

    n_body = len(pp)

    m_max = 0.
    m_ave = 0.
    ecc_rms = 0.
    inc_rms = 0.
    for ptcl in pp.values():
        ax, ecc, inc = ptcl.getOrbitalElement(m_sun)
        ecc_rms = ecc_rms + ecc*ecc
        inc_rms = inc_rms + inc*inc
        m_ave = m_ave + ptcl.mass
        if ptcl.mass > m_max :
            m_max = ptcl.mass

    ecc_rms = m.sqrt(ecc_rms / n_body)
    inc_rms = m.sqrt(inc_rms / n_body)
    m_ave = m_ave / n_body

    return [m_max, m_ave, ecc_rms, inc_rms]


def calcEnergy(pp, m_sun=1., eps=0.):
    """
    Calculate energy of particle system

    Arguments
    ----------
    pp :     dictionary {int : Particle}. particle system. key indicats particle id
    m_sun :  float. mass of a central star (default 1.0)
    eps      float. softening parameter (default 0.0)
    
    Returns
    ----------
    etot :        float. total energy
    ekin :        float. kinetic energy
    ephi_sun :    float. potential energy by central star
    ephi_planet : float. potential energy by interaction of particles
    """

    ekin = 0.
    ephi_sun = 0.
    ephi_planet = 0.
    
    for id, ptcl in pp.items():
        dr2 = ptcl.pos2() + eps**2
        dv2 = ptcl.vel2()
        
        ekin     += ptcl.mass * dv2
        ephi_sun += ptcl.mass / sqrt(dr2)

        for id2, ptcl2 in pp.items():
            if id < id2 :
                dr2 = ptcl.pos2(ptcl2) + eps**2               
                ephi_planet += ptcl.mass * ptcl2.mass / sqrt(dr2)
                
    ekin *= 0.5
    ephi_sun *= m_sun
    etot = ekin + ephi_sun + ephi_planet

    return [etot, ekin, ephi_sun, ephi_planet]


def calcNumberOfCluster(pp):
    """
    Make dictionary of the number of paritcles in each neighbor cluster

    Arguments
    ----------
    pp :     dictionary {int : Particle}. particle system. key indicats particle id
    
    Returns
    ----------
    numberOfCluster : list of [int, int]. set of the number of particles in a cluster 
                      and the number of clusters.
    numberOfCluster : list of [int, int]. set of the number of particles and 
                      and the number of neighbors of a particle
    """

    n_cluster = []
    n_neighbor = []
    for ptcl in pp.values():
        n_cluster.append(ptcl.cluster)
        n_neighbor.append(ptcl.neighbor)

    n_cluster.sort()
    n_neighbor.sort()

    numberOfCluster = []
    numberOfNeighbor = []

    n_c = 1
    i = 0
    while n_c < n_cluster[-1]+1 :
        n_p = 0
        while i < len(n_cluster) and n_c == n_cluster[i] :
            i = i + 1
            n_p = n_p + 1
            
        numberOfCluster.append([n_c, n_p/n_c])
        n_c = n_c + 1
        
    #numberOfCluster.append([n_c, 0])

    n_c = 0
    i = 0
    while n_c < n_neighbor[-1]+1 :
        n_p = 0
        while i < len(n_neighbor) and n_c == n_neighbor[i] :
            i = i + 1
            n_p = n_p + 1
            
        numberOfNeighbor.append([n_c, n_p])
        n_c = n_c + 1

    #numberOfNeighbor.append([n_c, 0])

    return [numberOfCluster, numberOfNeighbor]


def calcCumulativeNumber(pp):
    """
    Make dictionary of the cumulative number distributuion.

    Arguments
    ----------
    pp :     dictionary {int : Particle}. particle system. key indicats particle id.
    
    Returns
    ----------
    cumulativeNumber : list of [float, int]. set of the mass of particles and 
                       the cumulative number
    """

    mass = []
    for ptcl in pp.values():
        mass.append(ptcl.mass)

    mass.sort()
    mass.reverse()

    cumulativeNumber = []
    n = 1
    for m in mass :
        cumulativeNumber.append([m,n])
        n = n + 1

    return cumulativeNumber

def countNumberOfParticleSomeCondition(pp, condition, param) :
    """
    Count number of particle satisfy some condition

    Arguments
    ----------
    pp :        dictionary {int : Particle}. particle system. key indicats particle id
    condition : function(Particle, list of parameters). function to be return boolian
    param :     list. list of parameter
    
    Returns
    ----------
    numberOfParticle :     int. number of particle satisfy the condition
    totalMassOfParticle :  float. total mass of  particle satisfy the condition
    """

    numberOfParticle = 0
    totalMassOfParticle = 0.
    for ptcl in pp.values():
        if condition(ptcl, param) == True :
            numberOfParticle = numberOfParticle + 1
            totalMassOfParticle = totalMassOfParticle + ptcl.mass

    return [numberOfParticle, totalMassOfParticle]

def calcEccentricityDistribution(pp, \
                                 m_min=4.e21/M, \
                                 m_sun=1.0) :
    """
    Calculate root mean square of eccentricities and inclinations against mass

    Arguments
    ----------
    pp :      dictionary {int : Particle}. particle system. key indicats particle id
    m_min :  float. the minimum mass (default 4.e21 g)
    m_sun :  float. mass of a central star (default 1.0)
    
    Returns
    ----------
    eccinc :  list of [float, float, float, int]. the mass, root mean square of 
              eccentricities and inclination and the number of particle in the mass bin
    """

    Ptcl = []
    for ptcl in pp.values():
        ax, ecc, inc = ptcl.getOrbitalElement(m_sun)
        Ptcl.append([ptcl.mass, ecc, inc])

    Ptcl.sort()
    
    m_max = Ptcl[-1][0]
    m0 = m_min
    m1 = m0*2.
    
    eccinc = []
    i = 0
    while m0 < m_max :
        rms_ecc = 0.
        rms_inc = 0.
        n_p = 0
        mass = 0.
        while i < len(Ptcl) and Ptcl[i][0] < m1 :
            rms_ecc = rms_ecc + Ptcl[i][1]**2
            rms_inc = rms_inc + Ptcl[i][2]**2
            mass = mass + Ptcl[i][0]
            n_p = n_p + 1
            i = i + 1

        if n_p > 0 :
            rms_ecc = m.sqrt(rms_ecc/n_p)
            rms_inc = m.sqrt(rms_inc/n_p)
            mass = mass / n_p
        else :
            rms_ecc = None
            rms_inc = None
            mass = None
        
        eccinc.append([mass, rms_ecc, rms_inc, n_p])

        m0 = m1
        m1 = m0*2.

    return eccinc
        
    
        

    
    

    

    
        
        
    

    

    
