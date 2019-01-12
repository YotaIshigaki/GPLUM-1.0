from analysis import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def makeSnap(datfile, pngfile, \
             bHeader = True, \
             m_sun = 1.0, \
             xrange = [0.88, 1.12], \
             yrange = [-0.002, 0.047], \
             axis = "a-e", \
             m_cut  = 1.e21/M, \
             m_size = 1.e23/M, \
             size = 8., \
             color = "blue") :
    """
    Make a snapshot figure from a snapshot file.
    
    Arguments
    ----------
    datfile : character strings. snapshot filename
    pngfile : character string. snapshot figure filename
    bHeader : boolian. flag as to whether header exists in snapshot (default True)
    m_sun :   float. mass of a central star (default 1.0)
    xrange :  [float,float]. x direction range of a snapshot figure (default [0.89, 1.11])
    yrange :  [float,float]. y direction range of a snapshot figure (default [-0.002, 0.025])
    axis :    character string. elements of x,y in a snapshot figure (default "a-e")
              "a-e"  -- semimajor axis - eccentricity
              "a-i"  -- semimajor axis - inclination
              "a-ei" -- semimajor axis - root mean squared of eccentricity 
                        and inclination
    m_cut :   float. minimum particle mass to draw in a snapshot figure (default 1.e23 g)
    m_size :  float. mass as a reference for determining sizes of markers (default 1.e20 kg)
    size :    float. marker size of particles with mass m_size (default 10.)
    color :   character string. marker color (default "blue")

    Returns
    ---------
    m_max :   float. largest particle mass
    m_ave :   float. mean particle mass
    ecc_rms : float. root mean square of eccentricities
    inc_rms : float. root mean square of inclinations
    """

    header, pp = readSnap(datfile, bHeader)
    
    time_str = "%0.f yr"%(header.time/(2*PI))

    m_max, m_ave, ecc_rms, inc_rms = clacOrbitalElements(pp, m_sun)

    siz = []
    ax = []
    ecc = []
    inc = []
    for ptcl in pp.values() :
        if ptcl.mass >= m_cut :
            siz.append(((ptcl.mass/m_size)**(1./2.))*size)
            ax.append(ptcl.ax)
            ecc.append(ptcl.ecc)
            inc.append(ptcl.inc)

    # Make Figure
    if axis in ["a-e", "a-i", "a-ei"] and pngfile is not None: 
        plt.rcParams['font.family'] = "Times New Roman"
        plt.rcParams['font.size'] = 17
        
        plt.figure()
        if axis == "a-e" :
            plt.scatter(ax, ecc, s=siz,\
                        c=color, alpha=0.15, linewidths=0.5)
            plt.xlabel("Semi-major Axis (AU)", fontsize=25)
            plt.ylabel("Eccentricity", fontsize=25)
            
        elif axis == "a-i" :
            plt.scatter(ax, inc, s=siz,\
                        c=color, alpha=0.17, linewidths=0.5)
            plt.xlabel("Semi-major Axis (AU)", fontsize=25)
            plt.ylabel("Inclination", fontsize=25)
            
        elif axis == "a-ei" :
            plt.scatter(ax, [sqrt(ecci**2 + inci**2) for (ecci, inci) in zip(ecc, inc)],\
                        s=siz,\
                        c=color, alpha=0.17, linewidths=0.5)
            plt.xlabel("$a$ (AU)", fontsize=25)
            plt.ylabel("$\sqrt{e^2+i^2}$", fontsize=25)

        plt.scatter([],[],alpha=0,label=' ')
        plt.legend(fontsize=15, loc='upper right', frameon=False, title=time_str)
        plt.subplots_adjust(bottom=0.13, left=0.15)

        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.xlim(xrange)
        plt.ylim(yrange)
    
        plt.savefig(pngfile, dpi=130)
        plt.close()

    return [header, m_max, m_ave, ecc_rms, inc_rms]


def makeCumulativeNumber(datfiles, pngfile, \
                         bHeader = True, \
                         xrange = [3.e21, 2.e26], \
                         yrange = [0.5, 1.e5], \
                         labels=None, \
                         colors=["mediumorchid", "blueviolet", "blue", "midnightblue"], \
                         styles=[":","-.","--","-"]) :
    """
    Make a figure of cumulative number distributuion
    
    Arguments
    ----------
    datfiles: list of character strings. snapshot filenames
    pngfile:  character string. snapshot figure filename
    bHeader : boolian. flag as to whether header exists in snapshot (default True)
    xrange:   [float,float]. x direction range of a figure (default [3.e21, 1.26])
    yrange:   [float,float]. y direction range of a figure (default [1.e5, 0.5])
    labels:   list of character strings. labels (default None)
    colors:   list of character strings.  marker colors 
              (default ["mediumorchid", "blueviolet", "blue", "midnightblue"])
    styles:   list of character strings. line styles (default ["-.","--","-"])

    Returns
    ---------
    """

    plt.rcParams['font.family'] = "Times New Roman"
    plt.rcParams['font.size'] = 17

    plt.figure()
    plt.xscale("log")
    plt.yscale("log")

    i = 0
    
    for datfile in datfiles:
        header, pp = readSnap(datfile, bHeader)

        label = " "
        if labels == None :
            label = "%0.f yr"%(header.time/(2*PI))
        else :
            label = labels[i]

        cumulativeNumber = calcCumulativeNumber(pp)

        mass = []
        Nc = []
        for m, nc in cumulativeNumber :
            mass.append(m*M)
            Nc.append(nc-1)
            mass.append(m*M)
            Nc.append(nc)
        mass.append(0)
        Nc.append(cumulativeNumber[-1][1])
        

        plt.plot(mass, Nc, linewidth=1, label=label, \
                 color=colors[i], ls=styles[i])

        i = i + 1

    plt.xlabel("Mass (g)", fontsize=25)
    plt.ylabel("Cumulative Number", fontsize=25)
    plt.legend(fontsize=17, loc='upper right', frameon=False)
    plt.subplots_adjust(bottom=0.13, left=0.15)

    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlim(xrange)
    plt.ylim(yrange)

    plt.savefig(pngfile, dpi=100)
    plt.close()

    return
    

def makeNumberOfCluster(datfiles, pngfile1, pngfile2, \
                        bHeader = True, \
                        xrange1 = [0.7, 5.e2], \
                        yrange1 = [0.7, 1.e5], \
                        xrange2 = [0.7, 5.e1], \
                        yrange2 = [0.7, 1.e5], \
                        labels=None, \
                        colors=["mediumorchid", "blueviolet", "blue", "midnightblue"], \
                        styles=[":","-.","--","-"],\
                        markers=["^","s","o","D","h","*"]) :
    """
    Make a figure of number of cluster size distribution
    
    Arguments
    ----------
    datfiles : list of character strings. snapshot filenames
    pngfile1 :  character string. snapshot figure filename of number of cluster 
    pngfile2 :  character string. snapshot figure filename of number of neighbor 
    bHeader :  boolian. flag as to whether header exists in snapshot (default True)
    xrange1 :  [float,float]. x direction range of a figure of number of cluster 
               (default [0.5, 1.e5])
    yrange1 :  [float,float]. y direction range of a figure of number of cluster 
               (default [0.5, 1.e3])
    xrange2 :  [float,float]. x direction range of a figure of number of neighbor 
               (default [0.5, 1.e2])
    yrange2 :  [float,float]. y direction range of a figure of number of neighbor 
               (default [0.5, 1.e5])
    labels :   list of character strings. labels (default None)
    colors :   list of character strings.  marker colors 
               (default ["mediumorchid", "blueviolet", "blue", "midnightblue"])
    styles:    list of character strings. line styles (default ["-.","--","-"])
    markers :  list of character strings. mark styles 
               (default ["^","s","o","D","h","*"])

    Returns
    ---------
    """

    i = 0

    X1 = []
    Y1 = []
    X2 = []
    Y2 = []
    time = []
    
    for datfile in datfiles:
        header, pp = readSnap(datfile, bHeader)
        time.append(header.time)

        numberOfCluster, numberOfNeighbor = calcNumberOfCluster(pp)

        N_c = []
        N_p = []
        for n_c, n_p in numberOfCluster :
            N_c.append(n_c)
            N_p.append(n_p)

        X1.append(N_c)
        Y1.append(N_p)

        N_c = []
        N_p = []
        for n_c, n_p in numberOfNeighbor :
            N_c.append(n_c)
            N_p.append(n_p)

        X2.append(N_c)
        Y2.append(N_p)


    plt.rcParams['font.family'] = "Times New Roman"
    plt.rcParams['font.size'] = 17

    plt.figure()
    plt.xscale("log")
    plt.yscale("log")

    for i in range(len(datfiles)) :
        label = " "
        if labels == None :
            label = "%0.f yr"%(time[i]/(2*PI))
        else :
            label = labels[i]
            
        plt.plot(X1[i], Y1[i], \
                 linewidth=1, label=label, color=colors[i], \
                 marker=markers[i], ls=styles[i], markersize=8, fillstyle="none")

    plt.xlabel("Size of Neighbor Clusters", fontsize=25)
    plt.ylabel("The Number of Clusters", fontsize=25)
    plt.legend(fontsize=17, loc='upper right', frameon=False)
    plt.subplots_adjust(bottom=0.14, left=0.14)

    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlim(xrange1)
    plt.ylim(yrange1)

    plt.savefig(pngfile1, dpi=100)
    plt.close()

    
    plt.rcParams['font.family'] = "Times New Roman"
    plt.rcParams['font.size'] = 17

    plt.figure()
    plt.xscale("log")
    plt.yscale("log")

    for i in range(len(datfiles)) :
        label = " "
        if labels == None :
            label = "%0.f yr"%(time[i]/(2*PI))
        else :
            label = labels[i]
            
        plt.plot(X2[i], Y2[i], \
                 linewidth=1, label=label, color=colors[i], \
                 marker=markers[i], ls=styles[i], markersize=8, fillstyle="none")

    plt.xlabel("The Number of Neighbors", fontsize=25)
    plt.ylabel("The Number of Particles", fontsize=25)
    plt.legend(fontsize=17, loc='upper right', frameon=False)
    plt.subplots_adjust(bottom=0.14, left=0.14)

    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlim(xrange2)
    plt.ylim(yrange2)

    plt.savefig(pngfile2, dpi=100)
    plt.close()
        
    return


def makeCumulativeNumberOfCluster(datfiles, pngfile1, pngfile2, \
                                  bHeader = True, \
                                  xrange1 = [0.7, 5.e2], \
                                  yrange1 = [0.7, 1.e5], \
                                  xrange2 = [0.7, 5.e1], \
                                  yrange2 = [0.7, 1.e5], \
                                  labels=None, \
                                  colors=["mediumorchid", "blueviolet", "blue", "midnightblue"], \
                                  styles=[":","-.","--","-"], \
                                  markers=["^","s","o","D","h","*"]) :
    """
    Make a figure of cumulative number of cluster size distribution
    
    Arguments
    ----------
    datfiles : list of character strings. snapshot filenames
    pngfile1 :  character string. snapshot figure filename of number of cluster 
    pngfile2 :  character string. snapshot figure filename of number of neighbor 
    bHeader :  boolian. flag as to whether header exists in snapshot (default True)
    xrange1 :  [float,float]. x direction range of a figure of number of cluster 
               (default [0.5, 1.e5])
    yrange1 :  [float,float]. y direction range of a figure of number of cluster 
               (default [0.5, 1.e3])
    xrange2 :  [float,float]. x direction range of a figure of number of neighbor 
               (default [0.5, 1.e2])
    yrange2 :  [float,float]. y direction range of a figure of number of neighbor 
               (default [0.5, 1.e5])
    labels :   list of character strings. labels (default None)
    colors :   list of character strings.  marker colors 
               (default ["mediumorchid", "blueviolet", "blue", "midnightblue"])
    styles:    list of character strings. line styles (default ["-.","--","-"])
    markers :  list of character strings. mark styles 
               (default ["^","s","o","D","h","*"])

    Returns
    ---------
    """

    i = 0

    X1 = []
    Y1 = []
    X2 = []
    Y2 = []
    time = []
    
    for datfile in datfiles:
        header, pp = readSnap(datfile, bHeader)
        time.append(header.time)

        numberOfCluster, numberOfNeighbor = calcNumberOfCluster(pp)
        numberOfCluster.reverse()
        numberOfNeighbor.reverse()

        N_c = []
        N_p = []
        n = 0
        for n_c, n_p in numberOfCluster :
            #N_c.append(n_c)
            #N_p.append(n)
            #n = n + n_p * n_c
            n = n + n_p
            if n_c > 0 :
                N_c.append(n_c)
                N_p.append(n)
        #N_c.append(0)
        #N_p.append(n)

        X1.append(N_c)
        Y1.append(N_p)

        N_c = []
        N_p = []
        n = 0
        for n_c, n_p in numberOfNeighbor :
            #N_c.append(n_c)
            #N_p.append(n)
            n = n + n_p
            if n_c > 0 :
                N_c.append(n_c)
                N_p.append(n)
        #N_c.append(0)
        #N_p.append(n)

        X2.append(N_c)
        Y2.append(N_p)


    plt.rcParams['font.family'] = "Times New Roman"
    plt.rcParams['font.size'] = 17

    plt.figure()
    plt.xscale("log")
    plt.yscale("log")

    for i in range(len(datfiles)) :
        label = " "
        if labels == None :
            label = "%0.f yr"%(time[i]/(2*PI))
        else :
            label = labels[i]
            
        plt.plot(X1[i], Y1[i], \
                 linewidth=1, label=label, color=colors[i], ls=styles[i],\
                 marker=markers[i], markersize=6, fillstyle="none")

    plt.xlabel("Size of Neighbor Cluster", fontsize=25)
    plt.ylabel("Cumulative Number of Clusters", fontsize=25)
    plt.legend(fontsize=17, loc='upper right', frameon=False)
    plt.subplots_adjust(bottom=0.14, left=0.14)

    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlim(xrange1)
    plt.ylim(yrange1)

    plt.savefig(pngfile1, dpi=100)
    plt.close()

    
    plt.rcParams['font.family'] = "Times New Roman"
    plt.rcParams['font.size'] = 17

    plt.figure()
    plt.xscale("log")
    plt.yscale("log")

    for i in range(len(datfiles)) :
        label = " "
        if labels == None :
            label = "%0.f yr"%(time[i]/(2*PI))
        else :
            label = labels[i]
            
        plt.plot(X2[i], Y2[i], \
                 linewidth=1, label=label, color=colors[i], ls=styles[i],\
                 marker=markers[i], markersize=6, fillstyle="none")

    plt.xlabel("The Number of Neighbors", fontsize=25)
    plt.ylabel("Cumulative Number", fontsize=25)
    plt.legend(fontsize=17, loc='upper right', frameon=False)
    plt.subplots_adjust(bottom=0.14, left=0.14)

    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlim(xrange2)
    plt.ylim(yrange2)

    plt.savefig(pngfile2, dpi=100)
    plt.close()
        
    return


def makeEccentricityDistribution(datfiles, pngfile, \
                                 bHeader = True, \
                                 xrange = [3.e21, 2.e26], \
                                 yrange = [0.0005, 0.07], \
                                 labels=None, \
                                 colors=["mediumorchid", "blueviolet", "blue", "midnightblue"], \
                                 styles=[":","-.","--","-"],\
                                 markers=["^","s","o","D","h","*"]) :
    """
    Make a figure of root mean square of eccentricities and inclinations against mass
    Use mass bin whose width m_bin(m/m_bin)^2.5
    
    Arguments
    ----------
    datfiles: list of character strings. snapshot filenames
    pngfile:  character string. snapshot figure filename
    bHeader : boolian. flag as to whether header exists in snapshot (default True)
    xrange:   [float,float]. x direction range of a figure (default [3.e21, 1.26])
    yrange:   [float,float]. y direction range of a figure (default [1.e5, 0.5])
    labels:   list of character strings. labels (default None)
    colors:   list of character strings.  marker colors 
              (default ["mediumorchid", "blueviolet", "blue", "midnightblue"])
    styles:   list of character strings. line styles 
              (default ["-.","--","-"])
    marks :   list of character strings. mark styles 
              (default ["^","s","o","D","h","*"])

    Returns
    ---------
    """

    plt.rcParams['text.usetex'] = True
    plt.rcParams['text.latex.unicode'] = True
    plt.rcParams['font.family'] = "Times"
    plt.rcParams['font.size'] = 17

    plt.figure()
    plt.xscale("log")
    plt.yscale("log")
    plt.yticks([0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, \
                0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1], \
               [0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, \
                0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1])

    i = 0
    
    for datfile in datfiles:
        header, pp = readSnap(datfile, bHeader)

        label = " "
        if labels == None :
            label = "%0.f yr"%(header.time/(2*PI))
        else :
            label = labels[i]

        eccinc = calcEccentricityDistribution(pp)

        mass = []
        ecc = []
        for m, e, dmmy1, dmmy2 in eccinc :
            if e != None :
                mass.append(m*M)
                ecc.append(e)

        plt.plot(mass, ecc, linewidth=1, label=label, \
                 color=colors[i], marker=markers[i], ls=styles[i], fillstyle="none")

        i = i + 1

    plt.xlabel("Mass (g)", fontsize=25)
    plt.ylabel(r"$\langle e^2\rangle^{1/2}$", fontsize=25)
    plt.legend(fontsize=17, loc='lower left', frameon=False)
    plt.subplots_adjust(bottom=0.13, left=0.18)

    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlim(xrange)
    plt.ylim(yrange)

    plt.savefig(pngfile, dpi=100)

    plt.rcParams['text.usetex'] = False
    plt.rcParams['text.latex.unicode'] = False
    plt.close()

    return

def makeEccentricityDistribution2(datfile, pngfile, \
                                  bHeader = True, \
                                  xrange = [3.e21, 2.e26], \
                                  yrange = [0.0002, 0.07], \
                                  labels=None, \
                                  colors=["blue", "red"], \
                                  styles=["-","--"],\
                                  markers=["o","s"]) :
    """
    Make a figure of root mean square of eccentricities and inclinations against mass
    
    Arguments
    ----------
    datfile:  character strings. snapshot filename
    pngfile:  character string. snapshot figure filename
    bHeader : boolian. flag as to whether header exists in snapshot (default True)
    xrange:   [float,float]. x direction range of a figure (default [3.e21, 1.26])
    yrange:   [float,float]. y direction range of a figure (default [1.e5, 0.5])
    labels:   list of character strings. labels (default None)
    colors:   list of character strings.  marker colors 
              (default ["blue", "red"])
    styles:   list of character strings. line styles 
              (default ["-","--"])
    marks :   list of character strings. mark styles 
              (default ["o","s"])

    Returns
    ---------
    """

    plt.rcParams['text.usetex'] = True
    plt.rcParams['text.latex.unicode'] = True
    plt.rcParams['font.family'] = "Times"
    plt.rcParams['font.size'] = 17

    plt.figure()
    plt.xscale("log")
    plt.yscale("log")
    plt.yticks([0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, \
                0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1], \
               [0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, \
                0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1])
    
    header, pp = readSnap(datfile, bHeader)

    label1 = " "
    label2 = " "
    if labels == None :
        label1 = r"$\langle e^2\rangle^{1/2}$"
        label2 = r"$\langle i^2\rangle^{1/2}$"
    else :
        label1 = labels[0]
        label2 = labels[1]
        
    eccinc = calcEccentricityDistribution(pp)

    mass = []
    ecc = []
    inc = []
    for m, e, i, dmmy2 in eccinc :
        if e != None :
            mass.append(m*M)
            ecc.append(e)
            inc.append(i)
            
    plt.plot(mass, ecc, linewidth=1, label=label1, \
             color=colors[0], marker=markers[0], ls=styles[0], fillstyle="none")
    plt.plot(mass, inc, linewidth=1, label=label2, \
             color=colors[1], marker=markers[1], ls=styles[1], fillstyle="none")


    plt.xlabel("Mass (g)", fontsize=25)
    plt.ylabel(r"$\langle e^2\rangle^{1/2}, \langle i^2\rangle^{1/2}$", fontsize=25)
    plt.legend(fontsize=17, loc='lower left', frameon=False)
    plt.subplots_adjust(bottom=0.13, left=0.18)

    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlim(xrange)
    plt.ylim(yrange)

    plt.savefig(pngfile, dpi=100)

    plt.rcParams['text.usetex'] = False
    plt.rcParams['text.latex.unicode'] = False
    plt.close()

    return
