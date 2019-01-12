from makefigure import *
import sys
import os


def main() :
    dirname = sys.argv[1]
    Range = [int(sys.argv[2]), int(sys.argv[3])]

    makeFigureAll(dirname, Range)

    return


def makeFigureAll(dirname,\
                  Range = None,\
                  axis = "a-e",
                  bSnap = False):
    """
    Make snapshot figures and other figures in directory

    Arguments
    ----------
    dirname :  character string. directory name
    Rnage :    [int, int]. range of snapshot file number (default None)
    axis :     character string. elements of x,y in a snapshot figure (default "a-e")
    bSnap :    boolian. flag as to whether to output Snapshot (default False)


    Returns
    ----------
    """

    figure_dir = dirname + "/figure"
    
    if not os.path.exists(figure_dir):
        os.makedirs(figure_dir)

    snaps = []
    if Range == None :
        files_all = os.listdir(dir)
        files = [f for f in files_all if os.path.isfile(os.path.join(dir, f))]
        snaps = [f for f in files if f[0:4]=="snap"]
        snaps.sort()
    else :
        for i in range(Range[0],Range[1]+1) :
            snaps.append("snap"+format(i, '0>6')+".dat")

    midsnaps = []
    for i in range(4) :
        midsnaps.append(snaps[min(int(round(float(len(snaps))/4*(i+1))), len(snaps)-1)])
    midsnaps.sort()

    largestid = []
    if True :
        header, pp = readSnap(dirname + "/" + midsnaps[-1])
        massid = []
        for id, ptcl in pp.items():
            massid.append([ptcl.mass, id])
        massid.sort()
        massid.reverse()
            
        largestid.append(massid[0][1])
        largestid.append(massid[1][1])
        largestid.append(massid[2][1])
        largestid.append(massid[3][1])
        largestid.append(massid[4][1])
        print(largestid[0],largestid[1],largestid[2],largestid[3],largestid[4])

    tail = ""
    if axis in ["a-e", "a-i", "a-ei", "e-i"] :
        if axis == "a-e" : tail = "a_e"
        elif axis == "a-i" : tail = "a_i"
        elif axis == "a-ei" : tail = "a_ei"
        
    for snap in midsnaps :
        snapfile = dirname + "/" + snap
        figurefile = figure_dir + "/" + snap[0:-4] + tail + ".png"

        makeSnap(snapfile, figurefile)
        
    makeCumulativeNumber([dirname+"/"+snap for snap in midsnaps], \
                          figure_dir +"/"+"CumulativeNumber.png")
    makeCumulativeNumberOfCluster([dirname+"/"+snap for snap in midsnaps], \
                                  figure_dir +"/"+"NumberOfCluster.png", \
                                  figure_dir +"/"+"NumberOfNeighbor.png")
    makeEccentricityDistribution([dirname+"/"+snap for snap in midsnaps], \
                                 figure_dir +"/"+"EccentricityDistribution.png")

    #"""
    datfile  = dirname + "/outcomes.dat"
    datfile1 = dirname + "/migration.dat" 
    fout  = open(datfile, 'w')
    fout1 = open(datfile1, 'w')

    time_v = []
    m_max_v = []
    m_ave_v = []
    ecc_rms_v = []
    inc_rms_v = []
    n_body_v = []
    for snap in snaps :
        
        snapfile = dirname + "/" + snap
        figurefile = figure_dir + "/" + snap[0:-4] + tail + ".png"

        if bSnap == False : figurefile = None
        
        header, m_max, m_ave, ecc_rms, inc_rms = makeSnap(snapfile, figurefile)
        
        time_v.append(header.time)
        m_max_v.append(m_max)
        m_ave_v.append(m_ave)
        ecc_rms_v.append(ecc_rms)
        inc_rms_v.append(inc_rms)
        n_body_v.append(header.n_body)

        fout.write(str(header.time) + "\t" \
                   + str(header.n_body) + "\t" \
                   + str(m_max) + "\t" \
                   + str(m_ave) + "\t" \
                   + str(ecc_rms) + "\t" \
                   + str(inc_rms) + "\n" )

        header, pp = readSnapForParticles(snapfile, largestid)
        for ptcl in pp.values(): ptcl.getOrbitalElement()
        fout1.write(str(header.time) + "\t" \
                    + str(pp[largestid[0]].mass) + "\t" \
                    + str(pp[largestid[0]].ax) + "\t" \
                    + str(pp[largestid[0]].ecc) + "\t" \
                    + str(pp[largestid[0]].inc) + "\t" \
                    + str(pp[largestid[1]].mass) + "\t" \
                    + str(pp[largestid[1]].ax) + "\t" \
                    + str(pp[largestid[1]].ecc) + "\t" \
                    + str(pp[largestid[1]].inc) + "\t" \
                    + str(pp[largestid[2]].mass) + "\t" \
                    + str(pp[largestid[2]].ax) + "\t" \
                    + str(pp[largestid[2]].ecc) + "\t" \
                    + str(pp[largestid[2]].inc) + "\t" \
                    + str(pp[largestid[3]].mass) + "\t" \
                    + str(pp[largestid[3]].ax) + "\t" \
                    + str(pp[largestid[3]].ecc) + "\t" \
                    + str(pp[largestid[3]].inc) + "\t" \
                    + str(pp[largestid[4]].mass) + "\t" \
                    + str(pp[largestid[4]].ax) + "\t" \
                    + str(pp[largestid[4]].ecc) + "\t" \
                    + str(pp[largestid[4]].inc) + "\n" )

        print(snap)

    fout.close()
    fout1.close()

    plt.rcParams['text.usetex'] = True
    plt.rcParams['text.latex.unicode'] = True
    plt.rcParams['font.family'] = "Times"
    
    plt.figure()
    plt.xscale("log")
    plt.yscale("log")
    plt.plot([time/(2.*PI) for time in time_v], [mass*M for mass in m_max_v], \
             linewidth=1, \
             label=r"$M_{\mathrm{max}}$", color='b', ls='-')
    plt.plot([time/(2.*PI) for time in time_v], [mass*M for mass in m_ave_v], \
             linewidth=1, \
             label=r"$\langle m\rangle$", color='b', ls='--')
    plt.ylabel(r"$M_{\mathrm{max}}, \langle m\rangle$ (kg)", fontsize=25)
    plt.xlabel(r"Time (year)", fontsize=25)
    plt.legend(fontsize=17, loc='upper left', frameon=False)
    plt.xlim([min(time_v)/(2.*PI),max(time_v)/(2.*PI)])
    plt.ylim([3.e21,2.e26])
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.subplots_adjust(bottom=0.13, left=0.15)
    
    pngfile = figure_dir + "/massEvolution.png"
    plt.savefig(pngfile, dpi=100)
    plt.close()

    plt.figure()
    plt.xscale("log")
    plt.yscale("log")    
    #plt.yticks([1.e-4, 2.e-4, 5.e-4, 1.e-3, 2.e-3, 5.e-3, \
    #            1.e-2, 2.e-2, 5.e-2, 1.e-1, 2.e-1, 5.e-1, 1.e0], \
    #            [r"$1\times 10^{-4}$",r"$2\times 10^{-4}$",r"$5\times 10^{-4}$",\
    #             r"$1\times 10^{-3}$",r"$2\times 10^{-3}$",r"$5\times 10^{-3}$",\
    #             r"$1\times 10^{-2}$",r"$2\times 10^{-2}$",r"$5\times 10^{-2}$",\
    #             r"$1\times 10^{-1}$",r"$2\times 10^{-1}$",r"$5\times 10^{-1}$",r"$1\times 10^{0}$"])
    plt.yticks([0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, \
                0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1], \
               [0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.005, \
                0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1])
    plt.plot([time/(2.*PI) for time in time_v], ecc_rms_v, \
             linewidth=1, \
             label=r"$\langle e^2\rangle^{1/2}$", color='b', ls='-')
    plt.plot([time/(2.*PI) for time in time_v], inc_rms_v, \
             linewidth=1, \
             label=r"$\langle i^2\rangle^{1/2}$", color='r', ls='--')
    plt.ylabel(r"$\langle e^2\rangle^{1/2}, \langle i^2\rangle^{1/2}$", fontsize=25)
    plt.xlabel(r"Time (year)", fontsize=25)
    plt.legend(fontsize=17, loc='upper left', frameon=False)
    plt.xlim([min(time_v)/(2.*PI),max(time_v)/(2.*PI)])
    plt.ylim([min(ecc_rms_v+inc_rms_v)/2,max(ecc_rms_v+inc_rms_v)*2])
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    #plt.subplots_adjust(bottom=0.13, left=0.2)
    plt.subplots_adjust(bottom=0.13, left=0.18)
    
    pngfile = figure_dir + "/eccincEvolution.png"
    plt.savefig(pngfile, dpi=100)
    plt.close()
    #"""

    return


if __name__ == '__main__':
    main()
