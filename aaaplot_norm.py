#This program uses files with chromosomes sorted by dictionary order
import sys, os, string, random

Gwidth = 100; # of bins for gene region
Twidth = (Gwidth * 3)

def splitline(rrpp):
    rrrppp = list(rrpp.split())
    return [str(rrrppp[0]), int(rrrppp[1]), int(rrrppp[2]), str(rrrppp[3]), float(rrrppp[4]), str(rrrppp[5])]

def iofile(rfile, pfile, signal, coverage, signalcoverage, olist):
     #//////////////////////////////////////////////////////////////////////////
     #/////////                                                      ///////////
     #/////////             left    interesting   right              ///////////
     #/////////      2*a-b ----a-1|a-----------b-1|b-------2*b-a-1   ///////////
     #/////////                                                      ///////////
     #//////////////////////////////////////////////////////////////////////////
         
    clr = 0    #count the lines checked in region file
    clp = 0    #count the lines checked in position file

    for round in [-1, 0, 1]:
        rf = open(rfile, "r")
        pf = open(pfile, "r")
        coverage[round + 1] = 0
        clr = 0
        clp = 0
        rr = rf.readline()
        pp = pf.readline()
        if rr in ['','\n']:
            print "rf end of line"
            break
        if pp in ['','\n']:
            print "pf end of line"
            if round == 0:
                while rr not in ['','\n']:
                    rr = rf.readline()
                    olist.append(0)
            break

        clr += 1
        clp += 1
        
        (rchr, ra, rb, rid, rscore, rstrand) = splitline(rr)
        (pchr, pa, pb, pid, pscore, pstrand) = splitline(pp)
        if rstrand not in ['-', '+'] or pstrand not in ['-', '+']:
            print("r: " + rstrand + " p: " + pstrand)
        
        coverage[round + 1] += rb - ra
        covered = [0] * (rb - ra)

        while True:
            if rchr < pchr or (rchr == pchr and rb + round * (rb - ra)<= pa):  #if the position is on the right of the region
                rr = rf.readline()
                if round == 0:
                    olist.append(covered.count(1))
                signalcoverage[round + 1] += covered.count(1)
                if rr in ['','\n']:
                    #print "rf end of line"
                    break
                else:
                    (rchr, ra, rb, rid, rscore, rstrand) = splitline(rr)
                    if rstrand not in ['-', '+'] or pstrand not in ['-', '+']:
                        print("r: " + rstrand + " p: " + pstrand)
            
                    coverage[round + 1] += rb - ra
                    covered = [0] * (rb - ra)
                    clr += 1
								
            else: #if the position falls in or locates on the left of the region
                if rchr == pchr and max(pa, ra + round * (rb - ra)) < min(pb, rb + round * (rb - ra)): #if the position falls in
                    for pos in range(max(pa, ra + round * (rb - ra)), min(pb, rb + round * (rb - ra))):
                        LOCTN = int((pos - 2 * ra + rb + 0.5) * Gwidth/(rb - ra))
                        covered[pos - (ra + round * (rb - ra))] = 1
                        if rstrand == '+':
                            signal[rstrand == pstrand][LOCTN] += 1 #increase the corresponding bin's signal strength
                        elif rstrand == '-':
                            signal[rstrand == pstrand][Twidth - 1 - LOCTN] += 1
                            if rstrand not in ['-', '+'] or pstrand not in ['-', '+']:
                                print("r: " + rstrand + " p: " + pstrand)
										
                pp = pf.readline()
                if pp in ['','\n']:
                    #print "pf end of line"
                    if round == 0:
                        olist.append(covered.count(1))
                        rr = rf.readline()
                    signalcoverage[round + 1] += covered.count(1)
                    if round == 0:
                        while rr not in ['','\n']:
                            rr = rf.readline()
                            olist.append(0)
                    break
                else:
                    (pchr, pa, pb, pid, pscore, pstrand) = splitline(pp)
                    
                    clp += 1
        rf.close()
        pf.close()

        #print("went through " + str(clr) + rfile + str(clp) + pfile)
        if round == 0:
            print(str(sys.argv[1]) + "\t" + str(signalcoverage[round + 1]) + " of " + str(coverage[round + 1]) + " bp are covered")
    return [olist, signalcoverage[1]]

def permtest(olist, signalcoverage):
    distrilist151 = []
    distrilist243 = []
    for i in range(1000):
        random.shuffle(olist)
        distrilist151.append(sum(olist[0:92]))
        distrilist243.append(sum(olist[0:243]))
    distrilist151.append(signalcoverage)
    distrilist243.append(signalcoverage)
    print signalcoverage
    print (sorted(distrilist151)).index(signalcoverage), (sorted(distrilist151)).index(signalcoverage)



if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("\nThis program plots metaplot of the signal of pf file over rf regions")
        print("Usage: python\t" + sys.argv[0] + "\tregionfile positionfile")
        exit(1)

    coverage = [0, 0, 0]
    signalcoverage = [0, 0, 0]
    signal = [[0] * Twidth, [0] * Twidth]   #define an array. We can't write [[0.0] * Twidth]*2 which links the two subarray
    [distri, SC] = iofile(sys.argv[1], sys.argv[2], signal, coverage, signalcoverage, []) #for each combination of region file and position file, plot
#    permtest(distri, SC)
    oo1 = open(sys.argv[1] + sys.argv[2] + "S", "w")
    oo0 = open(sys.argv[1] + sys.argv[2] + "AS", "w")
    for i in range(Twidth):
        outs1 = str(i) + "\t" + str(signal[1][i]) + "\n"
        outs0 = str(i) + "\t" + str(signal[0][i]) + "\n"
        oo1.write(outs1)
        oo0.write(outs0)
    oo1.write(str(signalcoverage[1]) + "\t" + str(coverage[1]) + "\n")
    oo1.close();
    oo0.close();
