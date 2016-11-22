#This program uses files with chromosomes sorted by dictionary order
import sys, os, string, random

Gwidth = 100 # of bins
Twidth = (Gwidth * 2)
bp = 1000

def splitline(rrpp):
    rrrppp = list(rrpp.split())
    return [str(rrrppp[0]), int(rrrppp[1]), int(rrrppp[2]), str(rrrppp[3]), float(rrrppp[4]), str(rrrppp[5])]

def iofile(rfile, pfile, signal):
     #//////////////////////////////////////////////////////////////////////////
     #/////////                                                      ///////////
     #/////////                                                      ///////////
     #/////////     a-Gwidth----a-1|a------a+Gwidth-1  -> a-Gwidth, a+ Gwidth   ///////////
     #/////////     b-Gwidth----b-1|b------b+Gwidth-1                ///////////
     #//////////////////////////////////////////////////////////////////////////
         
    clr = 0    #count the lines checked in region file
    clp = 0    #count the lines checked in position file

    for round in [0, 1]:
        rf = open(rfile, "r")
        pf = open(pfile, "r")
        clr = 0
        clp = 0
        rr = rf.readline()
        pp = pf.readline()
        if rr in ['','\n']:
            print "rf end of line"
            break
        if pp in ['','\n']:
            print "pf end of line"
            break

        clr += 1
        clp += 1
        
        (rchr, ra, rb, rid, rscore, rstrand) = splitline(rr)
        (pchr, pa, pb, pid, pscore, pstrand) = splitline(pp)
        if rstrand not in ['-', '+'] or pstrand not in ['-', '+']:
            print("r: " + rstrand + " p: " + pstrand)

        while True:
            if rchr < pchr or (rchr == pchr and [ra, rb][round] + bp <= pa):  #if the position is on the right of the region
                rr = rf.readline()
                if rr in ['','\n']:
                    #print "rf end of line"
                    break
                else:
                    (rchr, ra, rb, rid, rscore, rstrand) = splitline(rr)
                    if rstrand not in ['-', '+'] or pstrand not in ['-', '+']:
                        print("r: " + rstrand + " p: " + pstrand)
                    clr += 1
								
            else: #if the position falls in or locates on the left of the region
                if rchr == pchr and max(pa, [ra, rb][round] - bp) < min(pb, [ra, rb][round] + bp): #if the position falls in
                    for pos in range(max(pa, [ra, rb][round] - bp), min(pb, [ra, rb][round] + bp)):
                        LOCTN = int((pos - [ra, rb][round] + bp + 0.5) * Twidth/(2 * bp))
                        if rstrand == ['+', '-'][round]:
                            if rstrand == '+':
                                signal[rstrand == pstrand][LOCTN] += 1 #increase the corresponding bin's signal strength
                            elif rstrand == '-':
                                signal[rstrand == pstrand][Twidth - 1 - LOCTN] += 1
                            if rstrand not in ['-', '+'] or pstrand not in ['-', '+']:
                                print("r: " + rstrand + " p: " + pstrand)
										
                pp = pf.readline()
                if pp in ['','\n']:
                    #print "pf end of line"
                    break
                else:
                    (pchr, pa, pb, pid, pscore, pstrand) = splitline(pp)
                    clp += 1
        rf.close()
        pf.close()
        print("went through " + str(clr) + rfile + str(clp) + pfile)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("\nThis program plots metaplot of the signal of pf file over rf regions")
        print("Usage: python\t" + sys.argv[0] + "\tregionfile positionfile")
        exit(1)

    signal = [[0] * Twidth, [0] * Twidth]   #define an array. We can't write [[0.0] * Twidth]*2 which links the two subarray
    distri = iofile(sys.argv[1], sys.argv[2], signal) #for each combination of region file and position file, plot

    oo1 = open(sys.argv[1] + sys.argv[2] + "S", "w")
    oo0 = open(sys.argv[1] + sys.argv[2] + "AS", "w")
    for i in range(Twidth):
        outs1 = str(i) + "\t" + str(signal[1][i]) + "\n"
        outs0 = str(i) + "\t" + str(signal[0][i]) + "\n"
        oo1.write(outs1)
        oo0.write(outs0)
    oo1.close();
    oo0.close();
