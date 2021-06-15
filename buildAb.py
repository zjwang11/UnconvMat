import os, sys, time, getopt 
import importlib
import numpy as np 
from copy import deepcopy

var = '??'
bilbaodata = '/home/jcgao/work/BilBaoData/'

def read_SiteK(nsg):
    filename = bilbaodata + "BRlist2_A/SiteK_"+str(nsg)+".cht"
    f = open(filename, "r")
    Nsite = int(f.readline().strip().split()[-1])
    Nktp  = int(f.readline().strip().split()[-1])
    name_siteg     = ['' for i in range(Nsite)]
    name_irr_siteg = [[] for i in range(Nsite)]
    name_irr_littg = [[] for i in range(Nktp)]
    f.readline()

    BRdata = [[[] for j in range(Nktp)] for i in range(Nsite)]
    for isite in range(Nsite):
        for ikvec in range(Nktp):
            tmp = f.readline().strip().split()
            name_siteg[isite] = tmp[-2]
            nirr_site = int(tmp[-1])
            nirr_kvec = int(f.readline().strip().split()[-1])
            name_irr_siteg[isite] = ['' for i in range(nirr_site)]
            name_irr_littg[ikvec] = ['' for i in range(nirr_kvec)]

            tmp1= f.readline()
            tmp = tmp1.strip().split()
            for i in range(nirr_kvec):
                name_irr_littg[ikvec][i] = tmp1[7+5*i:12+5*i].strip()
                #name_irr_littg[ikvec][i] = tmp[2+i]
            
            for irr_site in range(nirr_site):
                tmp = f.readline().strip().split()
                name_irr_siteg[isite][irr_site] = tmp[0]
                data = list(map(int,tmp[1:]))
                BRdata[isite][ikvec].append(data)
            f.readline()
    f.close()
    return BRdata, name_siteg, name_irr_siteg, name_irr_littg


def writecode(A,b,lowbound,uppbound,isEBR):
    num_var = len(A[0])
    var_str = ['a'+str(i) for i in range(num_var)]
    with open('runZ3.py','w') as f:
        f.write("from z3 import *\n")
        f.write("\n")
        f.write("def decomBR():\n")
        f.write("    s = Solver()\n")
        for i in range(num_var):
            f.write("    a{} = Int('a{}')\n".format(i,i))
        f.write("    vars = [")
        for i in range(num_var):
            f.write("a{}".format(i))
            if i!=num_var-1:
                f.write(",")
        f.write("]\n")

        f.write("\n")
        for i in range(num_var):
            f.write("    s.add(a{} >= {})\n".format(i,lowbound[i]))
            if not isEBR:
                f.write("    s.add(a{} <= {})\n".format(i,uppbound[i]))
        for ieq, eq in enumerate(A):
            tmp_str = ''
            for icoeff, coeff in enumerate(eq):
                tmp_str = tmp_str + "{}*a{}+".format(coeff,icoeff)
            tmp_str = tmp_str + "0 == {}".format(b[ieq])
            f.write("    s.add("+tmp_str+")\n")
        #f.write(f"solve({tmp_str})")
        f.write("\n")
        f.write("    alldata = []\n")
        f.write("    #f = open('decompose','a+')\n")
        f.write("    while s.check() == sat:\n")
        f.write("        m = s.model()\n")
        f.write("        if not m:\n")
        #f.write("            print('No solutions any more')\n")
        f.write("            break\n")
        #f.write("        print(m)\n")
        #f.write("    tmp_str = '          '.join([str(m[i]) for i in vars])\n")
        #f.write("    f.write('    ')\n")
        f.write("        tmp_str = '   '.join([str(m[i]) for i in vars])\n")
        f.write("        data = list(map(int,tmp_str.strip().split()))\n")
        f.write("        alldata.append(data)\n")
        f.write("        if len(alldata) > 100:\n")
        f.write("            break\n")
        f.write("        #for i in range(len(data)):\n")
        f.write("        #    f.write('{:>7d}'.format(data[i]))\n")
        f.write("        s.add(Not(And([v() == m[v] for v in m])))\n")
        #f.write("    else:\n")
        #f.write("        print('No solutions any more')\n")
        f.write("    #f.close()\n")
        f.write("\n")
        f.write("    return alldata\n")
        f.write("\n")
        f.write("\n")
        f.write("if __name__=='__main__':\n")
        f.write("    alldata = decomBR()\n")


def buildAb_oneshot(inputfile, BRdata, soc=0):

    isEBR = 0
    
    BR = []
    BRnum = []
    BRset = []
    BRlabel = []
    BRlabelset = []
    with open("atomwp.in", "r") as f:
        isEBR = int(f.readline().strip())
        nBR = int(f.readline().strip())
        for i in range(nBR):
            line = f.readline().strip().split()
            BR.append(tuple(map(int, line[0:2])))
            BRlabel.append(line[2])
        for idx, iBR in enumerate(BR):
            if iBR not in BRset:
                BRset.append(iBR)
                BRlabelset.append(BRlabel[idx])
        for i in BRset:
            BRnum.append(BR.count(i))

    maxk = []
    irre_onk = []
    print(inputfile)
    with open(inputfile,"r") as f:
        (nsg, nmaxk, nband) = tuple(map(int,f.readline().strip().split()))
        for i in range(nmaxk):
            data = list(map(int,f.readline().strip().split()))
            maxk.append(data[0])
            irre_onk.append(data[1:])

    nirr_eachk = [] 
    kname = []
    with open(bilbaodata + "kvec_list_A.txt","r") as f:
        while True:
            (isg, numk) = tuple(map(int,f.readline().strip().split()))
            if (isg==nsg):
                for i in range(numk):
                    tmp = f.readline().strip().split()
                    nirr_eachk.append(tuple(map(int,tmp[0:4])))
                    kname.append(tmp[5])
                break
            else:
                for i in range(numk):
                    f.readline()

    for ik in range(len(maxk)):
        if maxk[ik] < 0:
            absk = abs(maxk[ik])
            abskname = kname[absk-1]
            knameA = abskname+'A'
            kaindex = kname.index(knameA) + 1
            maxk[ik] = kaindex 
    
    nkir = []
    for ik in maxk:
        for jdata in nirr_eachk:
            if ik==jdata[0]:
                if soc==0:
                    nkir.append(jdata[2])
                else:
                    nkir.append(jdata[3])

    #BRdata, name_siteg, name_irr_siteg, name_irr_littg  = read_SiteK(nsg)

    # construct A matrix
    A = np.zeros((sum(nkir), len(BRset)),dtype=int)

    for iw in range(len(BRset)):
        isite = BRset[iw][0]
        iirre = BRset[iw][1]
        icol = iw 
        irow = 0
        for ik in range(len(maxk)):
            for ir in range(nkir[ik]):
                irk = ir if soc==0 else (ir-nkir[ik])
                A[irow, icol] = A[irow, icol] + BRdata[isite-1][maxk[ik]-1][iirre-1][irk]
                irow = irow + 1

    # construct b array
    b = np.array([0 for i in range(sum(nkir))])
    for ik in range(len(maxk)):
        for rep in irre_onk[ik]:
            if soc!=0: rep = rep-nirr_eachk[maxk[ik]-1][2]
            b[sum(nkir[:ik])+rep-1] += 1
           
    lowbound = np.array([0 for i in range(len(A[0]))])
    uppbound = np.array(BRnum)

    writecode(A,b,lowbound,uppbound,isEBR)
    
    #from runZ3 import decomBR
    import runZ3
    if (sys.version_info.major == 3):
        importlib.reload(runZ3)
    elif (sys.version_info.major == 2):
        reload(runZ3)
    time.sleep(1)
    alldata = runZ3.decomBR()
    del runZ3

    num_decom = len(alldata)
    print('Number of solutions:', num_decom)
    alldata = np.array(alldata).T 

    #with open('decompose','w') as f:
    #    f.write("#")
    #    for idx, ibr in enumerate(BRset):
    #        f.write(" {:>2d}@{:<2d}({:>2d}) ".format(ibr[1],ibr[0],uppbound[idx]))
    #    f.write("\n")
    dataformat = "{:>7d}"
    with open('decompose_'+inputfile,'w') as f:
        f.write('{:>7d}\n'.format(num_decom))
        for idx, ibr in enumerate(BRset):
            f.write("{:>4d}".format(idx+1)+
                    "  {:>2d}@{:<2d}".format(ibr[1],ibr[0])+
                    "{:>10}".format(BRlabelset[idx])+
                    "    ({:3d}) :".format(uppbound[idx]))
            if len(alldata)!=0:
                for jdx, jbr in enumerate(alldata[idx]):
                    f.write("{:>3d};".format(alldata[idx][jdx]))
            f.write("\n")
        #f.write("# ")
        #for idx in range(len(BRset)):
        #    f.write(" ({:3d}) ".format(uppbound[idx]))
        #f.write("\n")

    return num_decom 


def check_input(inputfile):
    num_vars = 0
    f = open(inputfile,'r')
    (nsg, numk, nband) = tuple(map(int, f.readline().strip().split()))

    k_withvar = []
    for i in range(numk):
        data = f.readline().strip().split()
        for d in data:
            if d == var:
                num_vars += 1
                k_withvar.append(int(data[0]))
    return nsg, num_vars, k_withvar 


def run_buildAb(inputfile='fort.67', soc=0):

    nsg, num_vars, k_withvar = check_input(inputfile)
    BRdata, name_siteg, name_irr_siteg, name_irr_littg = read_SiteK(nsg)
    nirr_eachk = [] 
    with open(bilbaodata + "kvec_list_A.txt","r") as f:
        while True:
            (isg, numk) = tuple(map(int,f.readline().strip().split()))
            if (isg==nsg):
                for i in range(numk):
                    nirr_eachk.append(tuple(map(int,f.readline().strip().split()[0:4])))
                break
            else:
                for i in range(numk):
                    f.readline()

    irr_trial = []
    for ik in range(len(k_withvar)):
        if (soc == 0):
            irr_trial.append([i+1 for i in range(nirr_eachk[k_withvar[ik]-1][2])])
        else:
            irr_trial.append([nirr_eachk[k_withvar[ik]-1][2]+i+1 for i in range(nirr_eachk[k_withvar[ik]-1][3])])

    if (num_vars == 0):
        num_decom = buildAb_oneshot(inputfile, BRdata, soc)
    elif (num_vars == 1):
        summary = open('solutions','w')
        summary.write('# trial_irr  num_solutions\n')
        for irr1 in irr_trial[0]:
            new_inputfile = write_input(inputfile,[irr1])
            num_decom = buildAb_oneshot(new_inputfile, BRdata, soc)
            if (num_decom > 0):
                summary.write('{:>5d}:{:>5d}\n'.format(irr1, num_decom))
            else:
                summary.write('{:>5d}:\n'.format(irr1))
        summary.close()
    elif (num_vars == 2):
        summary = open('solutions','w')
        summary.write('# trial_irr1  trial_irr2  num_solutions\n')
        for irr1 in irr_trial[0]:
            for irr2 in irr_trial[1]:
                new_inputfile = write_input(inputfile,[irr1, irr2])
                num_decom = buildAb_oneshot(new_inputfile, BRdata, soc)
                if (num_decom > 0):
                    summary.write('{:>5d}{:>5d}:{:>5d}\n'.format(irr1,irr2, num_decom))
                else:
                    summary.write('{:>5d}{:>5d}:\n'.format(irr1,irr2))
        summary.close()
    elif (num_vars == 3):
        summary = open('solutions','w')
        summary.write('# trial_irr1  trial_irr2  trial_irr3  num_solutions\n')
        for irr1 in irr_trial[0]:
            for irr2 in irr_trial[1]:
                for irr3 in irr_trial[2]:
                    new_inputfile = write_input(inputfile,[irr1, irr2, irr3])
                    num_decom = buildAb_oneshot(new_inputfile, BRdata, soc)
                    if (num_decom > 0):
                        summary.write('{:>5d}{:>5d}{:>5d}:{:>5d}\n'.format(irr1,irr2,irr3, num_decom))
                    else:
                        summary.write('{:>5d}{:>5d}{:>5d}:\n'.format(irr1,irr2,irr3))
        summary.close()
    else:
        sys.exit('{}(>3) variables in input file'.format(num_vars))
    return


def write_input(inputfile, irr_trial):
    newfile = inputfile
    for irr in irr_trial:
        newfile = newfile + '_' + str(irr)
    with open(inputfile, 'r') as f, open(newfile, 'w') as w:
        tmp = f.readline()
        numk = int(tmp.split()[1])
        w.write(tmp)
        for ik in range(numk):
            tmp = f.readline().strip().split()
            for d in tmp:
                if d == var:
                    w.write('{:>3d}'.format(irr_trial.pop(0)))
                else:
                    w.write('{:>3s}'.format(d))
            w.write('\n')
    return newfile 


def run_rsi(decompose_ebr='decompose.601', brdatafile='fort.67', soc=0, trs=1):

    wfile = open('rsidelta','w')

    replist = np.load(bilbaodata+'pgirr.npy')
    wplist = np.load(bilbaodata+'wp_list.npy')

    nsg, num_vars, k_withvar = check_input(brdatafile)
    if soc == 0 and trs == 1:
        data = np.load(bilbaodata+'src_rsi/nosoc_trs_code/'+str(nsg)+'.npy').item()

    f = open(decompose_ebr, 'r')
    fdata = f.readlines()
    f.close()
    num_decom = int(fdata[0].strip())
    if num_decom == 0:
        sys.exit('No decomposition in the file '+decompose_ebr)
        
    decomT = []
    ebrname = []
    ebrname_str = []
    ebrname_ind = []
    for line in fdata[1:]:
        tmp0 = line.strip().split(':')[0]
        tmp = line.strip().split(':')[1]
        decomT.append(list(map(int, tmp.strip().split(';')[:-1])))
        ebrname.append(tmp0.strip().split()[1])
        ebrname_str.append(tmp0.strip().split()[2])
        ebrname_ind.append(tmp0.strip().split()[0])
    decomT = np.array(decomT)
    decom = decomT.T

    ebrname_GM = []
    for ebr in ebrname_str:
        nameirr = ebr.split('@')[0]
        w = ebr.split('@')[1]
        pgname = wplist[nsg][w]
        for pg in replist[1:]:
            if pg['name'] == pgname:
                break 
        for ind, irname in enumerate(pg['rep3']):
            if irname == nameirr:
                break
        ebrname_GM.append(pg['rep1'][ind]+'@'+w)

    for w in data.keys():
        data[w]['row_A'] = []

        pgname = wplist[nsg][w]
        for pg in replist[1:]:
            if pg['name'] == pgname:
                break 
        pgtmp = deepcopy(pg)
        pgtmp['rep1'].reverse()
        pgtmp['rep3'].reverse()
        GMlist = pgtmp['rep1']
        Alist = pgtmp['rep3']
        for rep in data[w]['row']:
            reptmp = rep
            for repGM, repA in zip(GMlist, Alist):
                reptmp = reptmp.replace(repGM, repA)
            data[w]['row_A'].append(reptmp)

    existwyck = []
    for ebr in ebrname_GM:
        w = ebr.split('@')[1]
        if w not in existwyck:
            existwyck.append(w)

    wfile.write('{:>4d}\n'.format(len(existwyck)))

    for w in existwyck:
        wfile.write('{:>5s}'.format(w))
        if w not in data.keys():
            sys.exit('wyck in essebr {:>5s} not in rsi formula'.format(w))
        lenrow = len(data[w]['row'])
        lencol = len(data[w]['col'])

        constraint = []
        ebrind = []
        for ebr in data[w]['row']:
            if '&' in ebr:
                e1 = ebr.split('&')[0]
                e2 = ebr.split('&')[1]
                ebrind.append(ebrname_GM.index(e1))
                constraint.append([ebrname_GM.index(e1),ebrname_GM.index(e2)])
            else:
                ebrind.append(ebrname_GM.index(ebr))
        
        decomW = np.zeros((lenrow, num_decom), dtype=int)
        for isol in range(num_decom):
            for constr in constraint:
                c1 = constr[0]
                c2 = constr[1]
                if decomT[c1,isol] != decomT[c2,isol]:
                    decomW[:,isol] = 0
                    break 
            else:
                for ii, ind in enumerate(ebrind):
                    decomW[ii,isol] = decomT[ind,isol]

        delta = np.dot(data[w]['L'], decomW)
        
        module = []
        for irank in range(lenrow):
            if irank < min(lenrow, lencol):
                module.append(data[w]['ind'][irank,irank])
            else:
                module.append(0)
        
        wfile.write('{:5d}{:5d}{:5d}\n'.format(lenrow, num_decom, len(module)-module.count(1)))
        for imod in module:
            wfile.write('{:>4d}'.format(imod))
        wfile.write('\n')

        # rsi formula
        for ind, imod in enumerate(module):
            if imod != 1:
                for icol in range(lenrow):
                    if data[w]['L'][ind, icol] != 0:
                        ebrname_A = data[w]['row_A'][icol]
                        if   data[w]['L'][ind, icol] == 1:
                            wfile.write('+')
                        elif data[w]['L'][ind, icol] ==-1:
                            wfile.write('-')
                        elif data[w]['L'][ind, icol] >  0:
                            wfile.write('+'+str(data[w]['L'][ind,icol]))
                        else:
                            wfile.write(str(data[w]['L'][ind,icol]))
                        wfile.write(r'n('+ebrname_A+')')
                wfile.write('\n')

        for irow in range(lenrow):
            for icol in range(num_decom):
                wfile.write('{:>3d}'.format(delta[irow,icol]))
            wfile.write('\n')
    wfile.close()    
    
    return 


def run_essebr(decompose_ebr='decompose.601', decompose_abr='decompose.603',brdatafile='fort.67', soc=0):

    f = open(decompose_ebr, 'r')
    fdata = f.readlines()
    f.close()
    num_decom = int(fdata[0].strip())
    if num_decom == 0:
        sys.exit('No decomposition in the file '+decompose_ebr)
        
    decom = []
    ebrname = []
    ebrname_str = []
    ebrname_ind = []
    for line in fdata[1:]:
        tmp0 = line.strip().split(':')[0]
        tmp = line.strip().split(':')[1]
        decom.append(list(map(int, tmp.strip().split(';')[:-1])))
        ebrname.append(tmp0.strip().split()[1])
        ebrname_str.append(tmp0.strip().split()[2])
        ebrname_ind.append(tmp0.strip().split()[0])

    f = open(decompose_abr, 'r')
    fdata = f.readlines()
    f.close()
    num_decom = int(fdata[0].strip())
    if num_decom >  0:
        sys.exit('There are decompositions in the file '+decompose_abr)
    
    essebr_list = []
    essebr_str_list = []
    essebr_ind_list = []
    num_decom = len(decom[0])
    num_ebr   = len(decom)
    #for ie in range(num_ebr):
    #    for n in range(num_decom):
    #        if decom[ie][n] == 0:
    #            break 
    #    else:
    #        essebr_list.append(ebrname[ie])
    #        essebr_str_list.append(ebrname_str[ie])
    if len(essebr_list) == 0:
        essebr_list = ebrname
        essebr_str_list = ebrname_str
        essebr_ind_list = ebrname_ind
    print(essebr_str_list)

    nsg, num_vars, k_withvar = check_input(brdatafile)
    BRdata, name_siteg, name_irr_siteg, name_irr_littg = read_SiteK(nsg)

    decom_collect = []
    essw = open('essebr','w')
    for ind, essebr in enumerate(essebr_list):
        new_inputfile = write_inputrmebr(brdatafile, essebr, BRdata)
        if new_inputfile is None:
            continue
        num_decom = buildAb_oneshot(new_inputfile, BRdata, soc)
        decom_collect.append(num_decom)
        essw.write('{:>4s} {:>10s} {:>10s} {:>10d}\n'.format(
                          essebr_ind_list[ind], essebr, essebr_str_list[ind], num_decom))
    essw.close()
    
    # can not find abr when just remove one ebr
    if sum(decom_collect) == 0:
        essw2 = open('essebr2','w')
        for ind1, essebr1 in enumerate(essebr_list):
            new_inputfile = write_inputrmebr(brdatafile, essebr1, BRdata)
            if new_inputfile is None:
                continue
            for ind2 in range(ind1, len(essebr_list)):
                essebr2 = essebr_list[ind2]
                new_inputfile2 = write_inputrmebr(new_inputfile, essebr2, BRdata)
                if new_inputfile2 is None:
                    continue
                num_decom = buildAb_oneshot(new_inputfile2, BRdata, soc)
                essw2.write('{:>4s} {:>10s} {:>10s} {:>4s} {:>10s} {:>10s} {:>10d}\n'.format(
                    essebr_ind_list[ind1], essebr1, essebr_str_list[ind1], 
                    essebr_ind_list[ind2], essebr2, essebr_str_list[ind2], num_decom))
        essw2.close()
    
    # can not find abr when just remove one and two ebr
    #if sum(decom_collect) == 0:
    #    essw3 = open('essebr3','w')
    #    for ind1, essebr1 in enumerate(essebr_list):
    #        new_inputfile = write_inputrmebr(brdatafile, essebr1, BRdata)
    #        if new_inputfile is None:
    #            continue
    #        for ind2 in range(ind1, len(essebr_list)):
    #            essebr2 = essebr_list[ind2]
    #            new_inputfile2 = write_inputrmebr(new_inputfile, essebr2, BRdata)
    #            if new_inputfile2 is None:
    #                continue
    #            for ind3 in range(ind2, len(essebr_list)):
    #                essebr3 = essebr_list[ind3]
    #                new_inputfile3 = write_inputrmebr(new_inputfile, essebr3, BRdata)
    #                if new_inputfile3 is None:
    #                    continue
    #                num_decom = buildAb_oneshot(new_inputfile3, BRdata, soc)
    #                essw3.write('{:>4s} {:>10s} {:>10s} {:>4s} {:>10s} {:>10s} {:>4s} {:>10s} {:>10s} {:>10d}\n'
    #                            .format(essebr_ind_list[ind1], essebr1, essebr_str_list[ind1], 
    #                                    essebr_ind_list[ind2], essebr2, essebr_str_list[ind2], 
    #                                    essebr_ind_list[ind3], essebr3, essebr_str_list[ind3], num_decom))
    #    essw3.close()
    #
    ## can not find abr when just remove one and two and three ebr
    #if sum(decom_collect) == 0:
    #    essw4 = open('essebr4','w')
    #    for ind1, essebr1 in enumerate(essebr_list):
    #        new_inputfile = write_inputrmebr(brdatafile, essebr1, BRdata)
    #        if new_inputfile is None:
    #            continue
    #        for ind2 in range(ind1, len(essebr_list)):
    #            essebr2 = essebr_list[ind2]
    #            new_inputfile2 = write_inputrmebr(new_inputfile, essebr2, BRdata)
    #            if new_inputfile2 is None:
    #                continue
    #            for ind3 in range(ind2, len(essebr_list)):
    #                essebr3 = essebr_list[ind3]
    #                new_inputfile3 = write_inputrmebr(new_inputfile, essebr3, BRdata)
    #                if new_inputfile3 is None:
    #                    continue
    #                for ind4 in range(ind3, len(essebr_list)):
    #                    essebr4 = essebr_list[ind4]
    #                    new_inputfile4 = write_inputrmebr(new_inputfile, essebr4, BRdata)
    #                    if new_inputfile4 is None:
    #                        continue
    #                    num_decom = buildAb_oneshot(new_inputfile3, BRdata, soc)
    #                    essw3.write('{:>4s} {:>10s} {:>10s} {:>4s} {:>10s} {:>10s} {:>4s} {:>10s} {:>10s} {:>4s} {:>10s} {:>10s} {:>10d}\n'
    #                                .format(essebr_ind_list[ind1], essebr1, essebr_str_list[ind1], 
    #                                        essebr_ind_list[ind2], essebr2, essebr_str_list[ind2], 
    #                                        essebr_ind_list[ind3], essebr3, essebr_str_list[ind3], 
    #                                        essebr_ind_list[ind4], essebr4, essebr_str_list[ind4], num_decom))
    #    essw4.close()

    return

def write_inputrmebr(inputfile, essebr, BRdata):

    ess_site = int(essebr.strip().split('@')[1])
    ess_irr  = int(essebr.strip().split('@')[0])
    removedict = dict()
    for ik in range(len(BRdata[0])):
        removedict[ik+1] = []
        for ind, ir in enumerate(BRdata[ess_site-1][ik][ess_irr-1]):
            for i in range(ir):
                removedict[ik+1].append(ind+1)

    maxklist = []
    newfile = inputfile
    newfile = newfile + '_' + essebr.replace('@','a')
    with open(inputfile, 'r') as f, open(newfile, 'w') as w:
        tmp = f.readline()
        numk = int(tmp.split()[1])
        w.write(tmp)
        for ik in range(numk):
            tmp = list(map(int, f.readline().strip().split()))
            tmpir = tmp[1:]
            indk = tmp[0]
            maxklist.append(indk)
            for ir in removedict[indk]:
                try:
                    tmpir.remove(ir)
                except ValueError:
                    return None
            w.write('{:>3d}'.format(indk))
            for r in tmpir:
                w.write('{:>3d}'.format(r))
            w.write('\n')
        w.write('\n')
        w.write(essebr+'\n')
        for key in removedict:
            if key in maxklist:
                w.write('{:>3d}'.format(key))
                for r in removedict[key]:
                    w.write('{:>3d}'.format(r))
                w.write('\n')
    return newfile 

def helpfunction():
    print("""
    python buildAb.py [-i inputfile] [--soc=1] [--ess] [-e decomposefile_ebr] [-a decomposefile_abr]
    """)
    return

if __name__ == '__main__':

    work = 'decom'
    inputfile = 'fort.67'
    decompose_ebr = 'decompose.601'
    decompose_abr = 'decompose.603'
    soc = 0
    trs = 1

    argv = sys.argv[1:]
    try:
        opts, args = getopt.getopt(argv, 'hi:e:a:', ['help','ess','rsi','soc=','trs='])
    except getopt.GetoptError:
        helpfunction()
        sys.exit()
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            helpfunction()
        elif opt == '--ess':
            work = 'ess'
        elif opt == '--rsi':
            work = 'rsi'
        elif opt == '-i':
            inputfile = arg
        elif opt == '-e':
            decompose_ebr = arg
        elif opt == '-a':
            decompose_abr = arg
        elif opt == '--soc':
            soc = int(arg)
        elif opt == '--trs':
            trs = int(arg)

    if work == 'decom':
        run_buildAb(inputfile, soc)
    elif work == 'ess':
        run_essebr(decompose_ebr, decompose_abr, inputfile, soc)
    elif work == 'rsi':
        run_rsi(decompose_ebr, inputfile, soc, trs)
