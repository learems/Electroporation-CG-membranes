'''
Copyright (c) 2019, Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory.
Written by Harsh Bhatia (hbhatia@llnl.gov) and Peer-Timo Bremer (bhatia4@llnl.gov)
LLNL-CODE-763493. All rights reserved.

This file is part of MemSurfer, Version 1.0.
Released under GNU General Public License 3.0.
For details, see https://github.com/LLNL/MemSurfer.
'''

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
import os, sys
import numpy as np
import argparse
import logging
LOGGER = logging.getLogger(__name__)

import memsurfer
from memsurfer import utils

import MDAnalysis
from MDAnalysis.analysis.distances import distance_array
from MDAnalysis.analysis.leaflet import LeafletFinder, optimize_cutoff

import mdreader
#from lipidType import *

ddir =  sys.argv[1]
PM = sys.argv[2]
mem = sys.argv[3]
grofile = sys.argv[4]
xtcfile = sys.argv[5]
folder = sys.argv[6]

if ("APM" in PM):
    from lipidType_APM import *
elif ("BPM" in PM):
    from lipidType_BPM import *


# -------------------------------------------------------------------------------
# CUSTOM FUNCTIONS for extracting bilayer properties
# -------------------------------------------------------------------------------
# Function for computing cosines of the angles between lipid tail bonds and membrane normal
def get_CosTailAngles(tp,mt,syst,lipidTypeDic):
    # Replace cosTA with P2t if you want to output order parameter directly
    ### P2t = np.empty((len(tp),12))
    ### P2t[:] = np.NaN
    cosTA = np.empty((len(tp),12))
    cosTA[:] = np.NaN

    for ii in range(0,len(tp)):
        tp_ii = tp[ii]
        lipid_ii = syst.select_atoms('resnum '+ str(tp_ii.resnum))
        lipid_fromDic = lipidTypeDic[lipid_ii.resnames[0]]

        tail1 = lipid_ii.select_atoms('name ' + lipid_fromDic.getTailBeads() + ' and (name **A or name GL1 AM1)')
        tail2 = lipid_ii.select_atoms('name ' + lipid_fromDic.getTailBeads() + ' and (name **B or name GL2 AM2)')
        pnorm = mt.memb_smooth.pnormals[ii]
        
        # Tail1
        bond_vectors = tail1.positions[range(0,len(tail1.positions)-1)][:] - tail1.positions[1:][:]
        for jj in range(0,len(bond_vectors)):
            # Correct bond vector if a bead crossed the periodic boundary. Note that the value 100 is arbitrary and might need to be adjusted for other systems. 
            if (np.linalg.norm(bond_vectors[jj][:]) > 100):
                bond_vectors[jj][:] = bond_vectors[jj][:] - np.multiply(bbox[1][:],np.multiply(np.abs(bond_vectors[jj][:])>100,np.sign(bond_vectors[jj][:])))
            # Check if the bond vector trully corresponds to two connected beads
            # (and not two beads that are not connected but further apart). 
            # If the beads are not connected (distance >7A), then write the bond length istead of cos_theta. 
            if (np.linalg.norm(bond_vectors[jj][:]) < 7):
                cos_theta = np.dot(bond_vectors[jj][:],pnorm)/np.linalg.norm(bond_vectors[jj][:])
                cosTA[ii][jj] = cos_theta
                ### P2t[ii][jj] = (0.5*(3*cos_theta**2 - 1))
            else:
                cosTA[ii][jj] = np.linalg.norm(bond_vectors[jj][:])
                
        # Tail2        
        bond_vectors = tail2.positions[range(0,len(tail2.positions)-1)][:] - tail2.positions[1:][:]
        for jj in range(0,len(bond_vectors)):
            if (np.linalg.norm(bond_vectors[jj][:]) > 100):
                bond_vectors[jj][:] = bond_vectors[jj][:] - np.multiply(bbox[1][:],np.multiply(np.abs(bond_vectors[jj][:])>100,np.sign(bond_vectors[jj][:])))
            if (np.linalg.norm(bond_vectors[jj][:]) < 7):
                cos_theta = np.dot(bond_vectors[jj][:],pnorm)/np.linalg.norm(bond_vectors[jj][:])
                cosTA[ii][6+jj] = -cos_theta
                ### P2t[ii][6+jj] = (0.5*(3*cos_theta**2 - 1)
            else:
                cosTA[ii][6+jj] = np.linalg.norm(bond_vectors[jj][:])
   
    ### return P2t
    return cosTA

# Function for computing the cosine of the dipole angles
def get_CosDipoleAngles(tp,mt,syst,lipidTypeDic):
    # Mn returns the cosine of the angle between the dipole and the local membrane normal.
    # dz is the z-component of the dipole vector. This should be multiplied with elementary charge and the applied electric field to get the potential elergy of dipoles. 
    Mn = np.empty((len(tp),1))
    Mn[:] = np.NaN
    dz = np.empty((len(tp),1))
    dz[:] = np.NaN

    for ii in range(0,len(tp)):
        tp_ii = tp[ii]
        lipid_ii = syst.select_atoms('resnum '+ str(tp_ii.resnum))
        lipid_fromDic = lipidTypeDic[lipid_ii.resnames[0]]

        headN = lipid_ii.select_atoms('name ' + lipid_fromDic.getHeadBeads() + ' and (name NC3 NH3)')
        headP = lipid_ii.select_atoms('name ' + lipid_fromDic.getHeadBeads() + ' and (name PO4)')
        pnorm = mt.memb_smooth.pnormals[ii]
        zaxis = [0,0,1]

        if len(headN):
            dipole_vector = headN.positions - headP.positions
            # Correct the dipole vector for cases where beads crossed the periodic boundaries
            if (np.linalg.norm(dipole_vector) > 100):
                dipole_vector = dipole_vector - np.multiply(bbox[1][:],np.multiply(np.abs(dipole_vector)>100,np.sign(dipole_vector)))
            cos_theta = np.dot(dipole_vector,pnorm)/np.linalg.norm(dipole_vector)
            cos_theta_zaxis = np.dot(dipole_vector,zaxis)/np.linalg.norm(dipole_vector)
            Mn[ii] = cos_theta
            dz[ii] = dipole_vector[0][2]
        
    return Mn, dz

# Function for extracting lipid charge
def get_charge(tp,mt,syst,lipidTypeDic):
    Q = np.zeros((len(tp),1))

    for ii in range(0,len(tp)):
        tp_ii = tp[ii]
        lipid_ii = syst.select_atoms('resnum '+ str(tp_ii.resnum))
        lipid_fromDic = lipidTypeDic[lipid_ii.resnames[0]]
        Q[ii] = lipid_fromDic.ldef[4]

    return Q

# ------------------------------------------------------------------------------
# MEMSURFER: compute membrane surfaces
# ------------------------------------------------------------------------------
if __name__ == '__main__':

    print ('using memsurfer from ({})'.format(memsurfer.__file__))

    # --------------------------------------------------------------------------
    # Path to input data.
    args = dict()
    args['traj'] = [os.path.join(ddir, xtcfile)]
    args['topo'] = [os.path.join(ddir, grofile)]

    # List of lipids (all lipids within APM and BPM)
    lipids = ['POPX','PIPX','DPPX','DAPC', 'DOPC', 'DPPC', 'OIPC', 'OUPC', 'PAPC', 'PEPC', 'PFPC', 'PIPC', 'POPC', 'PUPC', 'DAPE', 'DOPE', 'DUPE', 'OAPE', 'OIPE', 'OUPE', 'PAPE', 'PIPE', 'POPE', 'PQPE', 'PUPE', 'BNSM', 'DBSM', 'DPSM', 'DXSM', 'PBSM', 'PGSM', 'PNSM', 'POSM', 'XNSM', 'DAPS', 'DOPS', 'DPPS', 'DUPS', 'OUPS', 'PAPS', 'PIPS', 'POPS', 'PQPS', 'PUPS', 'PAPI', 'PIPI', 'POPI', 'PUPI', 'PAP1', 'PAP2', 'PAP3', 'POP1', 'POP2', 'POP3', 'PAPA', 'PIPA', 'POPA', 'PUPA', 'PADG', 'PIDG', 'PODG', 'PUDG', 'DBCE', 'DPCE', 'DXCE', 'PNCE', 'POCE', 'XNCE', 'APC', 'IPC', 'OPC', 'PPC', 'UPC', 'IPE', 'PPE', 'CHOL', 'DBG1', 'DPG1', 'DXG1', 'PNG1', 'POG1', 'XNG1', 'DBG3', 'DPG3', 'DXG3', 'PNG3', 'POG3', 'XNG3', 'DBGS', 'DPGS', 'PNGS', 'POGS']

    # --------------------------------------------------------------------------
    # Create a logger
    memsurfer.utils.create_logger(2,1,0,'','')

    # The prefix we will use to name the output files
    outprefix = args['traj'][0]
    outprefix = outprefix[:outprefix.rfind('.')]

    # --------------------------------------------------------------------------
    # Use MDAnalysis to read the data
    fargs = ['-f', args['traj'][0], '-s', args['topo'][0]]
    LOGGER.info('arguments = {}'.format(fargs))

    syst = mdreader.MDreader(fargs)
    syst.add_argument('-sframe', metavar='SELFRAME', type=int, dest='selframe', default=-1, help='int \tframe to save')
    syst.do_parse()

    LOGGER.info('Number of frames in sim {}'.format(len(syst)))
    LOGGER.info('System dimensions: {}'.format(syst.dimensions[0:3]))

    # Select all frames
    selFrame = range(0,len(syst))

    # Our domain is periodic (in xy)
    periodic = True

    # --------------------------------------------------------------------------
    # Use a lipid master list to identify headgrups
    lipidTypeList = lipid_masterlist()

    # Get list of all suported lipid types
    lipidTypeNames = np.array([l[0] for l in lipidTypeList])

    # Get all resnames in this simulation
    resnames = np.unique(syst.atoms.resnames)

    # Make a lipid type dictionary
    LOGGER.warning('Following resnames either not lipids or not supported {}'
                    .format(resnames[np.logical_not(np.in1d(resnames, lipidTypeNames))]))

    lipidTypeDic = {}
    for i in np.where(np.in1d(lipidTypeNames, resnames))[0]:
        lipidTypeDic[lipidTypeNames[i]] = LipidType(lipidTypeList[i], syst)


    # -------------------------------------------------------------------------
    # Define headgroups atoms

    # Headgroup beads of flip-floping lipids (**CE, **DG, CHOL)
    # The headgroup bead is the last element in a given lipidType described in the lipid_masterlist()    
    defFlipFlopHeadgroups = MDAnalysis.core.groups.AtomGroup([], syst)
    for i in list(lipidTypeDic.values()):
        defFlipFlopHeadgroups += i.getNoneLeafletSelection(syst)

    # Headgroup beads of not-flip-floping lipids (all except for the ones above)
    defHeadgroups = MDAnalysis.core.groups.AtomGroup([], syst)
    for i in list(lipidTypeDic.values()):
        defHeadgroups += i.getLeafletSelection(syst)

    # All headgroups in the system (not-flip-floping and flip-floping)
    allHeadgroups = defHeadgroups + defFlipFlopHeadgroups

    # get leaflets
    rcutoff, n = optimize_cutoff(syst, defHeadgroups)
    lfls = LeafletFinder(syst, defHeadgroups, cutoff=rcutoff, pbc=True)

    grps = lfls.groups()

    if len(grps) == 2:
        LOGGER.info('Found {} leaflet groups'.format(len(grps)))

    # check if they're even
    top_head = grps[0]
    bot_head = grps[1]

    rt = float(len(top_head))/len(bot_head)
    if rt > 1.3 or rt < 0.77:
        raise ValueError('Found uneven leaflets. top = {}, bot = {}'.format(len(top_head), len(bot_head)))

    LOGGER.info('Found leaflets of size {}, {}'.format(len(top_head), len(bot_head)))

    # --------------------------------------------------------------------------
    # Compute membranes for all frames
    for i in selFrame:
        
        print('\n')
        # Select frame 0 to len(syst)
        ts = syst.trajectory[i]
        LOGGER.info('Frame: %5d, Time: %8.3f ps' % (ts.frame, syst.trajectory.time))

        # And this is the bounding box
        bbox=np.zeros((2,3))
        bbox[1,:] = syst.dimensions[:3]


        # Get all lipids in top/bot leaflets (including flip-flop lipids - therefore has to be done for each frame)
        tp = top_head + defFlipFlopHeadgroups.select_atoms('around 12 global group topsel', topsel=top_head)
        bt = bot_head + defFlipFlopHeadgroups.select_atoms('around 12 global group botsel', botsel=bot_head)
        # Check if there is an atom selected in both leaflets and remove this atom from tp and bt 
        if len(tp.select_atoms('group bt', bt=bt)):
            bothLeaflets = tp.select_atoms('group bt', bt=bt)
            tp = tp - bothLeaflets
            bt = bt - bothLeaflets

        # Determine which out of all headgroup atoms are not selected in bt and tp. 
        # Loop over nonselected atoms and calculate its minimum distance from tp and bt. Place the atoms into the closest group (bt or tp)
        foundHeadgroups = tp + bt
        notfoundHeadgroupsIDs = np.setdiff1d(allHeadgroups.ids,foundHeadgroups.ids)
        for ii in range(0,len(notfoundHeadgroupsIDs)):
            a = allHeadgroups.select_atoms('bynum '+str(notfoundHeadgroupsIDs[ii]))
            min_dist_tp = np.min(MDAnalysis.analysis.distances.distance_array(a.positions, tp.positions))
            min_dist_bt = np.min(MDAnalysis.analysis.distances.distance_array(a.positions, bt.positions))
    
            if (min_dist_tp <= min_dist_bt):
                tp = tp + allHeadgroups.select_atoms('bynum '+str(notfoundHeadgroupsIDs[ii]))
            else:
                bt = bt + allHeadgroups.select_atoms('bynum '+str(notfoundHeadgroupsIDs[ii]))

        # Check if there are still some headgroup atoms found in both bt and tp groups
        if len(tp.select_atoms('group bt', bt=bt)):
            errmsg = 'Frame {}: {} common atoms between leaflets.'.format(syst.trajectory.frame, len(tp.select_atoms('group bt', bt=bt)))
            LOGGER.warning(errmsg)
            raise ValueError(errmsg)
        # Check if there are still some headgroup atoms not selected
        if (len(tp)+len(bt)) < len(allHeadgroups):
            errmsg = 'Frame {}: {} headgroup atoms not selected!'.format(syst.trajectory.frame, len(allHeadgroups)-(len(tp)+len(bt)))
            LOGGER.warning(errmsg)
            raise ValueError(errmsg)

        ## Select one bead per lipid (could do this better by getting the residues - and selecting e.g. first of each
        LOGGER.info('We have {} lipids in upper and {} in lower leaflets'.format(len(tp), len(bt)))

        # ----------------------------------------------------------------------
        # This is where we use MemSurfer
        # ----------------------------------------------------------------------
        mt = memsurfer.Membrane.compute(tp.positions, tp.resnames, bbox, periodic)
        mb = memsurfer.Membrane.compute(bt.positions, bt.resnames, bbox, periodic)

        # compute total density
        ####sigmas = [5,10,15]
        ####sigmas = [5,7,9]
        ####get_nlipdis = True
        ####memsurfer.Membrane.compute_densities([mt, mb], [2], sigmas, get_nlipdis, 'all')
        
        # compute density of each type of lipid
        ####for l in lipids:
        ####    memsurfer.Membrane.compute_densities([mt, mb], [2], sigmas, get_nlipdis, l)

        # compute thickness
        memsurfer.Membrane.compute_thickness(mt, mb, 'smooth')

	# compute cosines of tail bond angles
        cosTAt = get_CosTailAngles(tp,mt,syst,lipidTypeDic)
        cosTAb = get_CosTailAngles(bt,mb,syst,lipidTypeDic)

        # compute cosines of headgroup dipole angles
        Mnt, dzt = get_CosDipoleAngles(tp,mt,syst,lipidTypeDic)
        Mnb, dzb = get_CosDipoleAngles(bt,mb,syst,lipidTypeDic)

        # get lipid charge
        Qt = get_charge(tp,mt,syst,lipidTypeDic)
        Qb = get_charge(bt,mb,syst,lipidTypeDic)

        # Write arrays into txt files
        fname = PM+'/'+mem+'/'+folder+'/'
        np.savetxt(fname+'box_'+str(i)+'.txt',mt.bbox,fmt='%.4f',delimiter=' ',newline='\n')  
        # top membrane        
        np.savetxt(fname+'pnormt_'+str(i)+'.txt',mt.memb_smooth.pnormals,fmt='%.4f',delimiter=' ',newline='\n')  
        np.savetxt(fname+'pareat_'+str(i)+'.txt',mt.memb_smooth.pareas,fmt='%.4f',delimiter=' ',newline='\n')        
        np.savetxt(fname+'mcurvt_'+str(i)+'.txt',mt.memb_smooth.mean_curv,fmt='%.4f',delimiter=' ',newline='\n')
        np.savetxt(fname+'gcurvt_'+str(i)+'.txt',mt.memb_smooth.gaus_curv,fmt='%.4f',delimiter=' ',newline='\n')
        np.savetxt(fname+'pointt_'+str(i)+'.txt',mt.points,fmt='%.4f',delimiter=' ',newline='\n')
        np.savetxt(fname+'vertit_'+str(i)+'.txt',mt.memb_smooth.vertices,fmt='%.4f',delimiter=' ',newline='\n')
        np.savetxt(fname+'thickt_'+str(i)+'.txt',mt.properties['thickness'],fmt='%.4f',delimiter=' ',newline='\n')
        np.savetxt(fname+'ordert_'+str(i)+'.txt',cosTAt,fmt='%.4f',delimiter=' ',newline='\n')        
        np.savetxt(fname+'labelt_'+str(i)+'.txt',mt.labels,fmt='%s',delimiter=' ',newline='\n')
        np.savetxt(fname+'diplnt_'+str(i)+'.txt',Mnt,fmt='%.4f',delimiter=' ',newline='\n') 
        np.savetxt(fname+'dipvzt_'+str(i)+'.txt',dzt,fmt='%.4f',delimiter=' ',newline='\n') 
        np.savetxt(fname+'chargt_'+str(i)+'.txt',Qt,fmt='%.4f',delimiter=' ',newline='\n') 
        np.savetxt(fname+'resnumt_'+str(i)+'.txt',tp.resnums,fmt='%.0f',delimiter=' ',newline='\n')
	
        # bottom membrane        
        np.savetxt(fname+'pnormb_'+str(i)+'.txt',mb.memb_smooth.pnormals,fmt='%.4f',delimiter=' ',newline='\n')  
        np.savetxt(fname+'pareab_'+str(i)+'.txt',mb.memb_smooth.pareas,fmt='%.4f',delimiter=' ',newline='\n')                        
        np.savetxt(fname+'mcurvb_'+str(i)+'.txt',mb.memb_smooth.mean_curv,fmt='%.4f',delimiter=' ',newline='\n')
        np.savetxt(fname+'gcurvb_'+str(i)+'.txt',mb.memb_smooth.gaus_curv,fmt='%.4f',delimiter=' ',newline='\n')
        np.savetxt(fname+'pointb_'+str(i)+'.txt',mb.points,fmt='%.4f',delimiter=' ',newline='\n')
        np.savetxt(fname+'vertib_'+str(i)+'.txt',mb.memb_smooth.vertices,fmt='%.4f',delimiter=' ',newline='\n')
        np.savetxt(fname+'thickb_'+str(i)+'.txt',mb.properties['thickness'],fmt='%.4f',delimiter=' ',newline='\n')
        np.savetxt(fname+'orderb_'+str(i)+'.txt',cosTAb,fmt='%.4f',delimiter=' ',newline='\n')
        np.savetxt(fname+'labelb_'+str(i)+'.txt',mb.labels,fmt='%s',delimiter=' ',newline='\n')
        np.savetxt(fname+'diplnb_'+str(i)+'.txt',Mnb,fmt='%.4f',delimiter=' ',newline='\n') 
        np.savetxt(fname+'dipvzb_'+str(i)+'.txt',dzb,fmt='%.4f',delimiter=' ',newline='\n') 
        np.savetxt(fname+'chargb_'+str(i)+'.txt',Qb,fmt='%.4f',delimiter=' ',newline='\n') 	
        np.savetxt(fname+'resnumb_'+str(i)+'.txt',bt.resnums,fmt='%.0f',delimiter=' ',newline='\n')

        ####for s in sigmas:  
            ####name0 = 'density_type2_all_k{0:.1f}'.format(s)   
            ####np.savetxt(fname+'densS'+str(s)+'t_'+str(i)+'.txt',mt.properties[name0],fmt='%.4f',delimiter=' ',newline='\n')
            ####np.savetxt(fname+'densS'+str(s)+'b_'+str(i)+'.txt',mb.properties[name0],fmt='%.4f',delimiter=' ',newline='\n')

            ####denslt = []
            ####denslb = []
            ####for l in lipids:   
            ####    name = 'density_type2_'+l+'_k{0:.1f}'.format(s)
            ####    tmpt = mt.properties[name]
            ####    tmpb = mb.properties[name]
            ####    denslt.append(tmpt)
            ####    denslb.append(tmpb)                
            ####np.savetxt(fname+'densLipt_'+str(i)+'_s'+str(s)+'.txt',denslt,fmt='%.4f',delimiter=' ',newline='\n')
            ####np.savetxt(fname+'densLipb_'+str(i)+'_s'+str(s)+'.txt',denslb,fmt='%.4f',delimiter=' ',newline='\n')
                
        # Write vtp files
        #if (i == selFrame[0]):
        if False:
            mt.write_all(outprefix+'_f{}-top'.format(ts.frame),
                    {'frame': ts.frame, 'time': syst.trajectory.time})

            mb.write_all(outprefix+'_f{}-bot'.format(ts.frame),
                    {'frame': ts.frame, 'time': syst.trajectory.time})

        # --------------------------------------------------------------------------
        # --------------------------------------------------------------------------


