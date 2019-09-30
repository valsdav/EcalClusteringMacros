from __future__ import print_function
from collections import defaultdict
from operator import itemgetter, attrgetter

def get_hits_maps(event, cluster_type="pfCluster", debug=False):
    ncalo = event.caloParticle_pt.size()
    calo_simE = event.caloParticle_simEnergy
    cluster_energy = getattr(event, cluster_type +"_energy")
    cluster_hits_energy = getattr(event, cluster_type +"Hit_energy")
    cluster_hits_nocalo_energy = getattr(event, cluster_type +"Hit_noCaloPart_energy")
    cluster_hits_nocalo_ieta = getattr(event, cluster_type +"Hit_noCaloPart_ieta")
    cluster_hits_nocalo_iphi = getattr(event, cluster_type +"Hit_noCaloPart_iphi")
    # map (ieta,iphi,icalo, simhit):(iclu, clhit)
    xtal_calo = defaultdict(list)
    # map (ieta, iphi, iclu, clhit):(icalo, simhit)
    xtal_cluster = defaultdict(dict)
    # map (ieta, iphi):(iclu, noisehit)
    xtal_cluster_noise = defaultdict(list)

    for icalo in range(ncalo):
        if debug: 
            print("--- icalo: ", icalo)
            print ("ieta iphi simhit [ cluster index , cluster hit]")
        for i, (ieta, iphi, simhit, clhit) in enumerate(zip(event.simHit_ieta[icalo],  event.simHit_iphi[icalo],
                        event.simHit_energy[icalo], cluster_hits_energy[icalo])):
            if clhit.size() > 0:
                if debug: print(ieta, iphi, "{:.5f}".format(simhit), [(hit.first, '{:.5f}'.format(hit.second)) for hit in clhit]) 
                for chit in clhit:
                    xtal_cluster[(ieta, iphi, chit.first, chit.second)][icalo] = simhit
                    xtal_calo[(ieta, iphi, icalo, simhit)].append((chit.first, chit.second))
        #
    # Check the noise hits (not overlapping with caloparticle sihimt)
    for nclus , (energys, ietas,iphis) in enumerate(zip(cluster_hits_nocalo_energy,   
                    cluster_hits_nocalo_ieta, cluster_hits_nocalo_iphi)):
        #print nclus , [(ieta, iphi, en) for ieta,iphi, en in zip(energys, ietas, iphis)]
        for en, ieta, iphi in zip(energys, ietas, iphis):
            xtal_cluster_noise[(ieta, iphi)].append((nclus, en))

    return xtal_cluster, xtal_calo, xtal_cluster_noise

########################################################################
# Association strategies 

def sim_fraction_stragegy(xtal_cluster, xtal_calo, event, cluster_type="pfCluster", min_fraction= 0):
    '''
     - Associate to each pfCluster the calo with the greatest fraction of 
       its simEnergy. (save the ordered list of them by fraction)
    - Save a list of pfcluster associated with the method before to 
       a caloparticle (there can be more than one)
    '''
    calo_simE = event.caloParticle_simEnergy
    cluster_energy = getattr(event, cluster_type+"_energy")
    ####################################################
   
    #######################################################
    cluster_calo_fraction = defaultdict(dict)
    for (ieta,iphi,clid, clhit), (caloinfo) in xtal_cluster.items():
        for icalo, simhit in caloinfo.items():
            if icalo not in cluster_calo_fraction[clid]:
                cluster_calo_fraction[clid][icalo] = 0.
            cluster_calo_fraction[clid][icalo] += simhit / calo_simE[icalo]
    
    # Filter out calo with less than 5% fraction 
    clean_cluster_calo_fraction = defaultdict(dict)
    if min_fraction == 0: 
        clean_cluster_calo_fraction = cluster_calo_fraction
    else:
        for clid, calos in cluster_calo_fraction.items():
            for icalo, frac in calos.items():
                if frac > min_fraction:
                    clean_cluster_calo_fraction[clid][icalo] = frac

    
    calo_cluster_assoc = defaultdict(list)
    cluster_calo_assoc = {}

    for clid, calofrac in clean_cluster_calo_fraction.items():
        # order caloparticle by fraction
        caloids =  list(sorted(calofrac.items(), key=itemgetter(1), reverse=True))
        # Associate to che cluster the list of caloparticles ordered by fraction 
        cluster_calo_assoc[clid] = caloids
        # save for the calocluster in the caloparticle if it is the one with more fraction
        # This is necessary in case one caloparticle is linked with more than one cluster
        calo_cluster_assoc[caloids[0][0]].append((clid, caloids[0][1] ))
        
    # Now sort the clusters associated to a caloparticle with the fraction 
    sorted_calo_cluster_assoc = {}
    for caloid, clinfo in calo_cluster_assoc.items():
        sorted_calo_cluster_assoc[caloid] = list(map(itemgetter(0), sorted(clinfo, key=itemgetter(1), reverse=True)))

    # Return cluster:caloparticle and  caloparticle_cluster maps
    return cluster_calo_assoc, sorted_calo_cluster_assoc




def sim_rechit_fractions_strategy(xtal_cluster, xtal_calo, event, cluster_type="pfCluster"):
    calo_simE = event.caloParticle_simEnergy
    cluster_energy = getattr(event, cluster_type+"_energy")

    ncalos_hits = defaultdict(dict)
    cluster_calo_fraction = defaultdict(dict)
    for (ieta,iphi,clid, clhit), (caloinfo) in xtal_cluster.items():
        for icalo in range(len(calo_simE)):
            if icalo not in cluster_calo_fraction[clid]:
                cluster_calo_fraction[clid][icalo] = 0.
                ncalos_hits[clid][icalo] = 0
            if icalo in caloinfo:
                cluster_calo_fraction[clid][icalo] += abs( (caloinfo[icalo] / calo_simE[icalo])-(clhit / cluster_energy[clid]) )
                ncalos_hits[clid][icalo] +=1
            else:
                cluster_calo_fraction[clid][icalo] += abs( (clhit / cluster_energy[clid]) )

    calo_cluster_assoc = defaultdict(list)
    cluster_calo_assoc = {}

    for clid, calofrac in cluster_calo_fraction.items():
        # order caloparticle. 
        # calculate new fraction as 1 - fraction
        caloids =  list(sorted(
                            map(lambda (k,f): (k, 1-f), 
                                filter(lambda (k,v): ncalos_hits[clid][k]>0, calofrac.items())
                            ),
                        key=itemgetter(1), reverse=True))
        # Associate to che cluster the list of caloparticles ordered by fraction 
        cluster_calo_assoc[clid] = caloids
        # save for the calocluster in the caloparticle if it is the one with more fraction
        # This is necessary in case one caloparticle is linked with more than one cluster
        calo_cluster_assoc[caloids[0][0]].append((clid, caloids[0][1] ))
        
    # Now sort the clusters associated to a caloparticle with the fraction (greater the better because we have done 1 -frac)
    sorted_calo_cluster_assoc = {}
    for caloid, clinfo in calo_cluster_assoc.items():
        sorted_calo_cluster_assoc[caloid] = list(map(itemgetter(0), sorted(clinfo, key=itemgetter(1), reverse=True)))

    # Return cluster:caloparticle and  caloparticle_cluster maps
    return cluster_calo_assoc, sorted_calo_cluster_assoc



def sim_rechit_global_fraction_strategy(xtal_cluster, xtal_calo, event,cluster_type="pfCluster"):
    calo_simE = event.caloParticle_simEnergy
    cluster_energy = getattr(event, cluster_type+"_energy")

    ncalos_hits = defaultdict(dict)
    cluster_calo_fraction = defaultdict(dict)
    for (ieta,iphi,clid, clhit), (caloinfo) in xtal_cluster.items():
        for icalo, simhit in caloinfo.items():
            if icalo not in cluster_calo_fraction[clid]:
                cluster_calo_fraction[clid][icalo] = 0.
            cluster_calo_fraction[clid][icalo] +=  (simhit/ calo_simE[icalo])-(clhit / cluster_energy[clid]) 
            

    calo_cluster_assoc = defaultdict(list)
    cluster_calo_assoc = {}

    for clid, calofrac in cluster_calo_fraction.items():
        # order caloparticle. 
        # calculate new fraction as 1 - fraction
        caloids =  list(sorted(
                            map(lambda (k,f): (k, 1-abs(f)), calofrac.items()),
                        key=itemgetter(1), reverse=True))
        # Associate to che cluster the list of caloparticles ordered by fraction 
        cluster_calo_assoc[clid] = caloids
        # save for the calocluster in the caloparticle if it is the one with more fraction
        # This is necessary in case one caloparticle is linked with more than one cluster
        calo_cluster_assoc[caloids[0][0]].append((clid, caloids[0][1] ))
        
    # Now sort the clusters associated to a caloparticle with the fraction (greater the better because we have done 1 -frac)
    sorted_calo_cluster_assoc = {}
    for caloid, clinfo in calo_cluster_assoc.items():
        sorted_calo_cluster_assoc[caloid] = list(map(itemgetter(0), sorted(clinfo, key=itemgetter(1), reverse=True)))

    # Return cluster:caloparticle and  caloparticle_cluster maps
    return cluster_calo_assoc, sorted_calo_cluster_assoc


def sim_rechit_diff_strategy(xtal_cluster, xtal_calo, event,cluster_type="pfCluster"):
    calo_simE = event.caloParticle_simEnergy
    cluster_energy = getattr(event, cluster_type+"_energy")

    ncalos_hits = defaultdict(dict)
    cluster_calo_fraction = defaultdict(dict)
    for (ieta,iphi,clid, clhit), (caloinfo) in xtal_cluster.items():
        for icalo, simhit in caloinfo.items():
            if icalo not in cluster_calo_fraction[clid]:
                cluster_calo_fraction[clid][icalo] = 0.
                ncalos_hits[clid][icalo] = 0
            cluster_calo_fraction[clid][icalo] +=  abs(simhit - clhit)
            ncalos_hits[clid][icalo] += 1
            
    mean_cluster_calo_fraction = defaultdict(dict)
    for clid, calos in cluster_calo_fraction.items():
        for icalo, frac in calos.items():
            mean_cluster_calo_fraction[clid][icalo] = frac / ncalos_hits[clid][icalo]
        
    calo_cluster_assoc = defaultdict(list)
    cluster_calo_assoc = {}

    for clid, calofrac in mean_cluster_calo_fraction.items():
        # order caloparticle. 
        # The score is 1 - abs(fraction). fraction = simhit / cluster energy
        caloids =  list(sorted(
                            map (lambda (k,v): (k, 1- v), calofrac.items()),
                         key=itemgetter(1), reverse=True))
        # Associate to che cluster the list of caloparticles ordered by fraction 
        cluster_calo_assoc[clid] = caloids
        # save for the calocluster in the caloparticle if it is the one with more fraction
        # This is necessary in case one caloparticle is linked with more than one cluster
        calo_cluster_assoc[caloids[0][0]].append((clid, caloids[0][1] ))
        
    # Now sort the clusters associated to a caloparticle with the fraction 
    sorted_calo_cluster_assoc = {}
    for caloid, clinfo in calo_cluster_assoc.items():
        sorted_calo_cluster_assoc[caloid] = list(map(itemgetter(0), sorted(clinfo, key=itemgetter(1), reverse=True)))

    # Return cluster:caloparticle and  caloparticle_cluster maps
    return cluster_calo_assoc, sorted_calo_cluster_assoc

def sim_rechit_tot_diff_strategy(xtal_cluster, xtal_calo, event,cluster_type="pfCluster"):
    calo_simE = event.caloParticle_simEnergy
    cluster_energy = getattr(event, cluster_type+"_energy")
    
    cluster_calo_fraction = defaultdict(dict)
    for (ieta,iphi,clid, clhit), (caloinfo) in xtal_cluster.items():
        for icalo, simhit in caloinfo.items():
            if icalo not in cluster_calo_fraction[clid]:
                cluster_calo_fraction[clid][icalo] = 0.
            cluster_calo_fraction[clid][icalo] +=  simhit - clhit
            
    calo_cluster_assoc = defaultdict(list)
    cluster_calo_assoc = {}

    for clid, calofrac in cluster_calo_fraction.items():
        # order caloparticle. 
        # The score is 1 - abs(fraction). fraction = simhit / cluster energy
        caloids =  list(sorted(
                            map (lambda (k,v): (k, 1- abs(v)), calofrac.items()),
                         key=itemgetter(1), reverse=True))
        # Associate to che cluster the list of caloparticles ordered by fraction 
        cluster_calo_assoc[clid] = caloids
        # save for the calocluster in the caloparticle if it is the one with more fraction
        # This is necessary in case one caloparticle is linked with more than one cluster
        calo_cluster_assoc[caloids[0][0]].append((clid, caloids[0][1] ))
        
    # Now sort the clusters associated to a caloparticle with the fraction 
    sorted_calo_cluster_assoc = {}
    for caloid, clinfo in calo_cluster_assoc.items():
        sorted_calo_cluster_assoc[caloid] = list(map(itemgetter(0), sorted(clinfo, key=itemgetter(1), reverse=True)))

    # Return cluster:caloparticle and  caloparticle_cluster maps
    return cluster_calo_assoc, sorted_calo_cluster_assoc


####################################################################################


strategies = {
    "sim_fraction": sim_fraction_stragegy,
    "sim_rechit_fractions": sim_rechit_fractions_strategy,
    "sim_rechit_global_fraction": sim_rechit_global_fraction_strategy,
    "sim_rechit_diff": sim_rechit_diff_strategy,
    #"sim_rechit_tot_diff": sim_rechit_tot_diff_strategy
}


def get_all_associations(event, cluster_type="pfCluster", debug=False):
    xtal_cluster, xtal_calo, xtal_cluster_noise = get_hits_maps(event, cluster_type=cluster_type, debug=debug)
    assoc = {}
    for strategy, f in strategies.items():
        assoc[strategy] = f(xtal_cluster, xtal_calo, event, cluster_type=cluster_type)
    return assoc, (xtal_cluster, xtal_calo, xtal_cluster_noise)