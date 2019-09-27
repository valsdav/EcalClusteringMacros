from __future__ import print_function
from collections import defaultdict
from operator import itemgetter, attrgetter

def get_hits_maps(event, debug=False):
    ncalo = event.caloParticle_pt.size()
    calo_simE = event.caloParticle_simEnergy
    pfCluster_energy = event.pfCluster_energy
    # map (ieta,iphi,icalo, simhit):(iclu, clhit)
    xtal_calo = defaultdict(list)
    # map (ieta, iphi, iclu, clhit):(icalo, simhit)
    xtal_cluster = defaultdict(list)
    # map (ieta, iphi):(iclu, noisehit)
    xtal_cluster_noise = defaultdict(list)

    for icalo in range(ncalo):
        if debug: 
            print("--- icalo: ", icalo)
            print ("ieta iphi simhit [ pfcluster index , pfcluster hit]")
        for i, (ieta, iphi, simhit, clhit) in enumerate(zip(event.simHit_ieta[icalo],  event.simHit_iphi[icalo],
                        event.simHit_energy[icalo], event.pfClusterHit_energy[icalo])):
            if clhit.size() > 0:
                if debug: print(ieta, iphi, "{:.5f}".format(simhit), [(hit.first, '{:.5f}'.format(hit.second)) for hit in clhit]) 
                for chit in clhit:
                    xtal_cluster[(ieta, iphi, chit.first, chit.second)].append((icalo, simhit))
                    xtal_calo[(ieta, iphi, icalo, simhit)].append((chit.first, chit.second))
        #
    # Check the noise hits (not overlapping with caloparticle sihimt)
    for nclus , (energys, ietas,iphis) in enumerate(zip(event.pfClusterHit_noCaloPart_energy,   
                    event.pfClusterHit_noCaloPart_ieta, event.pfClusterHit_noCaloPart_iphi)):
        #print nclus , [(ieta, iphi, en) for ieta,iphi, en in zip(energys, ietas, iphis)]
        for en, ieta, iphi in zip(energys, ietas, iphis):
            xtal_cluster_noise[(ieta, iphi)].append((nclus, en))

    return xtal_cluster, xtal_calo, xtal_cluster_noise

########################################################################
# Association strategies 

def sim_fraction_stragegy(xtal_cluster, xtal_calo, event, min_fraction= 0):
    '''
     - Associate to each pfCluster the calo with the greatest fraction of 
       its simEnergy. (save the ordered list of them by fraction)
    - Save a list of pfcluster associated with the method before to 
       a caloparticle (there can be more than one)
    '''
    calo_simE = event.caloParticle_simEnergy
    pfCluster_energy = event.pfCluster_energy

    ####################################################
   
    #######################################################
    cluster_calo_fraction = defaultdict(dict)
    for (ieta,iphi,clid, clhit), (caloinfo) in xtal_cluster.items():
        for icalo, simhit in caloinfo:
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



def sim_rechit_fractions_strategy(xtal_cluster, xtal_calo, event):
    calo_simE = event.caloParticle_simEnergy
    pfCluster_energy = event.pfCluster_energy

    cluster_calo_fraction = defaultdict(dict)
    for (ieta,iphi,clid, clhit), (caloinfo) in xtal_cluster.items():
        for icalo, simhit in caloinfo:
            if icalo not in cluster_calo_fraction[clid]:
                cluster_calo_fraction[clid][icalo] = 0.
            cluster_calo_fraction[clid][icalo] += abs( (simhit / calo_simE[icalo])-(clhit / pfCluster_energy[clid]) )
               
    calo_cluster_assoc = defaultdict(list)
    cluster_calo_assoc = {}

    for clid, calofrac in cluster_calo_fraction.items():
        # order caloparticle. The closer to 0 the better
        caloids =  list(sorted(calofrac.items(), key=itemgetter(1)))
        # Associate to che cluster the list of caloparticles ordered by fraction 
        cluster_calo_assoc[clid] = caloids
        # save for the calocluster in the caloparticle if it is the one with more fraction
        # This is necessary in case one caloparticle is linked with more than one cluster
        calo_cluster_assoc[caloids[0][0]].append((clid, caloids[0][1] ))
        
    # Now sort the clusters associated to a caloparticle with the fraction (closer to 0 the better)
    sorted_calo_cluster_assoc = {}
    for caloid, clinfo in calo_cluster_assoc.items():
        sorted_calo_cluster_assoc[caloid] = list(map(itemgetter(0), sorted(clinfo, key=itemgetter(1))))

    # Return cluster:caloparticle and  caloparticle_cluster maps
    return cluster_calo_assoc, sorted_calo_cluster_assoc



def sim_rechit_global_fraction_strategy(xtal_cluster, xtal_calo, event):
    calo_simE = event.caloParticle_simEnergy
    pfCluster_energy = event.pfCluster_energy

    cluster_calo_fraction = defaultdict(dict)
    for (ieta,iphi,clid, clhit), (caloinfo) in xtal_cluster.items():
        for icalo, simhit in caloinfo:
            if icalo not in cluster_calo_fraction[clid]:
                cluster_calo_fraction[clid][icalo] = 0.
            # no absolute value! We will calculate the global one
            cluster_calo_fraction[clid][icalo] +=  (simhit / calo_simE[icalo])-(clhit / pfCluster_energy[clid])
            
    
    calo_cluster_assoc = defaultdict(list)
    cluster_calo_assoc = {}

    for clid, calofrac in cluster_calo_fraction.items():
        # order caloparticle. The closer to 0 the better
        # Calculate not the absolute value of the fraction
        caloids =  list(sorted(map (lambda (k,v): (k, abs(v)), calofrac.items()), key=itemgetter(1)))
        # Associate to che cluster the list of caloparticles ordered by fraction 
        cluster_calo_assoc[clid] = caloids
        # save for the calocluster in the caloparticle if it is the one with more fraction
        # This is necessary in case one caloparticle is linked with more than one cluster
        calo_cluster_assoc[caloids[0][0]].append((clid, caloids[0][1] ))
        
    # Now sort the clusters associated to a caloparticle with the fraction (closer to 0 the better)
    sorted_calo_cluster_assoc = {}
    for caloid, clinfo in calo_cluster_assoc.items():
        sorted_calo_cluster_assoc[caloid] = list(map(itemgetter(0), sorted(clinfo, key=itemgetter(1))))

    # Return cluster:caloparticle and  caloparticle_cluster maps
    return cluster_calo_assoc, sorted_calo_cluster_assoc

####################################################################################


strategies = {
    "sim_fraction": sim_fraction_stragegy,
    "sim_rechit_fractions": sim_rechit_fractions_strategy,
    "sim_rechit_global_fraction": sim_rechit_global_fraction_strategy
}


def get_all_associations(event, debug=False):
    xtal_cluster, xtal_calo, xtal_cluster_noise = get_hits_maps(event, debug)
    assoc = {}
    for strategy, f in strategies.items():
        assoc[strategy] = f(xtal_cluster, xtal_calo, event)
    return assoc, (xtal_cluster, xtal_calo, xtal_cluster_noise)