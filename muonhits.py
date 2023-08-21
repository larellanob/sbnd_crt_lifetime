#!/usr/bin/python3

# imports
import uproot
import matplotlib.pyplot as plt
import numpy as np


nevents = None # number of events to read (set to None for all events)

# open root file
input_root_file = "muon_hitdumper_NS.root"
try:
    uproot_file = uproot.open(input_root_file)
except OSError:
    print("Can't find input root file. Terminating")
    print(input_root_file)
    exit()

# grab tree and relevant branches
tree = uproot_file['hitdumper/hitdumpertree']
muonhit_branches = [
    "run",
    "subrun",
    "event",
    "t0",
    "nmuontrks",   # muon tracks
    "muontrk_tpc", # muon track tpc
    "nmhits",
    "mhit_trk",
    "mhit_tpc",
    "mhit_plane",
    "mhit_wire",
    "mhit_channel",
    "mhit_peakT",
    "mhit_charge"
]

test_only = False
output_dir = "img/"
if test_only:
    output_dir += "test/"
    nevents = 10
else:
    output_dir += "full/"

# arrays is to read everything in a single step, don't need it, I
# would like to iterate through events

#muonhits = tree.arrays(muonhit_branches, library="np", entry_stop=nevents)
#mhits = tree.arrays(muonhit_branches, entry_stop=nevents)
#print("branches",mhits)


mhits_it = tree.iterate(muonhit_branches,step_size=1,entry_stop=nevents)

# loop over events
e_counter=0
for e in mhits_it:
    #print(e["nmhits"],e["mhit_trk"],len(e["mhit_trk"][0]))
    fig, axs = plt.subplots(3,2,sharex=True,sharey=True)
    plt.subplots_adjust(left=0.11, bottom=0.1, right=0.96, top=0.88, wspace=0.05, hspace=0.0)
    run, subrun, event = e["run"][0],e["subrun"][0],e["event"][0]
    nmuontracks = e["nmuontrks"][0]
    fig.suptitle(f"mhits per plane, Run: {run}, Subrun: {subrun}, Event: {event} - Muon tracks: {nmuontracks}")
    nhits = e["nmhits"][0]
    # loop over hits

    '''
    # in this loop, goes hit by hit and depending on which plane, adds
    # point to specific subplot
    for i in range(0,nhits):
        if i % 50 != 0:
            continue
        plane  = e["mhit_plane"][0][i]
        wire   = e["mhit_wire"][0][i]
        time   = e["mhit_peakT"][0][i]
        charge = e["mhit_charge"][0][i]
        #print(plane,wire,time,charge)
        axs[plane].scatter(wire,time,c=charge,marker=',',s=1)
        #axs[plane].scatter(wire,time,c=charge)
    plt.show()

    '''

    # in this loop we do a first pass to create different lists
    # according to hit plane, so we can more easily assign color using
    # charge
    # first 3 lists: tpc 0, last 3: tpc 1
    wires = [[],[],[],[],[],[]]
    times = [[],[],[],[],[],[]]
    charges = [[],[],[],[],[],[]]
    for i in range(0,nhits):
        if test_only and i % 5 != 0:
            continue
        plane  = e["mhit_plane"][0][i]
        #print(event,plane,i)
        tpc = e["mhit_tpc"][0][i]
        wires[3*tpc+plane].append(e["mhit_wire"][0][i])
        times[3*tpc+plane].append(e["mhit_peakT"][0][i])
        charges[3*tpc+plane].append(e["mhit_charge"][0][i])
    for tpc in [0,1]:
        for plane in [0, 1, 2]:
            axs[plane][1].set_ylabel(f"Plane {plane}")
            axs[plane][1].yaxis.set_label_position("right")
                #axs[plane][tpc].set_title(f"Plane {plane}",x = 1.04,y = 0.25,rotation = -90)
            axs[plane][tpc].scatter(wires[3*tpc+plane],times[3*tpc+plane],c=charges[3*tpc+plane],marker=',',s=1)
            axs[plane][tpc].set_xlim(-20,1900)
            axs[plane][tpc].set_ylim(-500,4500)
        plt.setp(axs[0][tpc].get_xticklabels(), visible=False)
        plt.setp(axs[1][tpc].get_xticklabels(), visible=False)
    #axs[2].set_xlabel("Wire")
    #axs[1].set_ylabel("Peak time")
    axs[2][0].set_xlabel("Wire")
    axs[2][1].set_xlabel("Wire")
    tpc_tracks = [0,0]
    for muontrk in range(0,e["nmuontrks"][0]):
        tpc_tracks[e["muontrk_tpc"][0][muontrk]] += 1
    axs[0][0].set_title(f"TPC 0 ({tpc_tracks[0]} muontracks)")
    axs[0][1].set_title(f"TPC 1 ({tpc_tracks[1]} muontracks)")
    axs[1][0].set_ylabel("Time")
    #axs[0][1].text(4,2.0,"random text",rotation=-90)
    
    
    
    plt.savefig(f"{output_dir}/mhits_evd_{e_counter}_{run}_{subrun}_{event}.png")
    plt.savefig(f"{output_dir}/mhits_evd_{e_counter}_{run}_{subrun}_{event}.pdf")
    e_counter+=1
    plt.close()
