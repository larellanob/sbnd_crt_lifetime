#!/usr/bin/python3

# imports
import uproot # for reading root files
import matplotlib.pyplot as plt # general figure plotting
import matplotlib.patches as pltpat # for drawing rectangle in collection plane less wires
import numpy as np # sign function and getting argument of minimum in list
import math # atan2 and radian-degree conversion
import os # make subdirs for run/subrun


test_only = True
draw_auxiliary_plane2 = True
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


if test_only:
    nevents = 5

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
    # full 2 TPC + 3 planes Event display
    fig, axs = plt.subplots(3,2,sharex=True,sharey=True)
    plt.subplots_adjust(left=0.11, bottom=0.1, right=0.96, top=0.88, wspace=0.0, hspace=0.0)

    # "declare" zoomed figure into "scaled" collection plane
    fig2, axs2 = plt.subplots(2,1,sharex = True)
    plt.subplots_adjust(left=0.0, bottom=0.1, right=1.0, top=0.88, wspace=0.05, hspace=0.0)
    # switch back to main figure for now
    plt.figure(fig)
        
    run, subrun, event = e["run"][0],e["subrun"][0],e["event"][0]
    nmuontracks = e["nmuontrks"][0]
    fig.suptitle(f"mhits per plane, Run: {run}, Subrun: {subrun}, Event: {event} - Muon tracks: {nmuontracks}")
    fig2.suptitle(f"Collection plane, Run: {run}, Subrun: {subrun}, Event: {event} - Muon tracks: {nmuontracks}")
    nhits = e["nmhits"][0]
    # loop over hits

    # in this loop we do a first pass to create different lists
    # according to hit plane, so we can more easily assign color using
    # charge
    # first 3 lists: tpc 0, last 3: tpc 1
    wires = [[],[],[],[],[],[]]
    times = [[],[],[],[],[],[]]
    charges = [[],[],[],[],[],[]]
    hit_track = [[],[],[],[],[],[]]
    for i in range(0,nhits):
        if test_only and i % 5 != 0:
            continue
        plane  = e["mhit_plane"][0][i]
        #print(event,plane,i)
        tpc = e["mhit_tpc"][0][i]
        wires[3*tpc+plane].append(e["mhit_wire"][0][i])
        times[3*tpc+plane].append(e["mhit_peakT"][0][i])
        charges[3*tpc+plane].append(e["mhit_charge"][0][i])
        hit_track[3*tpc+plane].append(e["mhit_trk"][0][i])
    for tpc in [0,1]:
        for plane in [0, 1, 2]:
            axs[plane][1].set_ylabel(f"Plane {plane}")
            axs[plane][1].yaxis.set_label_position("right")
                #axs[plane][tpc].set_title(f"Plane {plane}",x = 1.04,y = 0.25,rotation = -90)
            axs[plane][tpc].scatter(wires[3*tpc+plane],times[3*tpc+plane],c=charges[3*tpc+plane],marker=',',s=1)
            axs[plane][tpc].set_xlim(0,1983)
            axs[plane][tpc].set_ylim(0,3400)
            if plane == 2:
                rectangle = pltpat.Rectangle((1663,0),1983-1663,3400,fill=False,hatch='XX')
                axs[plane][tpc].add_patch(rectangle)

                # switch to plane2 figure
                '''

                Note: TPCs are switched vertically. Going from top to
                bottom, matplotlib.subplots are stored with indices
                0-1, but when drawing with increasing x-axis in the
                vertical, going top to bottom we have TPCs 1-0, so we
                will define drawn_tpc that inverts them and use this
                index when dealing with the plane2 figure

                '''
                drawn_tpc = 1-tpc
                plt.figure(fig2)
                axs2[drawn_tpc].set_xlim(0,1663*0.3)
                axs2[drawn_tpc].set_ylim(0,3400*0.5*0.16)
                scaled_z = [ j * 0.3 for j in wires[3*tpc+plane]] # in cm
                if tpc == 0:
                    axs2[drawn_tpc].set_ylim(-3400*0.5*0.16,0)
                    scaled_x= [ -272+j*0.5*0.16 for j in times[3*tpc+plane]] # in cm
                if tpc == 1:
                    scaled_x = [ 272.-j*0.5*0.16 for j in times[3*tpc+plane]] # in cm
                axs2[drawn_tpc].scatter(scaled_z,scaled_x,c=charges[3*tpc+plane],marker=',',s=1)
                axs2[drawn_tpc].set_aspect(aspect=1.0)

                '''
                for each track draw deltaX and deltaZ and calculate the angle
                '''
                if draw_auxiliary_plane2:
                    for track in range(0,nmuontracks):
                        # grab smaller and biggest x and z for each track
                        if tpc == 0:
                            scaled_x = [ -272+times[3*tpc+plane][j]*0.5*0.16
                                        for j in range(0,len(times[3*tpc+plane]))
                                        if hit_track[3*tpc+plane][j] == track
                                       ]
                        if tpc == 1:
                            scaled_x = [ 272.-times[3*tpc+plane][j]*0.5*0.16
                                         for j in range(0,len(times[3*tpc+plane]))
                                         if hit_track[3*tpc+plane][j] == track
                                        ]
                        '''
                        if len(scaled_x) != 0:
                            xmin = np.amin(scaled_x)
                            xmax = np.amax(scaled_x)
                            if tpc == 0:
                                xmin = np.amax(scaled_x)
                                xmax = np.amin(scaled_x)
                        '''
                        scaled_z = [ wires[3*tpc+plane][j]* 0.3
                                     for j in range(0,len(wires[3*tpc+plane]))
                                     if hit_track[3*tpc+plane][j] == track
                                    ]
                        if len(scaled_z) != 0:
                            zmin = np.amin(scaled_z)
                            zmax = np.amax(scaled_z)
                            xmin = scaled_x[np.argmin(scaled_z)]
                            xmax = scaled_x[np.argmax(scaled_z)]
                        if len(scaled_z) == 0:
                            continue
                        line_h_x = [zmax,zmax]
                        line_v_x = [xmin,xmax]
                        axs2[drawn_tpc].plot(line_h_x,line_v_x,linestyle='--', color='Orange')
                        line_h_z = [zmin,zmax]
                        line_v_z = [xmin,xmin]
                        axs2[drawn_tpc].plot(line_h_z,line_v_z,linestyle='--', color='Orange')
                        deltax = xmax-xmin
                        deltaz = zmax-zmin
                        angle_thetaxz = math.degrees(math.atan2(deltax,deltaz))

                        # angle in degrees
                        axs2[drawn_tpc].text(0.55*(zmax+zmin),
                                             xmin+np.sign(deltax)*20,
                                             f"{angle_thetaxz:.2f}°",
                                             color='Red',
                                             fontweight='bold'
                                             )
                        # deltaz
                        axs2[drawn_tpc].text(0.5*(zmax+zmin),xmin-np.sign(deltax)*20,f"Δz = {zmax-zmin:.2f} cm",color='Orange')
                        # deltax
                        axs2[drawn_tpc].text(zmax+10,0.5*(xmin+xmax),f"Δx = {xmax-xmin:.2f} cm",color='Orange')
                    
                # switch back to main figure
                plt.figure(fig)
                
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
    axs2[1].set_title(f"TPC 0 ({tpc_tracks[0]} muontracks)",y=1.0,pad=-14)
    axs2[0].set_title(f"TPC 1 ({tpc_tracks[1]} muontracks)",y=1.0,pad=-14)
    #axs2[0].set_ylabel("X axis (cm)")
    fig2.text(0.12,0.5,"X coordinate (cm)", va='center',rotation='vertical')
    axs2[1].set_xlabel("Z coordinate (cm)")
    #axs[0][1].text(4,2.0,"random text",rotation=-90)
    
    output_dir = ""
    if test_only:
        output_dir = "img/test"
    else:
        os.mkdirs(f"full/{run}/{subrun}")
        output_dir = f"img/full/{run}/{subrun}"
    print(f"{output_dir}/mhits_evd_{e_counter}_{run}_{subrun}_{event}.png")
    
    plt.savefig(f"{output_dir}/mhits_evd_{e_counter}_{run}_{subrun}_{event}.png")
    plt.savefig(f"{output_dir}/mhits_evd_{e_counter}_{run}_{subrun}_{event}.pdf")
    plt.figure(fig2)
    plt.savefig(f"{output_dir}/mhits_evd_{e_counter}_{run}_{subrun}_{event}_plane2.png")
    plt.savefig(f"{output_dir}/mhits_evd_{e_counter}_{run}_{subrun}_{event}_plane2.pdf")
    e_counter+=1
    plt.close('all')

