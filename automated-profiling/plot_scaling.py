#!/usr/bin/env python
# script to read scaling data from either a json or csv file and produce scaling plots
# sim,tally,nparticles,xdim,ydim,nthreads,walltime,pathstime,tallytime


from os.path import exists
from numpy import logspace
import pathlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
import argparse
from flatten_json import flatten


class eiron_profile:
    def __init__( self, args ):
        self.pdata = self.read_data(args)


    def read_data(self, args):
        # Read data from a csv or json file file into a panda dataframe
        # If file type is specified then default to that, otherwise decide from extension
        if args.filetype:
            filetype = args.filetype
        else:
            filetype = pathlib.Path(args.file).suffix[1:]
        
        if filetype=='json':
            print("Reading json data from",args.file)
            try:
                data = self.read_json_with_pandas(args.file)
            except:
                print("ERROR: Cannot read input file")
                exit()
        elif filetype=='csv':
            print("Reading csv data from",args.file)
            try:
                data = self.read_flat_csv_with_pandas(args.file)
            except:
                print("ERROR: Cannot read input file")
                exit()
        else:
            print("Error: Unknown file type",filetype)
            exit()

        # Add a number of processes column and sort by threads, processes and nodes 
        data['n_processes'] = data['n_nodes'] * data['n_mpi_ranks']            
        return data[data['job_state']=='COMPLETED'].sort_values(by=['n_omp_threads','n_processes','n_nodes'])
            
        
    # Read profile data from a csv file into a dataframe
    def read_flat_csv_with_pandas(self, file):
        self.filetype = "flat csv"
        data = pd.read_csv(file)
        return data

    # Read the data from a json file into a dataframe
    def read_json_with_pandas(self, file):
        self.filetype = 'json'
        self.filestruct = 'nested'
        data = pd.read_json(file)
        # Flatten
        df = pd.json_normalize(data['cases'])
        return df

    # Filters data on the given values from flat csv data, 0 means all values
    def filter_data_flat_csv(self,sim,tally,nparticles,xdim,ydim,nthreads):
        filtersim   = self.pdata['sim']==sim
        filtertally = self.pdata['tally']==tally
        filterpart  = (self.pdata['nparticles']==nparticles) | (nparticles == 0)
        filterxdim  = (self.pdata['xdim']==xdim) | (xdim == 0)
        filterydim  = (self.pdata['ydim']==ydim) | (ydim == 0)
        filternthreads  = (self.pdata['nthreads']==nthreads) | (nthreads == 0)
        return self.pdata[ filtersim & filtertally & filterpart & filterxdim & filterydim & filternthreads]


    # Filter data for plotting based on a dict containing the key values
    def filter_data(self, **filter_dict):
        keys = filter_dict.keys()
        filter = True        
        for key in keys:
            if key not in self.pdata.keys() :
                print("ERROR: Key value does not exist! Exiting")
                exit()
            filter = (self.pdata[key]==filter_dict[key]) & filter
        return self.pdata[ filter ];
    
    ###############################################################################################
    # Plotting routines
    ###############################################################################################
    
    # Set up the plotting environment with some default values
    def init_figure(self, title):
        self.fig, self.ax = plt.subplots(figsize=(6,4))#, layout='constrained')
        # Set some defaults before calling plot
        self.ax.set_xlabel('Number of cores')
        self.ax.set_ylabel('Execution time (s)')
        self.ax.set_xscale('log')
        self.ax.set_yscale('log')
        self.ax.set_xticks(logspace(0,7,8,base=2,dtype=int))
        self.ax.set_xticklabels(logspace(0,7,8,base=2,dtype=int))
        self.ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1f}"))
        self.ax.xaxis.grid()
        self.ax.yaxis.grid()
        self.ax.xaxis.grid(which='minor')
        self.ax.yaxis.grid(which='minor')
        self.ax.set_title(title)
        
    # Plot the data to the current figure
    def plot_data(self, xdata, ydata, param_dict):
        # Plot the data
        self.ax.plot(xdata,ydata,'o-',**param_dict)

    # Decide on the plot output format    
    def makeplot(self,filename):
        plt.savefig(filename)
        #plt.show()

    # Routines for producing particular figures

    # A simple figure for eirene strong thread scaling
    def fig_simple_strong_thread_scaling(self,nparticles):
        nparticles = 11000
        print('Plotting simple strong thread scaling:',nparticles,' particles')
        self.init_figure('Eirene Strong Thread Scaling')
        fdata = self.filter_data(n_mpi_ranks=1,n_nodes=1)
        self.plot_data(fdata['n_omp_threads'],fdata['timing.wall_time'],{'label':'walltime'})
        self.ax.set_xlabel('Number of OpenMP Threads')
        self.ax.legend()
        self.makeplot("plots/strongThreadScaling.png")

    # A simple figure for eirene strong mpi scaling
    def fig_simple_strong_mpi_scaling(self,nparticles,nthreads):
        nparticles = 11000
        print('Plotting simple strong MPI scaling:',nparticles,' particles')
        self.init_figure('Eirene Strong MPI Scaling')
        fdata = self.filter_data(n_omp_threads=nthreads).drop_duplicates('n_processes')
        self.plot_data(fdata['n_processes'],fdata['timing.wall_time'],{'label':'walltime'})
        self.ax.set_xlabel('Number of MPI processes')
        self.ax.legend()
        self.makeplot("plots/strongMPIScaling.png")

    # Eirene strong mpi scaling for different thread counts
    def fig_strong_mpi_scaling(self,nparticles):
        nparticles = 11000
        print('Plotting strong MPI scaling with threads:',nparticles,' particles')
        self.init_figure('Eirene Strong MPI Scaling')
        allthreads = self.pdata['n_omp_threads'].drop_duplicates()
        for threads in allthreads:
            fdata = self.filter_data(n_omp_threads=threads).drop_duplicates('n_processes')
            self.plot_data(fdata['n_processes'],fdata['timing.wall_time'],{'label':str(threads)+' threads'})
        self.ax.set_xlabel('Number of MPI processes')
        self.ax.legend()
        self.makeplot("plots/strongMPIScalingThreads.png")
        
    # Add strong scaling plots for the synchronous execution types        
    def plot_one_strong_scaling(self,nparticles,xgrid,ygrid,sim,tally,data):
        fdata = self.filter_data_flat_csv( sim, tally, nparticles, xgrid, ygrid, 0)
        self.plot_data(fdata['nthreads'],fdata[data],{'label':sim + ' ' + tally})

    # Add strong scaling plots showing total time of the 4 execution types        
    def plot_strong_scaling(self,nparticles,xgrid,ygrid):
        self.plot_one_strong_scaling(nparticles,xgrid,ygrid,'monolithic','private','walltime')
        self.plot_one_strong_scaling(nparticles,xgrid,ygrid,'monolithic','shared','walltime')
        self.plot_one_strong_scaling(nparticles,xgrid,ygrid,'synchronous','private','walltime')
        self.plot_one_strong_scaling(nparticles,xgrid,ygrid,'synchronous','shared','walltime')

    # Strong scaling showing total time of the 4 execution types
    def fig_strong_scaling(self,nparticles,xgrid,ygrid):
        print('Plotting strong scaling:',nparticles,' particles')
        self.init_figure('Strong scaling: ' + str(nparticles) + ' particles, ' + str(xgrid) + 'x' + str(ygrid) + ' grid')
        self.plot_strong_scaling(nparticles,xgrid,ygrid)
        self.ax.legend()
        self.makeplot("strongScaling.png")
    
    # Add strong scaling plots for the synchronous execution types        
    def plot_one_strong_scaling_split(self,nparticles,xgrid,ygrid,tally,data):
        fdata = self.filter_data_flat_csv( 'synchronous', tally, nparticles, xgrid, ygrid, 0)
        self.plot_data(fdata['nthreads'],fdata[data],{'label':data + ' ' + tally})       

    # Add strong scaling plots showing paths and tally time
    def plot_strong_scaling_split(self,nparticles,xgrid,ygrid):
        self.plot_one_strong_scaling_split(nparticles,xgrid,ygrid,'private','pathstime')
        self.plot_one_strong_scaling_split(nparticles,xgrid,ygrid,'shared','pathstime')
        self.plot_one_strong_scaling_split(nparticles,xgrid,ygrid,'private','tallytime')
        self.plot_one_strong_scaling_split(nparticles,xgrid,ygrid,'shared','tallytime')
        
    # Strong scaling showing total time of the 4 execution types
    def fig_strong_scaling_split(self,nparticles,xgrid,ygrid):
        print('Plotting strong scaling split:',nparticles,'particles')
        self.init_figure('Strong scaling split: ' + str(nparticles) + ' particles, ' + str(xgrid) + 'x' + str(ygrid) + ' grid')
        self.plot_strong_scaling_split(nparticles,xgrid,ygrid)
        self.ax.legend()
        self.makeplot("strongScalingSplit.png")
    
    
    # Get data for weak scaling of total execution time with particles
    # This relies on the data increasing as a power of 2
    def get_weak_scaling_particles(self,sim,tally,firstthread,xgrid,ygrid):
        thread = firstthread
        particles = self.filter_data_flat_csv( sim, tally, 0, xgrid, ygrid, 1)['nparticles'].sort_values()
        fdata = self.filter_data_flat_csv( sim, tally, particles.min(), xgrid, ygrid, thread)
        for npart in particles.iloc[1:]:
            thread = 2 * thread
            fdata1 = self.filter_data_flat_csv( sim, tally, npart, xgrid, ygrid, thread)            
            fdata = pd.concat([fdata,fdata1],axis=0)
        return fdata

    # Get data for weak scaling of total execution time with grid size
    def plot_weak_scaling_particles(self,firstthread,xgrid,ygrid):
        fdata = self.get_weak_scaling_particles('monolithic','private',firstthread,xgrid,ygrid)
        self.plot_data(fdata['nthreads'],fdata['walltime'],{'label':'Monolithic Private'})
        fdata = self.get_weak_scaling_particles('monolithic','shared',firstthread,xgrid,ygrid)
        self.plot_data(fdata['nthreads'],fdata['walltime'],{'label':'Monolithic Shared'})
        fdata = self.get_weak_scaling_particles('synchronous','private',firstthread,xgrid,ygrid)
        self.plot_data(fdata['nthreads'],fdata['walltime'],{'label':'Synchronous Private'})
        fdata = self.get_weak_scaling_particles('synchronous','shared',firstthread,xgrid,ygrid)
        self.plot_data(fdata['nthreads'],fdata['walltime'],{'label':'Synchronous Shared'})

    # Make figure of weak scaling with particles
    def fig_weak_scaling_particles(self,firstthread,xgrid,ygrid):
        print('Plotting weak scaling with particles:',xgrid,'x',ygrid,'grid')
        self.init_figure('Weak scaling with particles. '+str(xgrid)+'x'+str(ygrid)+' grid')
        self.plot_weak_scaling_particles(firstthread,xgrid,ygrid)
        minpart = self.pdata['nparticles'].min()
        self.ax.text(2,3.9,'From %8d particles' %(minpart))
        self.ax.legend()
        self.ax.set_yscale('linear')
        self.makeplot('weakScalingParticles.png')

    # Get data for weak scaling with grid size for pathstime and tallytime
    def plot_weak_scaling_split_particles(self,firstthread,xgrid,ygrid):
        fdata = self.get_weak_scaling_particles('synchronous','private',firstthread,xgrid,ygrid)
        self.plot_data(fdata['nthreads'],fdata['pathstime'],{'label':'pathstime private'})
        fdata = self.get_weak_scaling_particles('synchronous','shared',firstthread,xgrid,ygrid)
        self.plot_data(fdata['nthreads'],fdata['pathstime'],{'label':'pathstime shared'})
        fdata = self.get_weak_scaling_particles('synchronous','private',firstthread,xgrid,ygrid)
        self.plot_data(fdata['nthreads'],fdata['tallytime'],{'label':'tallytime private'})
        fdata = self.get_weak_scaling_particles('synchronous','shared',firstthread,xgrid,ygrid)
        self.plot_data(fdata['nthreads'],fdata['tallytime'],{'label':'tallytime shared'})

    # Make figure of weak scaling with particles for pathstime and tallytime
    def fig_weak_scaling_split_particles(self,firstthread,xgrid,ygrid):
        print('Plotting weak scaling split with particles:',xgrid,'x',ygrid,'grid')
        self.init_figure('Weak scaling split with particles. '+str(xgrid)+'x'+str(ygrid)+' grid')
        self.plot_weak_scaling_split_particles(firstthread,xgrid,ygrid)
        minpart = self.pdata['nparticles'].min()
        self.ax.text(2,2.5,'From %8d particles' %(minpart))
        self.ax.legend()
        self.ax.set_yscale('linear')
        self.makeplot('weakScalingSplitParticles.png')
    
    # Get data for weak scaling of total execution time with grid size
    # This relies on the data increasing as a power of 2
    def get_weak_scaling_grid(self,sim,tally,firstthread,nparticles):
        thread = firstthread
        xgrid = self.filter_data_flat_csv( sim, tally, nparticles, 0, 0, 1)['xdim'].sort_values()
        fdata = self.filter_data_flat_csv( sim, tally, nparticles, xgrid.min(), xgrid.min(), thread)
        for grid in xgrid.iloc[1:]:
            thread = 4 * thread
            fdata1 = self.filter_data_flat_csv( sim, tally, nparticles, grid, grid, thread)            
            fdata = pd.concat([fdata,fdata1],axis=0)
        return fdata

    # Add plots for weak scaling with grid size
    def plot_weak_scaling_grid(self,firstthread,nparticles):
        fdata = self.get_weak_scaling_grid('monolithic','private',firstthread,nparticles)
        self.plot_data(fdata['nthreads'],fdata['walltime'],{'label':'Monolithic Private'})
        fdata = self.get_weak_scaling_grid('monolithic','shared',firstthread,nparticles)
        self.plot_data(fdata['nthreads'],fdata['walltime'],{'label':'Monolithic Shared'})
        fdata = self.get_weak_scaling_grid('synchronous','private',firstthread,nparticles)
        self.plot_data(fdata['nthreads'],fdata['walltime'],{'label':'Synchronous Private'})
        fdata = self.get_weak_scaling_grid('synchronous','shared',firstthread,nparticles)
        self.plot_data(fdata['nthreads'],fdata['walltime'],{'label':'Synchronous Shared'})

    # Make figure of weak scaling with grid size
    def fig_weak_scaling_grid(self,firstthread,nparticles):
        print('Plotting weak scaling with grid size:', nparticles, ' particles' )
        self.init_figure('Weak scaling with grid size. '+ str(nparticles) + ' particles')
        self.plot_weak_scaling_grid(firstthread,nparticles)
        self.ax.set_yscale('linear')
        mingrid = self.pdata['xdim'].min()
        self.ax.annotate('Smallest grid=%4d' %(mingrid), xy=(1.7,3.4)) 
        self.ax.legend()
        self.makeplot('weakScalingGrid.png')

    # Add plots for weak scaling with grid size for pathstime and tallytime
    def plot_weak_scaling_split_grid(self,firstthread,nparticles):
        fdata = self.get_weak_scaling_grid('synchronous','private',firstthread,nparticles)
        self.plot_data(fdata['nthreads'],fdata['pathstime'],{'label':'pathstime private'})
        fdata = self.get_weak_scaling_grid('synchronous','shared',firstthread,nparticles)
        self.plot_data(fdata['nthreads'],fdata['pathstime'],{'label':'pathstime shared'})
        fdata = self.get_weak_scaling_grid('synchronous','private',firstthread,nparticles)
        self.plot_data(fdata['nthreads'],fdata['tallytime'],{'label':'tallytime private'})
        fdata = self.get_weak_scaling_grid('synchronous','shared',firstthread,nparticles)
        self.plot_data(fdata['nthreads'],fdata['tallytime'],{'label':'tallytime shared'})

    # Make figure of weak scaling with grid size for pathstime and tallytime
    def fig_weak_scaling_split_grid(self,firstthread,nparticles):
        print('Plotting weak scaling with grid size:', nparticles, ' particles' )
        self.init_figure('Weak scaling split with grid size. '+ str(nparticles) + ' particles')
        self.plot_weak_scaling_split_grid(firstthread,nparticles)
        self.ax.set_yscale('linear')
        # This won't work in general
        mingrid = self.pdata['xdim'].min()
        self.ax.text(2,0.6,'Smallest grid=%4d' %(mingrid)) 
        self.ax.legend()
        self.makeplot('weakScalingSplitGrid.png')

    # Add plots of speedup
    def plot_speedup(self,nparticles,xgrid,ygrid):
        fdata = self.filter_data_flat_csv('monolithic','private', nparticles, xgrid, ygrid, 0)
        speedup = fdata['walltime'].iloc[0]/fdata['walltime']
        self.plot_data(fdata['nthreads'],fdata['nthreads'],{'label':'ideal'})
        self.plot_data(fdata['nthreads'],speedup,{'label':'monolithic private'})
        fdata = self.filter_data_flat_csv('monolithic','shared', nparticles, xgrid, ygrid, 0)
        speedup = fdata['walltime'].iloc[0]/fdata['walltime']
        self.plot_data(fdata['nthreads'],speedup,{'label':'monolithic shared'})
        fdata = self.filter_data_flat_csv('synchronous','private', nparticles, xgrid, ygrid, 0)
        speedup = fdata['walltime'].iloc[0]/fdata['walltime']
        self.plot_data(fdata['nthreads'],speedup,{'label':'synchronous private'})
        fdata = self.filter_data_flat_csv('synchronous','shared', nparticles, xgrid, ygrid, 0)
        speedup = fdata['walltime'].iloc[0]/fdata['walltime']
        self.plot_data(fdata['nthreads'],speedup,{'label':'synchronous shared'})


    # Make figure of speedup
    def fig_speedup(self,nparticles,xgrid,ygrid):
        print('Plotting speedup:', nparticles, 'particles')
        self.init_figure('OpenMP speedup: ' + str(nparticles) + ' particles, ' + str(xgrid) + 'x' + str(ygrid) + ' grid')
        self.plot_speedup(nparticles,xgrid,ygrid)
        self.ax.set_ylabel('Speedup')
        self.ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1f}"))
        self.ax.legend()
        self.makeplot('openMPSpeedup.png')

    # Add plots of efficiency
    def plot_efficiency(self,nparticles,xgrid,ygrid):
        fdata = self.filter_data_flat_csv('monolithic','private', nparticles, xgrid, ygrid, 0)
        efficiency = fdata['walltime'].iloc[0] / (fdata['walltime'] * fdata['nthreads'])
        self.plot_data(fdata['nthreads'],efficiency,{'label':'monolithic private'})
        fdata = self.filter_data_flat_csv('monolithic','shared', nparticles, xgrid, ygrid, 0)
        efficiency = fdata['walltime'].iloc[0] / (fdata['walltime'] * fdata['nthreads'])
        self.plot_data(fdata['nthreads'],efficiency,{'label':'monolithic shared'})
        fdata = self.filter_data_flat_csv('synchronous','private', nparticles, xgrid, ygrid, 0)
        efficiency = fdata['walltime'].iloc[0] / (fdata['walltime'] * fdata['nthreads'])
        self.plot_data(fdata['nthreads'],efficiency,{'label':'synchronous private'})
        fdata = self.filter_data_flat_csv('synchronous','shared', nparticles, xgrid, ygrid, 0)
        efficiency = fdata['walltime'].iloc[0] / (fdata['walltime'] * fdata['nthreads'])
        self.plot_data(fdata['nthreads'],efficiency,{'label':'synchronous shared'})


    # Make figure of efficiency
    def fig_efficiency(self,nparticles,xgrid,ygrid):
        print('Plotting efficiency:', nparticles, 'particles,', xgrid ,'x',ygrid,'grid')
        self.init_figure('OpenMP efficiency: ' + str(nparticles) + ' particles, ' + str(xgrid) + 'x' + str(ygrid) + ' grid')
        self.plot_efficiency(nparticles,xgrid,ygrid)
        self.ax.set_ylabel('Efficiency')
        self.ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1f}"))
        self.ax.legend()
        self.makeplot('openMPEfficiency.png')

# Read command line arguments
def read_args():
    parser = argparse.ArgumentParser(description = 'Plot scalings from csv file containing scaling data.')
    parser.add_argument( "-f" ,"--file" , default="./SCALABILITY_REPORT/report_2D-D_slab.json" , help='csv file containing eiron profiling data (default=eiron_profile.csv)' )
    parser.add_argument( "-n" ,"--nparticles" , default=800000 , help='Number of particles for scaling plots' )
    parser.add_argument( "-g" ,"--gridsize" , default=256 , help='Grid size for scaling plots' )
    parser.add_argument( "--firstthread" , default=2 , help='First thread number for use in weak scaling plots' )
    parser.add_argument( "-t" ,"--filetype" , choices=['csv','json'], help='File type conataining profile data, either csv or json. Overrides extension checking' )    
    args = parser.parse_args()
    if not exists(args.file):
        print("Error:", args.file, "does not exist, exiting")
        exit()
    return args
    
        
print("Plotting scaling profiles")

args = read_args()

# Create object from from csv file
prof = eiron_profile( args )

pd.options.display.max_rows = 999
pd.options.display.max_columns = 999
#print(self.pdata)

# Make some plots
prof.fig_simple_strong_thread_scaling(args.nparticles)
prof.fig_simple_strong_mpi_scaling(args.nparticles, 1)
prof.fig_strong_mpi_scaling(args.nparticles)
#prof.fig_strong_scaling(args.nparticles,args.gridsize,args.gridsize)
#prof.fig_strong_scaling_split(args.nparticles,args.gridsize,args.gridsize)
#prof.fig_weak_scaling_particles(args.firstthread,args.gridsize,args.gridsize)
#prof.fig_weak_scaling_split_particles(args.firstthread,args.gridsize,args.gridsize)
#prof.fig_weak_scaling_grid(args.firstthread,args.nparticles)
#prof.fig_weak_scaling_split_grid(args.firstthread,args.nparticles)
#prof.fig_speedup(args.nparticles,args.gridsize,args.gridsize)
#prof.fig_efficiency(args.nparticles,args.gridsize,args.gridsize)





