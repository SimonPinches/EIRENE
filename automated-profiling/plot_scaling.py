#!/usr/bin/env python
# Authour H Leggate
# Script to read scaling data from either a json or csv file and produce scaling plots
# Data can be read in a flat csv file or a hierarchical json file
# There are a set of plots for data from either Eiron or Eirene

# Eiron data contains the following indexes
# sim,tally,nparticles,xdim,ydim,nthreads,walltime,pathstime,tallytime

# Eirene data contains the following indexes
# n_nodes,n_mpi_ranks,n_omp_threads 


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
        codetype = args.codetype
        self.pdata = self.read_data(args)


    def read_data(self, args):
        # Read data from a csv or json file file into a panda dataframe
        # If file type is specified then default to that, otherwise decide from extension
        if args.filetype:
            filetype = args.filetype
        else:
            filetype = pathlib.Path(args.file).suffix[1:]
        
        if filetype=='json':
            print('Reading json data from',args.file)
            try:
                data = self.read_json_with_pandas(args.file)
            except:
                print('ERROR: Cannot read json input file', args.file)
                exit()
        elif filetype=='csv':
            print('Reading csv data from',args.file)
            try:
                data = self.read_flat_csv_with_pandas(args.file)
            except:
                print('ERROR: Cannot read csv input file')
                exit()
        else:
            print('Error: Unknown file type',filetype)
            exit()

        return data

    def strplural(self,number,val):
        if number > 1:
            val=val+'s'
        return str(number) + ' ' + val
    
        
    # Read profile data from a csv file into a dataframe
    def read_flat_csv_with_pandas(self, file):
        self.filetype = "csv"
        self.filetype = "flat"
        data = pd.read_csv(file)
        return data

    # Read the data from a json file into a dataframe
    def read_json_with_pandas(self, file):
        self.filetype = 'json'
        self.filestruct = 'nested'
        data = pd.read_json(file)
        # Flatten
        data = pd.json_normalize(data['cases'])
        # Add a number of processes column and sort by threads, processes and nodes 
        data['n_processes'] = data['n_nodes'] * data['n_mpi_ranks']            
        return data[data['job_state']=='COMPLETED'].sort_values(by=['n_omp_threads','n_processes','n_nodes'])

    # Filter data for plotting based on a dict containing the key values
    def filter_data(self, **filter_dict):
        keys = filter_dict.keys()
        filter = True        
        for key in keys:
            if key not in self.pdata.keys() :
                print('ERROR: Key value',key,'does not exist! Exiting')
                exit()
            filter = (self.pdata[key]==filter_dict[key]) & filter
        return self.pdata[ filter ];
    
    ###############################################################################################
    # Generic plotting routines
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

    ################################################################################################
    # Routines for producing particular figures
    # There are different figures for Eirene and Eiron
    # it is possible to make generic routines but there seems little benefit
    # as they are so quick to write
    ################################################################################################

    ################################################################################################
    # Eirene figures
    ################################################################################################
    
    # A simple figure for eirene strong thread scaling
    def omp_simple_strong_scaling_eirene(self,nparticles):
        print('Plotting Eirene simple OpenMP strong scaling:',self.strplural(nparticles,'particle'))
        self.init_figure('Eirene OpenMP Strong Scaling')
        fdata = self.filter_data(n_mpi_ranks=1,n_nodes=1)
        self.plot_data(fdata['n_omp_threads'],fdata['timing.wall_time'],{'label':'walltime'})
        self.ax.set_xlabel('Number of OpenMP Threads')
        self.ax.legend()
        self.makeplot("plots/eireneOpenMPStrongScaling.png")

    # A simple figure for eirene strong mpi scaling
    def mpi_simple_strong_scaling_eirene(self,nparticles,nthreads):
        print('Plotting Eirene simple MPI strong scaling:',self.strplural(nparticles,'particle'),self.strplural(nthreads,'thread'))
        self.init_figure('Eirene Strong MPI Scaling')
        fdata = self.filter_data(n_omp_threads=nthreads).drop_duplicates('n_processes')
        self.plot_data(fdata['n_processes'],fdata['timing.wall_time'],{'label':'walltime'})
        self.ax.set_xlabel('Number of MPI processes')
        self.ax.legend()
        self.makeplot('plots/eireneMPIStrongScaling_'+str(nthreads)+'threads.png')

    # Eirene strong mpi scaling for different thread counts
    def hybrid_strong_scaling_eirene(self,nparticles):
        print('Plotting Eirene hybrid MPI OpenMP strong scaling with threads:',self.strplural(nparticles,'particle'))
        self.init_figure('Eirene Hybrid Strong Scaling')
        allthreads = self.pdata['n_omp_threads'].drop_duplicates()
        for threads in allthreads:
            fdata = self.filter_data(n_omp_threads=threads).drop_duplicates('n_processes')
            self.plot_data(fdata['n_processes'],fdata['timing.wall_time'],{'label':self.strplural(threads,'thread')})
        self.ax.set_xlabel('Number of MPI processes')
        self.ax.legend()
        self.makeplot("plots/eireneHybridStrongScaling.png")

    # Strong OpenMP scaling showing total time and the time for each species
    def omp_strong_scaling_split_eirene(self,nparticles):
        print('Plotting Eirene OpenMP strong scaling split:',self.strplural(nparticles,'particle'))
        self.init_figure('OpenMP strong scaling split: ' + self.strplural(nparticles,'particle'))
        fdata = self.filter_data(n_mpi_ranks=1,n_nodes=1)
        self.plot_data(fdata['n_omp_threads'],fdata['timing.wall_time'],{'label':'walltime'})
        self.plot_data(fdata['n_omp_threads'],fdata['timing.atom_time'],{'label':'atoms'})
        self.plot_data(fdata['n_omp_threads'],fdata['timing.mol_time'],{'label':'molecules'})
        self.plot_data(fdata['n_omp_threads'],fdata['timing.ion_time'],{'label':'ions'})
        self.ax.set_xlabel('Number of OpenMP threads')
        self.ax.legend()
        self.makeplot('plots/eireneOpenMPStrongScalingSplit.png')
        
    # Strong MPI scaling showing total time and the time for each species
    def mpi_strong_scaling_split_eirene(self,nparticles,nthreads):
        print('Plotting Eirene MPI strong scaling split:',self.strplural(nparticles,'particle'),self.strplural(nthreads,'thread'))
        self.init_figure('Strong scaling split: ' + self.strplural(nparticles,'particle') + ', ' + self.strplural(nthreads,'thread'))
        fdata = self.filter_data(n_omp_threads=nthreads).drop_duplicates('n_processes')
        self.plot_data(fdata['n_processes'],fdata['timing.wall_time'],{'label':'walltime'})
        self.plot_data(fdata['n_processes'],fdata['timing.atom_time'],{'label':'atoms'})
        self.plot_data(fdata['n_processes'],fdata['timing.mol_time'],{'label':'molecules'})
        self.plot_data(fdata['n_processes'],fdata['timing.ion_time'],{'label':'ions'})
        self.ax.set_xlabel('Number of MPI processes')
        self.ax.legend()
        self.makeplot('plots/eireneMPIStrongScalingSplit_'+str(nthreads)+'threads.png')

    # Make figure of Eirene MPI speedup
    def mpi_speedup_eirene(self,nparticles):
        print('Plotting Eirene MPI speedup:', self.strplural(nparticles,'particle'))
        self.init_figure('Eirene MPI speedup: ' + str(nparticles) + ' particles')
        allthreads = self.pdata['n_omp_threads'].drop_duplicates()
        for threads in allthreads:
            fdata = self.filter_data(n_omp_threads=threads).drop_duplicates('n_processes')
            speedup =  fdata['timing.wall_time'].iloc[0]/fdata['timing.wall_time']
            if threads == 1:
                self.plot_data(fdata['n_processes'],fdata['n_processes'],{'label':'ideal'})
            self.plot_data(fdata['n_processes'],speedup,{'label':self.strplural(threads,'thread')})
        self.ax.set_ylabel('Speedup')
        self.ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1f}"))
        self.ax.legend()
        self.makeplot('plots/eireneMPISpeedup.png')

    # Make figure of Eirene OpenMP speedup
    def omp_speedup_eirene(self,nparticles):
        print('Plotting Eirene OpenMP speedup:', self.strplural(nparticles,'particle'))
        self.init_figure('Eirene OpenMP speedup: ' + str(nparticles) + ' particles')
        fdata = self.filter_data(n_nodes=1,n_mpi_ranks=1)
        speedup = fdata['timing.wall_time'].iloc[0]/fdata['timing.wall_time']
        self.plot_data(fdata['n_omp_threads'],fdata['n_omp_threads'],{'label':'ideal'})
        self.plot_data(fdata['n_omp_threads'],speedup,{'label':'walltime'})
        self.ax.set_ylabel('Speedup')
        self.ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1f}"))
        self.ax.legend()
        self.makeplot('plots/eireneOpenMPSpeedup.png')

    # Make figure of Eirene MPI efficiency
    def mpi_efficiency_eirene(self,nparticles):
        print('Plotting Eirene MPI efficiency:', self.strplural(nparticles,'particle'))
        self.init_figure('MPI efficiency: ' + self.strplural(nparticles,'particle'))
        allthreads = self.pdata['n_omp_threads'].drop_duplicates()
        for threads in allthreads:
            fdata = self.filter_data(n_omp_threads=threads).drop_duplicates('n_processes')
            efficiency =  fdata['timing.wall_time'].iloc[0]/(fdata['timing.wall_time']*fdata['n_processes'])
            self.plot_data(fdata['n_processes'],efficiency,{'label':self.strplural(threads,'thread')})
        self.ax.set_ylabel('Efficiency')
        self.ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1f}"))
        self.ax.legend()
        self.makeplot('plots/eireneMPIEfficiency.png')

    # Make figure of Eirene OpenMP efficiency
    def omp_efficiency_eirene(self,nparticles):
        print('Plotting Eirene OpenMP efficiency:', self.strplural(nparticles,'particle'))
        self.init_figure('OpenMP efficiency: ' + self.strplural(nparticles,'particle'))
        fdata = self.filter_data(n_nodes=1,n_mpi_ranks=1)
        efficiency = fdata['timing.wall_time'].iloc[0]/(fdata['timing.wall_time']*fdata['n_omp_threads'])
        self.plot_data(fdata['n_omp_threads'],efficiency,{'label':'walltime'})
        self.ax.set_ylabel('Efficiency')
        self.ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1f}"))
        self.ax.legend()
        self.makeplot('plots/eireneOpenMPEfficiency.png')


    ##############################################################################################
    # Eiron figures
    ##############################################################################################

        
    # Add strong scaling plots for the synchronous execution types        
    def plot_one_strong_scaling_eiron(self,nparticles,xgrid,ygrid,sim,tally,data):
        fdata = self.filter_data( sim=sim, tally=tally, nparticles=nparticles, xdim=xgrid, ydim=ygrid)
        self.plot_data(fdata['nthreads'],fdata[data],{'label':sim + ' ' + tally})

    # Add strong scaling plots showing total time of the 4 execution types        
    def plot_strong_scaling_eiron(self,nparticles,xgrid,ygrid):
        self.plot_one_strong_scaling(nparticles,xgrid,ygrid,'monolithic','private','walltime')
        self.plot_one_strong_scaling(nparticles,xgrid,ygrid,'monolithic','shared','walltime')
        self.plot_one_strong_scaling(nparticles,xgrid,ygrid,'synchronous','private','walltime')
        self.plot_one_strong_scaling(nparticles,xgrid,ygrid,'synchronous','shared','walltime')

    # Strong scaling showing total time of the 4 execution types
    def fig_strong_scaling_eiron(self,nparticles,xgrid,ygrid):
        print('Plotting strong scaling:',nparticles,' particles')
        self.init_figure('Strong scaling: ' + str(nparticles) + ' particles, ' + str(xgrid) + 'x' + str(ygrid) + ' grid')
        self.plot_strong_scaling(nparticles,xgrid,ygrid)
        self.ax.legend()
        self.makeplot("strongScaling.png")
    

    # Add strong scaling plots for the synchronous execution types        
    def plot_one_strong_scaling_split_eiron(self,nparticles,xgrid,ygrid,tally,data):
        fdata = self.filter_data( sim='synchronous', tally=tally, nparticles=nparticles, xdim=xgrid, ydim=ygrid)
        self.plot_data(fdata['nthreads'],fdata[data],{'label':data + ' ' + tally})       

    # Add strong scaling plots showing paths and tally time
    def plot_strong_scaling_split_eiron(self,nparticles,xgrid,ygrid):
        self.plot_one_strong_scaling_split(nparticles,xgrid,ygrid,'private','pathstime')
        self.plot_one_strong_scaling_split(nparticles,xgrid,ygrid,'shared','pathstime')
        self.plot_one_strong_scaling_split(nparticles,xgrid,ygrid,'private','tallytime')
        self.plot_one_strong_scaling_split(nparticles,xgrid,ygrid,'shared','tallytime')
        
    # Strong scaling showing total time of the 4 execution types
    def fig_strong_scaling_split_eiron(self,nparticles,xgrid,ygrid):
        print('Plotting strong scaling split:',nparticles,'particles')
        self.init_figure('Strong scaling split: ' + str(nparticles) + ' particles, ' + str(xgrid) + 'x' + str(ygrid) + ' grid')
        self.plot_strong_scaling_split(nparticles,xgrid,ygrid)
        self.ax.legend()
        self.makeplot("strongScalingSplit.png")
    
    
    # Get data for weak scaling of total execution time with particles
    # This relies on the data increasing as a power of 2
    def get_weak_scaling_particles_eiron(self,sim,tally,firstthread,xgrid,ygrid):
        thread = firstthread
        particles = self.filter_data( sim=sim, tally=tally, xdim=xgrid, ydim=ygrid, nthreads=1)['nparticles'].sort_values()
        fdata = self.filter_data( sim=sim, tally=tally, nparticles=particles.min(), xdim=xgrid, ydim=ygrid, nthreads=thread)
        for npart in particles.iloc[1:]:
            thread = 2 * thread
            fdata1 = self.filter_data( sim=sim, tally=tally, nparticles=npart, xdim=xgrid, ydim=ygrid, nthreads=thread)            
            fdata = pd.concat([fdata,fdata1],axis=0)
        return fdata

    # Get data for weak scaling of total execution time with grid size
    def plot_weak_scaling_particles_eiron(self,firstthread,xgrid,ygrid):
        fdata = self.get_weak_scaling_particles('monolithic','private',firstthread,xgrid,ygrid)
        self.plot_data(fdata['nthreads'],fdata['walltime'],{'label':'Monolithic Private'})
        fdata = self.get_weak_scaling_particles('monolithic','shared',firstthread,xgrid,ygrid)
        self.plot_data(fdata['nthreads'],fdata['walltime'],{'label':'Monolithic Shared'})
        fdata = self.get_weak_scaling_particles('synchronous','private',firstthread,xgrid,ygrid)
        self.plot_data(fdata['nthreads'],fdata['walltime'],{'label':'Synchronous Private'})
        fdata = self.get_weak_scaling_particles('synchronous','shared',firstthread,xgrid,ygrid)
        self.plot_data(fdata['nthreads'],fdata['walltime'],{'label':'Synchronous Shared'})

    # Make figure of weak scaling with particles
    def fig_weak_scaling_particles_eiron(self,firstthread,xgrid,ygrid):
        print('Plotting weak scaling with particles:',xgrid,'x',ygrid,'grid')
        self.init_figure('Weak scaling with particles. '+str(xgrid)+'x'+str(ygrid)+' grid')
        self.plot_weak_scaling_particles(firstthread,xgrid,ygrid)
        minpart = self.pdata['nparticles'].min()
        self.ax.text(2,3.9,'From %8d particles' %(minpart))
        self.ax.legend()
        self.ax.set_yscale('linear')
        self.makeplot('weakScalingParticles.png')


    # Get data for weak scaling with grid size for pathstime and tallytime
    def plot_weak_scaling_split_particles_eiron(self,firstthread,xgrid,ygrid):
        fdata = self.get_weak_scaling_particles('synchronous','private',firstthread,xgrid,ygrid)
        self.plot_data(fdata['nthreads'],fdata['pathstime'],{'label':'pathstime private'})
        fdata = self.get_weak_scaling_particles('synchronous','shared',firstthread,xgrid,ygrid)
        self.plot_data(fdata['nthreads'],fdata['pathstime'],{'label':'pathstime shared'})
        fdata = self.get_weak_scaling_particles('synchronous','private',firstthread,xgrid,ygrid)
        self.plot_data(fdata['nthreads'],fdata['tallytime'],{'label':'tallytime private'})
        fdata = self.get_weak_scaling_particles('synchronous','shared',firstthread,xgrid,ygrid)
        self.plot_data(fdata['nthreads'],fdata['tallytime'],{'label':'tallytime shared'})

    # Make figure of weak scaling with particles for pathstime and tallytime
    def fig_weak_scaling_split_particles_eiron(self,firstthread,xgrid,ygrid):
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
    def get_weak_scaling_grid_eiron(self,sim,tally,firstthread,nparticles):
        thread = firstthread
        xgrid = self.filter_data( sim=sim, tally=tally, nparticles=nparticles, nthreads=1)['xdim'].sort_values()
        fdata = self.filter_data( sim=sim, tally=tally, nparticles=nparticles, xdim=xgrid.min(), ydim=xgrid.min(), nthreads=thread)
        for grid in xgrid.iloc[1:]:
            thread = 4 * thread
            fdata1 = self.filter_data( sim=sim, tally=tally, nparticles=nparticles, xdim=grid, ydim=grid, nthreads=thread)            
            fdata = pd.concat([fdata,fdata1],axis=0)
        return fdata

    # Add plots for weak scaling with grid size
    def plot_weak_scaling_grid_eiron(self,firstthread,nparticles):
        fdata = self.get_weak_scaling_grid('monolithic','private',firstthread,nparticles)
        self.plot_data(fdata['nthreads'],fdata['walltime'],{'label':'Monolithic Private'})
        fdata = self.get_weak_scaling_grid('monolithic','shared',firstthread,nparticles)
        self.plot_data(fdata['nthreads'],fdata['walltime'],{'label':'Monolithic Shared'})
        fdata = self.get_weak_scaling_grid('synchronous','private',firstthread,nparticles)
        self.plot_data(fdata['nthreads'],fdata['walltime'],{'label':'Synchronous Private'})
        fdata = self.get_weak_scaling_grid('synchronous','shared',firstthread,nparticles)
        self.plot_data(fdata['nthreads'],fdata['walltime'],{'label':'Synchronous Shared'})

    # Make figure of weak scaling with grid size
    def fig_weak_scaling_grid_eiron(self,firstthread,nparticles):
        print('Plotting weak scaling with grid size:', nparticles, ' particles' )
        self.init_figure('Weak scaling with grid size. '+ str(nparticles) + ' particles')
        self.plot_weak_scaling_grid(firstthread,nparticles)
        self.ax.set_yscale('linear')
        mingrid = self.pdata['xdim'].min()
        self.ax.annotate('Smallest grid=%4d' %(mingrid), xy=(1.7,3.4)) 
        self.ax.legend()
        self.makeplot('weakScalingGrid.png')

        
    # Add plots for weak scaling with grid size for pathstime and tallytime
    def plot_weak_scaling_split_grid_eiron(self,firstthread,nparticles):
        fdata = self.get_weak_scaling_grid('synchronous','private',firstthread,nparticles)
        self.plot_data(fdata['nthreads'],fdata['pathstime'],{'label':'pathstime private'})
        fdata = self.get_weak_scaling_grid('synchronous','shared',firstthread,nparticles)
        self.plot_data(fdata['nthreads'],fdata['pathstime'],{'label':'pathstime shared'})
        fdata = self.get_weak_scaling_grid('synchronous','private',firstthread,nparticles)
        self.plot_data(fdata['nthreads'],fdata['tallytime'],{'label':'tallytime private'})
        fdata = self.get_weak_scaling_grid('synchronous','shared',firstthread,nparticles)
        self.plot_data(fdata['nthreads'],fdata['tallytime'],{'label':'tallytime shared'})

    # Make figure of weak scaling with grid size for pathstime and tallytime
    def fig_weak_scaling_split_grid_eiron(self,firstthread,nparticles):
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
    def plot_speedup_eiron(self,nparticles,xgrid,ygrid):
        fdata = self.filter_data(sim='monolithic',tally='private', nparticles=nparticles, xdim=xgrid, ydim=ygrid)
        speedup = fdata['walltime'].iloc[0]/fdata['walltime']
        self.plot_data(fdata['nthreads'],fdata['nthreads'],{'label':'ideal'})
        self.plot_data(fdata['nthreads'],speedup,{'label':'monolithic private'})
        fdata = self.filter_data(sim='monolithic',tally='shared', nparticles=nparticles, xdim=xgrid, ydim=ygrid)
        speedup = fdata['walltime'].iloc[0]/fdata['walltime']
        self.plot_data(fdata['nthreads'],speedup,{'label':'monolithic shared'})
        fdata = self.filter_data(sim='synchronous',tally='private', nparticles=nparticles, xdim=xgrid, ydim=ygrid)
        speedup = fdata['walltime'].iloc[0]/fdata['walltime']
        self.plot_data(fdata['nthreads'],speedup,{'label':'synchronous private'})
        fdata = self.filter_data(sim='synchronous',tally='shared', nparticles=nparticles, xdim=xgrid, ydim=ygrid)
        speedup = fdata['walltime'].iloc[0]/fdata['walltime']
        self.plot_data(fdata['nthreads'],speedup,{'label':'synchronous shared'})

    # Make figure of speedup
    def fig_speedup_eiron(self,nparticles,xgrid,ygrid):
        print('Plotting speedup:', nparticles, 'particles')
        self.init_figure('OpenMP speedup: ' + str(nparticles) + ' particles, ' + str(xgrid) + 'x' + str(ygrid) + ' grid')
        self.plot_speedup_eiron(nparticles,xgrid,ygrid)
        self.ax.set_ylabel('Speedup')
        self.ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1f}"))
        self.ax.legend()
        self.makeplot('openMPSpeedup.png')

        
    # Add plots of efficiency
    def plot_efficiency_eirene(self,nparticles,xgrid,ygrid):
        fdata = self.filter_data(sim='monolithic',tally='private', nparticles=nparticles, xdim=xgrid, ydim=ygrid)
        efficiency = fdata['walltime'].iloc[0] / (fdata['walltime'] * fdata['nthreads'])
        self.plot_data(fdata['nthreads'],efficiency,{'label':'monolithic private'})
        fdata = self.filter_data(sim='monolithic',tally='shared', nparticles=nparticles, xdim=xgrid, ydim=ygrid)
        efficiency = fdata['walltime'].iloc[0] / (fdata['walltime'] * fdata['nthreads'])
        self.plot_data(fdata['nthreads'],efficiency,{'label':'monolithic shared'})
        fdata = self.filter_data(sim='synchronous',tally='private', nparticles=nparticles, xdim=xgrid, ydim=ygrid)
        efficiency = fdata['walltime'].iloc[0] / (fdata['walltime'] * fdata['nthreads'])
        self.plot_data(fdata['nthreads'],efficiency,{'label':'synchronous private'})
        fdata = self.filter_data(sim='synchronous',tally='shared', nparticles=nparticles, xdim=xgrid, ydim=ygrid)
        efficiency = fdata['walltime'].iloc[0] / (fdata['walltime'] * fdata['nthreads'])
        self.plot_data(fdata['nthreads'],efficiency,{'label':'synchronous shared'})

    # Add plots of efficiency
    def plot_efficiency_eiron(self,nparticles,xgrid,ygrid):
        fdata = self.filter_data(sim='monolithic',tally='private', nparticles=nparticles, xdim=xgrid, ydim=ygrid)
        efficiency = fdata['walltime'].iloc[0] / (fdata['walltime'] * fdata['nthreads'])
        self.plot_data(fdata['nthreads'],efficiency,{'label':'monolithic private'})
        fdata = self.filter_data(sim='monolithic',tally='shared', nparticles=nparticles, xdim=xgrid, ydim=ygrid)
        efficiency = fdata['walltime'].iloc[0] / (fdata['walltime'] * fdata['nthreads'])
        self.plot_data(fdata['nthreads'],efficiency,{'label':'monolithic shared'})
        fdata = self.filter_data(sim='synchronous',tally='private', nparticles=nparticles, xdim=xgrid, ydim=ygrid)
        efficiency = fdata['walltime'].iloc[0] / (fdata['walltime'] * fdata['nthreads'])
        self.plot_data(fdata['nthreads'],efficiency,{'label':'synchronous private'})
        fdata = self.filter_data(sim='synchronous',tally='shared', nparticles=nparticles, xdim=xgrid, ydim=ygrid)
        efficiency = fdata['walltime'].iloc[0] / (fdata['walltime'] * fdata['nthreads'])
        self.plot_data(fdata['nthreads'],efficiency,{'label':'synchronous shared'})

    # Make figure of efficiency
    def fig_efficiency_eiron(self,nparticles,xgrid,ygrid):
        print('Plotting efficiency:', nparticles, 'particles,', xgrid ,'x',ygrid,'grid')
        self.init_figure('OpenMP efficiency: ' + str(nparticles) + ' particles, ' + str(xgrid) + 'x' + str(ygrid) + ' grid')
        self.plot_efficiency(nparticles,xgrid,ygrid)
        self.ax.set_ylabel('Efficiency')
        self.ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1f}"))
        self.ax.legend()
        self.makeplot('openMPEfficiency.png')

######################################################################################
# main routines
######################################################################################

        
# Read command line arguments
def read_args():
    parser = argparse.ArgumentParser(description = 'Plot scalings from csv file containing scaling data.')
    parser.add_argument( "-f" ,"--file" , default="./SCALABILITY_REPORT/report_2D-D_slab.json" ,
                         help='csv file containing eiron profiling data (default=eiron_profile.csv)' )
    parser.add_argument( "-n" ,"--nparticles" , default=800000 ,
                         help='Number of particles for scaling plots' )
    parser.add_argument( "-g" ,"--gridsize" , default=256 ,
                         help='Grid size for scaling plots' )
    parser.add_argument( "--firstthread" , default=2 ,
                         help='First thread number for use in weak scaling plots' )
    parser.add_argument( "-t" ,"--filetype" , choices=['csv','json'],
                         help='File type conataining profile data, either csv or json. Overrides extension checking' )
    parser.add_argument( "-c" ,"--codetype" , default='eirene', choices=['eiron','eirene'],
                         help='Source of the profile data, either eiron or eirene.' )
    
    args = parser.parse_args()

    if not exists(args.file):
        print("Error:", args.file, "does not exist, exiting")
        exit()
                
    return args

# Read data and make the plots
def make_plots():
    print("Plotting scaling profiles")
    
    args = read_args()

    # Create object from from csv file
    prof = eiron_profile( args )

    #pd.options.display.max_rows = 999
    #pd.options.display.max_columns = 999
    #print(self.pdata)

    # Make some plots
    # Eirene
    if args.codetype == 'eirene':
        args.nparticles = 11000
        prof.omp_simple_strong_scaling_eirene(args.nparticles)
        prof.mpi_simple_strong_scaling_eirene(args.nparticles,1)
        prof.hybrid_strong_scaling_eirene(args.nparticles)
        prof.mpi_speedup_eirene(args.nparticles)
        prof.omp_speedup_eirene(args.nparticles)
        prof.mpi_efficiency_eirene(args.nparticles)
        prof.omp_efficiency_eirene(args.nparticles)
        prof.mpi_strong_scaling_split_eirene(args.nparticles,1)
        prof.omp_strong_scaling_split_eirene(args.nparticles)

    # Eiron
    if args.codetype == 'eiron':
        prof.fig_strong_scaling_eiron(args.nparticles,args.gridsize,args.gridsize)
        prof.fig_strong_scaling_split_eiron(args.nparticles,args.gridsize,args.gridsize)
        prof.fig_weak_scaling_particles_eiron(args.firstthread,args.gridsize,args.gridsize)
        prof.fig_weak_scaling_split_particles_eiron(args.firstthread,args.gridsize,args.gridsize)
        prof.fig_weak_scaling_grid_eiron(args.firstthread,args.nparticles)
        prof.fig_weak_scaling_split_grid_eiron(args.firstthread,args.nparticles)
        prof.fig_speedup_eiron(args.nparticles,args.gridsize,args.gridsize)
        prof.fig_efficiency_eiron(args.nparticles,args.gridsize,args.gridsize)

##################################################################################
# Do stuff from the cli
##################################################################################
if __name__ == '__main__':
    make_plots()



