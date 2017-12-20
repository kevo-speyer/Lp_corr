#!/bin/python
# It's necessary a dynamic lybrary to compute the correlations:
# $ f2py -c cos_corr.f90 --fcompiler=intelem -m fast_corr

def main(): # main_script
    """
    This program reads film_xmol and calculates the orientational correlation between bonds
    """
    import numpy as np

    # Read simulatino parameters
    Lx, Ly, Lz, frame_tot, n_brush, n_mon_d, n_chain_d, dt =  read_params() # read system paramets 
    
    # Read data and do data analysis
    #R_abs = get_norm_R_next(n_liq, n_mon, n_chain, frame_tot, Lx, Ly)
    corr, a = get_corr_Rnext(n_brush, n_mon_d, n_chain_d, frame_tot, Lx, Ly, Lz)
    
    # Estimate persistence length (l_p), from correlations of r_next vectors
    l_p = -a/np.log(corr[1]) # in length units

    # Output data
    save_data(corr,'mean_Rnext_corr.mide')
    #save_data(R_abs,'mean_abs_Rnext.mide')
    save_data(np.array(l_p),'mean_lp.mide')

    print "\n Persistence length estimation successfully calculated and saved in file 'mean_lp.mide' "
    print " Result is in distance units (i.e. sigma for Lennard-Jonnes Potentials) \n"


################### END MAIN PROGRAM ######################

def get_corr_Rnext(n_brush, n_mon_d, n_chain_d, frame_tot, Lx, Ly, Lz):
    """ calculate cos correlations"""
    import numpy as np 
    from fast_corr import cos_corr

    corr = np.zeros(n_mon_d-1)
    corr_tmp = np.zeros_like(corr)
    corr_tmp_py = np.zeros_like(corr)
    r0_d = np.empty([3,n_mon_d*n_chain_d])
    R_next = np.empty([3,n_mon_d-1])
    R_next_new = np.empty([3,n_mon_d-1])
    
    a = 0.
    j = 0
    for i_frame in range( 1,int(frame_tot+1) ): #  2): 
        r0_d = read_chain_pos(i_frame, n_brush, int(n_mon_d*n_chain_d) ) # Read Brush positions

        for i_chain in range(n_chain_d):             
            # Get R_next vectors from r0
            mon_init = i_chain * n_mon_d
            mon_fin = (i_chain+1) * n_mon_d
            #print r0_d[:,mon_init:mon_fin]
            R_next = r0_d[:,1+mon_init:mon_fin] - r0_d[:,mon_init:mon_fin-1]             
           
            # CORRECT FOR PBC    
            R_next_new = PBC_corr(R_next, Lx, Ly, Lz)
            #print R_next_new[:,0],R_next_new[:,2]
            
            #Get mean distance between monomers
            a_tmp = np.linalg.norm(R_next_new,axis=0)
            a += np.mean(a_tmp)    

            # Calculate auto Correlation for each chain          
            corr_tmp = cos_corr( R_next_new ) #, n_mon_d-1, 3 ) # Don't pass optional arguments
            #corr_tmp = r_next_c_py(R_next_new, n_mon_d-1) # python routine, to calculate the same as above
            # Python routin is 300 x slower than fortran         
            
            corr += corr_tmp
            j += 1 
    
    # Normalize corr and print results        
    corr = corr / float(j)
    a = a / float(j)
  
    return corr, a 

def get_var_Rnext(R_mean, n_liq, n_mon, n_chain, frame_tot, Lx, Ly): 
    """ Get Rnext Variance"""
    import numpy as np

    R_var = np.zeros(n_mon-1)  
    j=0
    for i_frame in range( 1, int(frame_tot+1) ):
        r0_b = read_chain_pos(i_frame, n_liq, int(n_mon*n_chain) ) # Read Brush positions

        for i_chain in range(n_chain):             
            # Get R_next vectors from r0
            mon_init = i_chain * n_mon
            mon_fin = (i_chain+1)*n_mon

            R_next = r0_b[:,1+mon_init:mon_fin] - r0_b[:,mon_init:mon_fin-1] 
            
            # INVERT CHAINS FOT TOP WALL 
            if R_next[2,0] < 0:
                R_next = -R_next

            # CORRECT FOR PBC    
            R_next = PBC_corr(R_next, Lx, Ly)
 
            R_var += ((R_mean - R_next)**2).sum(0)
            j += 1

    R_var = R_var / float(j)
    
    return R_var

def get_mean_Rnext(n_liq, n_mon, n_chain, frame_tot, Lx, Ly):
    """Calculate mean Rnext"""
    import numpy as np
    
    R_mean = np.zeros((3,n_mon-1))
    j = 0

    for i_frame in range( 1, int(frame_tot+1) ):
        r0_b = read_chain_pos(i_frame, n_liq, int(n_mon*n_chain) ) # Read Brush positions

        for i_chain in range(n_chain):             
            
            # Get R_next vectors from r0
            mon_init = i_chain * n_mon
            mon_fin = (i_chain+1)*n_mon

            R_next = r0_b[:,1+mon_init:mon_fin] - r0_b[:,mon_init:mon_fin-1] 
            
            # INVERT CHAINS FOT TOP WALL 
            if R_next[2,0] < 0:
                R_next = -R_next
           
            # CORRECT FOR PBC    
            R_next = PBC_corr(R_next, Lx, Ly)
            
            if R_next.max() > 2.:
                print "ERROR" 
                print i_frame,i_chain,R_next.argmax(),R_next.max()
                print "i_frame,i_chain,R_next.argmax(),R_next.max()"

            # INTENTO VER PORQUE DA MAL
            R_mean += R_next
            j += 1

    R_mean = R_mean / float(j)
    
    return R_mean 

def save_data(data, string):
    """" Save histogram to file"""
    import numpy as np

    f = open(string, 'w+')

    if np.size(data.shape) == 0:
        s = str(data)
        s += str('\n')
        f.write(s)
              
    elif np.size(data.shape) == 1: 
       for i in range(data.shape[0]):  
           s = str(str(i)+' '+ str(data[i]) + '\n')
           f.write(s)

    elif np.size(data.shape) == 2:
        for i in range(data.shape[1]):
            s = str(i)
            for j in range(data.shape[0]):
                s += str(' ' + str(data[j,i]))

            s += str('\n')
            f.write(s)

  
    f.close()
 
def PBC_corr(R_next, Lx, Ly, Lz):
    """" Correct R_next for PBC in x, y and z"""
    import numpy as np

    R_PBC = np.zeros_like(R_next)

    R_PBC[0,:] = R_next[0,:]-Lx*(R_next[0,:]*2/Lx).astype(int)
    R_PBC[1,:] = R_next[1,:]-Ly*(R_next[1,:]*2/Ly).astype(int)
    R_PBC[2,:] = R_next[2,:]-Lz*(R_next[2,:]*2/Lz).astype(int) 
   
    #FOR DEBUGGING PURPOSES
    #print R_PBC
    a = np.linalg.norm(R_PBC,axis=0)
    for norma in a:
        if (norma>1.5):
            print "ERROR, stretched bond:", norma
    
    return R_PBC

def autocorr(data):
    """Calculate autocorreltacion function, given a 1D numpy array"""
    import numpy as np
    
    yunbiased = data-np.mean(data)
    ynorm = np.sum(yunbiased**2)
    acor = np.correlate(yunbiased, yunbiased, "same")/ynorm
    acor = acor[len(acor)/2:]
    
    return acor

def read_params():
    """
    This script reads Lx, Ly, Lz, n_frames, nm, n_brush, delta_t from conf_old and mfa_input
    """
    import numpy as np
    import os.path
    import sys

    #Check if necessary files are present
    for fl in ["mfa_input","conf_old","film_xmol"]:
        if not os.path.isfile(fl):
            print "ERROR: File '",fl,"' missing, aborting run!"
            sys.exit()

    j=0
    with open('mfa_input') as f:
        for line in f:
            j += 1
            if(j==3):
                n_steps = float(line.split()[0])
            if(j==4):
                n_obs = float(line.split()[0])    
            if(j==5):
                dt = float(line.split()[0])                                     
            if(j==12):
                Lx = float(line.split()[0])
            if(j==13):
                Ly = float(line.split()[0])
            if(j==14):
                Lz = float(line.split()[0])
    delta_t = dt * float(n_obs)    
    n_frames = int(n_steps / n_obs)                      

    #Now read system_input
    j=0
    with open('conf_old') as f:
        for line in f:
            j += 1
            if(j==3):
                Lx = Lx * float(line.split()[0])
                Ly = Ly * float(line.split()[1])

            if(j==2):    
                n_mon = float(line.split()[0])
                n_chain = float(line.split()[1])
                n_mon_d = float(line.split()[2])
                n_chain_d = float(line.split()[3])

            if(j>3):
                break
                    
    return Lx,Ly,Lz,int(n_frames),int(n_chain*n_mon),int(n_mon_d),int(n_chain_d), delta_t

def read_chain_pos(frame_nr, n_brush, n_liq):
    """This routine reads the positions of the melt and brush for a specific frame from file 'film_xmol'.
    The outputs are two numpy arrays x0 (melt) and x0_b(brush) with the dimension in the first index, and particle number
    in the second x0[0][0] is the x position of the first particle.
    This function should be calles as: r0 = read_conf(1); 
    to obtain the first frame (frame 0 does not exist)"""
    import numpy as np
    import linecache
    x0_b = np.empty([3,n_liq])

    frame_lines = int( n_liq + n_brush + 2 ) # Select lines to read from film_xmol
    line_init_b = int( (frame_nr - 1) * frame_lines + 3 )
    line_end_b  = int( line_init_b + n_brush)
    line_init_m = int( line_end_b )
    line_end_m  = int( line_init_m + n_liq )
    j=0

    #DEBUG
    #print line_init_m,line_end_m
    for i in range(line_init_m,line_end_m): # Save brush positions in one array
        #print i
        x = linecache.getline('film_xmol', i).split()[1:4]
        x = np.array(x, dtype='|S16')
        x = x.astype(np.float)     
        x0_b[:,j] = x[:]
        j+=1
    
    return x0_b # , v0, v0_b, r_cm # return are two  numpy arrays with melt particles in 'em   

def r_next_c_py(R_next, n_mon_1):
    """ Python routine to calculate orientational correlations in a polymer. It is very Slow, better 
    use Fortran"""
    import numpy as np
  
    corr_tmp = np.zeros(n_mon_1)
    corr = np.zeros(n_mon_1)
    n_event = np.zeros(n_mon_1)

    corr_tmp[:] = 0.
    corr[:] = 0.

    for i in range(n_mon_1):
        for j in range(i,n_mon_1):

            n_event[j-i] += 1.
            dmy = 0.
            dmy2 = 0.
            dmy3 = 0.

            for i_dim in range(3):
                dmy += R_next[i_dim,i] * R_next[i_dim,j] 
                dmy2 += np.square(R_next[i_dim,i])
                dmy3 += np.square(R_next[i_dim,j])
    
            #corr_tmp[j-i] += float(dmy) / float(np.sqrt(dmy2 * dmy3))

            corr_tmp[j-i] += np.dot( R_next[:,i], R_next[:,j]) /( np.linalg.norm(R_next[:,i])* np.linalg.norm(R_next[:,j])) 

    for i in range(n_mon_1):
        corr[i] = corr_tmp[i] / float(n_event[i]) #.astype(int)        

    return corr

def get_norm_R_next(n_brush, n_mon_d, n_chain_d, frame_tot, Lx, Ly):
    """ Calculate mean distance between consequent monomers"""
    import numpy as np
    
    R_mean = np.zeros(n_mon-1)
    j = 0

    for i_frame in range( 1, int(frame_tot+1) ):
        r0_b = read_chain_pos(i_frame, n_liq, int(n_mon*n_chain) ) # Read Brush positions

        for i_chain in range(n_chain):             
            
            # Get R_next vectors from r0
            mon_init = i_chain * n_mon
            mon_fin = (i_chain+1)*n_mon

            R_next = r0_b[:,1+mon_init:mon_fin] - r0_b[:,mon_init:mon_fin-1] 
            
            # INVERT CHAINS FOT TOP WALL 
            if R_next[2,0] < 0:
                R_next = -R_next
           
            # CORRECT FOR PBC    
            R_next = PBC_corr(R_next, Lx, Ly)
            
            R_mean += np.sqrt((R_next**2).sum(0))            
            j += 1

    R_mean = R_mean / float(j)
    
    return R_mean 

 
########### END SUPPORT ROUTINES ##########

if __name__ == "__main__":
    main()

