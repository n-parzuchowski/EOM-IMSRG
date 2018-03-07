import sys
import os.path

email = "nathan.parz@gmail.com"

mem = ['500mb','500mb','500mb''500mb','1gb','2gb','3gb','4gb','5gb','6gb','20gb','20gb','40gb','40gb','50gb','50gb','60gb','60gb'] 
wtime = ['00:20:00','00:20:00','00:20:00', '00:20:00','00:40:00','01:00:00','02:00:00','02:00:00',\
          '04:00:00','04:00:00','10:00:00','10:00:00','24:00:00','24:00:00','24:00:00','24:00:00']
#wtime = [ '00:20:00','00:40:00','01:00:00','02:00:00','02:00:00',\
 #         '04:00:00','04:00:00','10:00:00','10:00:00','24:00:00','24:00:00','12:00:00','12:00:00']
#wtime = 20*["01:00:00"]
ompnum = ['8','8','8','8','8','8','8','8','8','8','8','8','8','8','8','8','8']


hwdefault = ["12","16","20","24","28","32","36","40"]
Rdefault = ["8","10","12","14"]
default_TBME = "chi2b_srg0625"

print 
print 
print "This script prepares job submission scripts and input files for the EOM-IMSRG code." 
print "There are defaults set throughout the script that are chosen if you press enter instead of entering values."
print "Entering walltimes, mem and ppn requirements is too painful, so those are hardcoded in the script. you can change them there."
print 
print 
GROUP = raw_input("group to charge resources to. Leave blank if not applicable: ")

nuc = raw_input('Enter nucleus name: (He4,O16,Ca40,etc...) ' ) 

if (nuc == 'H2'): 
    nprot = 0
    nneut = 2   # this is wrong but the progam only plays nice with closed shells.
    # don't worry, the deuteron module is separate and fixes this.
    print "Computing deuteron with FCI in lab-frame."
    nuc = "deut"
elif (nuc == 'He4'): 
    nprot = 2
    nneut = 2
elif (nuc == 'O16'): 
    nprot = 8
    nneut = 8
elif (nuc == 'O22'): 
    nprot = 8
    nneut = 14
elif (nuc == 'O24'): 
    nprot = 8
    nneut = 16
elif (nuc == 'Ca40'): 
    nprot = 20
    nneut = 20
elif (nuc == 'Ca48'): 
    nprot = 20
    nneut = 28
elif (nuc == 'Ni56'): 
    nprot = 28
    nneut = 28
elif ( nuc == 'C14'):
    nprot = 6
    nneut = 8
elif ( nuc == 'C12'):
    nprot = 6
    nneut = 6
elif ( nuc == 'Si28'):
    nprot = 14
    nneut = 14
else:
    print 'Invalid Entry'
    sys.exit()
    

hwstr = raw_input('\nEnter hw values, seperated by commas: ( 20,22,24 ) ') 
if (hwstr == ""):
    print "USING DEFAULT HW:", hwdefault
    print 
    hwlist = hwdefault
else:
    hwlist = hwstr.strip().split(',')

Rstr = raw_input('\nEnter eMax values, seperated by commas: (3,5,7,9) ' ) 

if (Rstr == ""):
    print "USING DEFAULT eMax:", Rdefault
    print 
    Rlist = Rdefault
else:
    Rlist = Rstr.strip().split(',') 
    

twobody_prefix = raw_input('\nEnter two-body interaction file prefix: (chi2b_srg0625) ') 
if (twobody_prefix==""):
    print "using default: chi2b_srg0625"
    twobody_prefix=default_TBME

srgloc = twobody_prefix.index("srg")+3

srg_param = twobody_prefix[srgloc:srgloc+4] 

Jramp_n2lo = "chi2b3b400cD-02cE0098_hwconv036_srg"+srg_param+"ho40J"
Cramp_n2lo = "chi2b3b_srg"+srg_param+"ho40C" 
n2loSAT= "chi2b3bSAT_J7666_hwconv022_JT3Nfull73_srg0000ho40C"
magic = "em1.8-2.0"
almostmagic = "em2.0-2.0"
print "\nOPTIONS FOR THREE-BODY FORCE:"
print "1: "+Jramp_n2lo
print "2: "+Cramp_n2lo
print "3: "+n2loSAT
print "4: "+magic
print "5: "+almostmagic
print "x: enter file prefix"
print "none:  no threebody"  

threebody_enter= raw_input("select and option: ")
if (threebody_enter=="1"):
    threebody_prefix = Jramp_n2lo
elif (threebody_enter=="2"):
    threebody_prefix = Cramp_n2lo
elif (threebody_enter=="3"):
    threebody_prefix = n2loSAT
elif (threebody_enter=="4"):
    threebody_prefix = magic
elif (threebody_enter=="5"):
    threebody_prefix = almostmagic
elif (threebody_enter.lower()=="x"):
    threebody_prefix = raw_input('\nenter a threebody file, or type "none": ')
elif (threebody_enter.lower()=="none"):
    threebody_prefix = "none"
else:
    threebody_prefix = "none"

    
lMax = raw_input('\nenter lMax: ')
if (lMax==""): lMax="10"
Rfile = raw_input('\nenter eMax for interaction file: ')
if (Rfile==""): Rfile="14"
lMaxfile = raw_input('\nenter lMax for interaction file: ')
if (lMaxfile==""): lMaxfile="10"

if threebody_prefix.lower()=="none": 
    E3Max = '0'
    E3Maxfile='0'
else:
    E3Max = raw_input('\nenter E3Max cutoff: ')
    E3Maxfile = raw_input('\nenter E3Max cutoff for interaction file: ')
    if (E3Max==""): E3Max="0"
    if (E3Maxfile==""): E3Maxfile="14"
 
hamtype = raw_input('\nFor intrinsic hamiltonian: "1" for harmonic trap: "2", full: "3": ') 
if (hamtype==""): hamtype="1"
hf = raw_input('\nFor HF: "HF". Otherwise type: "HO": ') 
if (hf==""): hf="HF"
mag = raw_input('\nFor magnus: "mag". Otherwise: "trad" or "disc": ' ) 
if (mag==""): mag="mag"
tda = raw_input('\nselect calculation: GS,EOM: ' )

if (tda.lower()=="gs"):
    tdaint = "0"
else:
    tdaint = "1"

hw_com = '0.0'
CMint = '0'
RRMSint = '0'
chkrestart=False
                 
print "\nTIP FOR COMPUTING OBSERVABLES"
print "If IMSRG calculation will not finish in one submission (checkpointing required), don't run observables."
print "Regenerate submission scripts to evolve observables after the IMSRG calculation finishes, with Omega written to file." 
print "THIS ONLY WORKS WITH MAGNUS. Be sure to write/read OMEGA (see option below)"

print '\nTo calculate c.m. expectation value, and factorization frequency'
CMint = raw_input('enter 1, otherwise enter 0: ') 
if (CMint==""): CMint="0"

print '\nTo calculate Point Nucleon RMS radius,' 
RRMSint = raw_input('enter 1, otherwise enter 0: ')         
if (RRMSint==""): RRMSint="0"

other= raw_input('\nenter other scalar observables (rho21 or none): ')
other = other.lower()

if other =="":
    other = "none"

print "\n.eom file specifies excited states and tensor observables" 
eomfile = raw_input("enter .eom file name, if applicable: ") 
if (eomfile==""): eomfile="standard.eom"

lawbeta = raw_input("\nenter lawson beta value: ")
if lawbeta == "":
    lawbeta = "0.0"

com_freq = raw_input("\nenter hw for Hcm: ")
if com_freq == "":
    com_freq = "0.0"
        

omega_string = '.false. , .false.'
bare_string = '.false. , .false.'
dc_string = '.false. , .false.'
checkpoint_str = '.false.'
write_human_str = '.false.'
try:
    print '\nDo you want to write/read some operators to/from file?' 
    writing_stuff = raw_input('\n(Y/N)?') 

    if (writing_stuff.lower() == 'y'):
         print 'answer "y" or "n":'  

         writehuman = raw_input('\nwrite to human readable form?: ')

         if writehuman.lower() == 'y':
             write_human_str = '.true.'

         writebare = raw_input('\nwrite bare operators: ')   
         readbare = raw_input('\nread bare operators: ')
   
         if writebare.lower() == 'y':
             
             if readbare.lower() == 'y':
                 bare_string = '.true. , .true.' 
             elif readbare.lower()=='n': 
                 bare_string = '.true. , .false.'
             else: 
                 print 'Try again.' 
                 raise IndexError 
         elif writebare.lower()=='n':
             if readbare.lower() == 'y':
                 bare_string = '.false. , .true.' 
             elif readbare.lower()=='n': 
                 bare_string = '.false. , .false.'
             else: 
                 print 'Try again.' 
                 raise IndexError
         else: 
             print 'Try again.' 
             raise IndexError
 
         writedc = raw_input('\nwrite decoupled operators: ')   
         readdc = raw_input('\nread decoupled operators: ')   

         if writedc.lower() == 'y':
             
             if readdc.lower() == 'y':
                 dc_string = '.true. , .true.' 
             elif readdc.lower()=='n': 
                 dc_string = '.true. , .false.'
             else: 
                 print 'Try again.' 
                 raise IndexError 
         elif writedc.lower()=='n':
             if readdc.lower() == 'y':
                 dc_string = '.false. , .true.' 
             elif readdc.lower()=='n': 
                 dc_string = '.false. , .false.'
             else: 
                 print 'Try again.' 
                 raise IndexError
         else: 
             print 'Try again.' 
             raise IndexError

         if (mag.lower() == 'mag'):
             writeomega = raw_input('\nwrite omega: ')   
             readomega = raw_input('\nread omega: ')   
             
             if writeomega.lower() == 'y':
             
                 if readomega.lower() == 'y':
                     omega_string = '.true. , .true.' 
                 elif readomega.lower()=='n': 
                     omega_string = '.true. , .false.'
                 else: 
                     print 'Try again.' 
                     raise IndexError 
             elif writeomega.lower()=='n':
                 if readomega.lower() == 'y':
                     omega_string = '.false. , .true.' 
                 elif readomega.lower()=='n': 
                     omega_string = '.false. , .false.'
                 else: 
                     print 'Try again.' 
                     raise IndexError
             else: 
                 print 'Try again.' 
                 raise IndexError
             checkpointing = raw_input('\ncheckpointing: ') 
             
             if checkpointing.lower() == 'y':
                 checkpoint_str = '.true.'
                 restart = raw_input('\n checkpoint restart? (y/n): ')
                 if (restart.lower()=='y'):
                     chkrestart=True
             elif checkpointing.lower() == 'n':
                 checkpoint_str = '.false.'
             else:
                 print 'Try again.' 
                 raise IndexError
             

             
except IndexError:
    print 'something else is wrong...' 
                  

    
fq = open('run_all.bat','w')

fq.write( '#!/bin/bash \n\n') 

if hf == 'HF':
    HFint = '1'
else:
    HFint = '2'
    
if mag == 'mag':
    magint = '1'
    quads = raw_input('\nDo you want to correct for quadrupoles? (Y/N) ')
    if (quads.lower() == 'y'):
        trips = raw_input('\nDo you want to correct for triples? (Y/N) ')
        if (trips.lower() == 'y'):
            magint='5'
            mag = mag+'TQ'
        else:
            magint='4'
            mag = mag+'Q'
    else:
        trips = raw_input('\nDo you want to correct for triples? (Y/N) ')
        if (trips.lower() == 'y'):
            magint='6'
            mag = mag+'T'

elif mag == 'disc':
    magint = '3'
else:
    magint = '2'

for R in Rlist:
    for hw in hwlist:
                    
                    hwx = int(hw) 
                    Rx = int(R) 

                    if magint=='5':
                        if Rx > 100: 
                            resubmit=True
                        else:
                            resubmit=False
                    else:
                        if Rx > 100: 
                            resubmit=True
                        else:
                            resubmit=False

                    memreq = mem[Rx] 
                    timreq = wtime[Rx]
                    threebody_file="none"
                    jobname = nuc+"_"+twobody_prefix+'_eMax'+(2-len(R))*'0'+R+'_hwHO0'+hw 
                    if (int(lMax) < int(R) ) :
                        jobname = nuc+"_"+twobody_prefix+'_eMax'+(2-len(R))*'0'+R+'_lMax'+(2-len(lMax))*'0'+lMax+'_hwHO0'+hw
                    TBMEfile = twobody_prefix+'_eMax'+(2-len(Rfile))*'0'+Rfile+'_hwHO0'+hw+'.me2j.gz'
                    spfile = 'hk'+Rfile+'.sps'
                    if (threebody_prefix != "none"):
                        threebody_file = threebody_prefix+'_eMax'+(2-len(Rfile))*'0'+Rfile+'_EMax'+(2-len(E3Maxfile))*'0'+E3Maxfile+'_hwHO0'+hw+'.me3j.bin'  
                        jobname=jobname.replace('2b_','2b3b_')
                    if (int(lMaxfile) < int(Rfile) ) :                    
                        TBMEfile = twobody_prefix+'_eMax'+(2-len(Rfile))*'0'+Rfile+'_lMax'+(2-len(lMaxfile))*'0'+lMaxfile+'_hwHO0'+hw+'.me2j.gz'
                        spfile = 'hk'+Rfile+"_lMax"+lMaxfile+'.sps'
                        if (threebody_prefix != "none"):
                            threebody_file = threebody_prefix+'_eMax'+(2-len(Rfile))*'0'+Rfile+\
                                             '_lMax'+(2-len(lMaxfile))*'0'+lMaxfile+'_EMax'+(2-len(E3Maxfile))*'0'+E3Maxfile+'_hwHO0'+hw+'.me3j.bin'  
                            jobname=jobname.replace('2b_','2b3b_')

                    initfile = jobname+'.ini'
                    # write pbs file ===========================        
                    fx = open('pbs_'+jobname,'w') 

                    fx.write('#!/bin/sh \n\n')
                    fx.write('#PBS -l walltime='+timreq+'\n')
                    fx.write('#PBS -l nodes=1:ppn='+ompnum[Rx]+'\n')
                    fx.write('#PBS -l mem='+memreq+'\n') 
                    fx.write('#PBS -j oe\n')
                    fx.write('#PBS -N '+jobname+'\n') 
                    fx.write('#PBS -M '+email+'\n')

                    if GROUP =="":                        
                        fx.write('#PBS -m a\n\n')
                    else:
                        fx.write('#PBS -m a\n')
                        fx.write('#PBS -A '+GROUP+'\n\n')

                    fx.write('cd $HOME/nuclear_IMSRG/src/im-SRG\n\n')
                    fx.write('export OMP_NUM_THREADS='+ompnum[Rx]+'\n\n')

                    if (chkrestart):
                        fx.write('./run_IMSRG '+initfile+' pbs_'+jobname+'\n')
                    else:
                        fx.write('./run_IMSRG '+initfile+'\n')

                    fx.write('qstat -f ${PBS_JOBID}\nexit 0\n')
                    fx.close()

                    fq.write('qsub pbs_'+jobname+'\n') 

                    # write ini file ==========================



                    if (magint=='4'):
                        magout = "1,'y','n'"
                    elif (magint=='5'):
                        magout = "1,'y','y'"
                    elif (magint=='6'):
                        magout = "1,'y','y'"
                    else:
                        magout = magint+",'n','n'"

                    fx = open('inifiles/'+initfile,'w') 

                    fx.write('##########################\n')
                    fx.write('#### IMSRG INPUT FILE ####\n')
                    fx.write('#### KEEP THIS FORMAT ####\n')
                    fx.write('##########################\n')
                    fx.write('##########################\n')
                    fx.write('# ENTER OUTPUT FILE PREFIX\n')
                    fx.write(jobname +'\n') 
                    fx.write('# ENTER INTERACTION FILE NAME\n')
                    fx.write(TBMEfile + '\n')
                    fx.write('# ENTER eMax, lMax\n')
                    fx.write(R+','+lMax+'\n')
                    fx.write('# ENTER 3B INTERACTION FILE NAME\n')
                    fx.write(threebody_file + '\n') 
                    fx.write('# ENTER E3Max (enter "0" for no three body force)\n')
                    fx.write(E3Max + '\n') 
                    fx.write('# ENTER SINGLE PARTICLE INPUT FILE NAME\n')
                    fx.write(spfile + '\n')
                    fx.write('# ENTER HAMILTONIAN TYPE\n')
                    fx.write('# 1: T-V - Tcm -Vcm  2: harmonic trap T+U+V  3. T+V\n')
                    fx.write(hamtype+'\n') 
                    fx.write('# ENTER HO SPACING hw\n')
                    fx.write(hw +'\n')
                    fx.write('# ENTER NUMBER OF PROTONS\n')
                    fx.write(str(nprot)+'\n')
                    fx.write('# ENTER NUMBER OF NEUTRONS\n')
                    fx.write(str(nneut)+'\n')
                    fx.write('# ENTER 1 for HF basis\n')
                    fx.write('# OR 2 for HO basis\n')
                    fx.write(HFint+'\n')
                    fx.write('# ENTER 1 for magnus method\n')
                    fx.write('# or 2 for traditional ode (y/n): "quads", "trips"\n') 
                    fx.write(magout+'\n') 
                    fx.write('# 0: gs only, 1: EOM, 2: TDA\n')
                    fx.write('# ENTER 0 for ground state only\n') 
                    fx.write(tdaint+'\n')
                    fx.write('# ENTER 1 TO CALCULATE Hcm, 0 otherwise\n')
                    fx.write(CMint+'\n')
                    fx.write('# ENTER 1 TO CALCULATE Rrms, 0 otherwise\n') 
                    fx.write(RRMSint+'\n') 
                    fx.write('# ENTER other observables or "none"\n') 
                    fx.write(other+'\n')
                    fx.write('# EOM file (standard.eom)\n')
                    fx.write(eomfile+'\n')            
                    fx.write('# Lawson beta value, c.m. hw\n') 
                    fx.write(lawbeta+','+com_freq+'\n') 
                    fx.write('# .true. for checkpointing (only applicable to magnus)\n' )
                    fx.write(checkpoint_str+'\n')
                    fx.write('# write normal ordered bare, read normal ordered bare\n')
                    fx.write(bare_string+'\n')                                     
                    fx.write('#  write normal ordered decoupled, read normal ordered decoupled\n')
                    fx.write(dc_string+'\n')
                    fx.write('#  write omega, read omega\n')
                    fx.write(omega_string+'\n')
                    fx.write('#  write human readable\n')
                    fx.write(write_human_str+'\n')
                    fx.write('########################################################\n')
                    fx.write('# NOTES \n')
                    fx.write('#\n')
                    fx.write("# 1. THIS FILE'S NAME SHOULD BE SUPPLIED AS THE ONLY\n") 
                    fx.write('# COMMAND ARGUMENT WITH THE EXECUTABLE "./run_IMSRG"\n')
                    fx.write('#\n')
                    fx.write('# 2. THIS FILE IS READ BY SUBROUTINE "read_main_input_file"\n')
                    fx.write('# which is found in "basic_IMSRG.f90" \n')
                    fx.write('# \n')
                    fx.write('# 3. FILENAMES ARE READ TO COMMON BLOCK "files" \n')
                    fx.write('########################################################\n')

                    fx.close()
                    
        
        
        
fq.close()

os.system("chmod 0755 run_all.bat")
