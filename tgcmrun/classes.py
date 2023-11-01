import re,os,subprocess,sys, stat
from interpic import interpic
import fnmatch

def find_file(pattern, path):
    """
    Find a file in the specified path that matches the given pattern. Assumes only one match.

    :param pattern: Pattern to look for in the file names.
    :param path: Path of the directory to search in.
    :return: File path if a match is found, else None.
    """
    for root, dirs, files in os.walk(path):  # Recursively go through all directories and subdirectories
        for name in files:
            if fnmatch.fnmatch(name, pattern):  # Check if file name matches the pattern
                return os.path.join(root, name)  # If so, return the file path immediately

    return None
#----------------------- Begin Class Version definition ---------------------
class Version:
  def name(self,model_name):
    if model_name == 'tiegcm':
#     return 'tiegcm1.95'
      return 'tiegcm2.0'
    elif model_name == 'timegcm':
      return 'timegcm1.42'
#
# Can specify a tag name for output history files, 
# (e.g. if tag=2.0, history file names would start with tiegcm2.0)
# The default is the null string (file names will start with tiegcm)
#
  def tag(self,model_name):
    if model_name == 'tiegcm':
      return ''      # can be null (e.g., testing or trunk version)
#     return '1.95'  # use previous tag to make startup history files
#     return '2.0'   
    elif model_name == 'timegcm':
      return '1.42'

#----------------------- Begin Class Namelist definition ---------------------
class Namelist:
 
  def make_namelist(self,default_file,new_file,list_mods,list_rm,run_fullname,tgcmdata):
#
# Make new namelist file with changes necessary for current run.
#
    default_list = self.file_to_list(default_file)                 # make 2d list from default namelist file
    newlist = self.modify_namelist(default_list,list_mods,list_rm) # modify list for this run
    self.list_to_file(newlist,new_file,run_fullname,tgcmdata)      # write modified list to new namelist file

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  def file_to_list(self,file):
#
# Create and return a list made from the provided namelist file.
#
# Cannot use these on read-only default namelist files:
#   replace_string(file,';','!')  # replace ';' with '!' (old namelist files)
#   replace_string(file,',',' ')  # replace ',' with ' '
#   remove_comments(file)         # remove commented text (after '!')

    f = open(file,'r')
    lines = f.readlines()
    lines = [line.rstrip() for line in lines]
    f.close()

    groupname = re.compile('^&\w+')
    endgroup = re.compile('^\/')
#
# Each line from the default namelist file will be a two-part list
# [keyword,value]. Multiple values will be concatenated to a single string.
#
    namelist = [] # new 2d list of [keyword = value]
    n = 0
    for line in lines:
      list = line.split('=')
      listlen = len(list)
      if listlen >= 2: 
        keyword = list[0]
        value = ''.join(list[1:])
        namelist.append([keyword,value])
        n = n+1
      elif listlen == 1:
        if groupname.match(list[0]):     # is group name (&tgcm_input)
          namelist.append([list[0],'']) # null value (could use "None"?)
          n = n+1
        elif endgroup.match(list[0]):
          continue
        else:                           # is continuation of values from previous line
          keyword = namelist[n-1][0]
          values = namelist[n-1][1]
          string = ''.join(values) + list[0]
          string = string.replace(" ","")
          string = string.replace("''","','")
          namelist.pop(n-1)
          namelist.append([keyword,string])

    for item in namelist:
      item[0] = item[0].strip()
      item[1] = item[1].strip()

    return namelist 
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  def modify_namelist(self,default_list,list_mods,list_rm):
#
# Modify values of default_list keywords found in list_mods.
# If keyword is not found in list_mods, use the default value.
# Each item in the lists is a 2-element list [keyword,value] pair.
#
    if len(list_mods) == 0: return default_list

    newlist = []
    for def_pair in default_list:
      keyword = def_pair[0]
      value = def_pair[1]
      for mods_pair in list_mods:
        if mods_pair[0] == def_pair[0]:
          value = mods_pair[1]
          break
      newlist.append([keyword,value])
#
# If any keywords in list_mods are not found in default_list, add them 
# to newlist with value from list_mods:
#
    for mods_pair in list_mods:
      found = 0
      for def_pair in default_list:
        if def_pair[0] == mods_pair[0]: 
          found = 1
          break
      if not found:
        keyword = mods_pair[0]
        value = mods_pair[1]
        newlist.append([keyword,value])
#
# Remove any keyword=value pairs from newlist as necessary:
# (pairs do not need values for this to work)
#
    for rm_pair in list_rm:
      n = 0
      for new_pair in newlist:
        if new_pair[0] == rm_pair[0]: 
          newlist.pop(n)
        n = n+1
    
    return newlist

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  def list_to_file(self,list,file,fullname,tgcmdata):
#
# Write a namelist file made from the provided 2d list.
#
    f = open(file,'w')
    groupname = re.compile('^&\w+')
    for pair in list:
      if groupname.match(pair[0]):
        f.write(pair[0]+'\n')
        comment = ( \
        "!\n",
        "! Namelist input file for run "+fullname+"\n"
        "!\n")
        for line in comment: f.write(line)
      else:
        value = pair[1]
        if '$TGCMDATA' in value:
          pair[1] = value.replace('$TGCMDATA',tgcmdata)
        line = pair[0] + ' = ' + pair[1]
        f.write(line+'\n')

    f.write('/\n')
    f.close()
#   print('Wrote namelist file ',file

#----------------------- Begin Class Job definition ---------------------
class Job:
  user = ''                # user name
  model_name = ''          # tiegcm or timegcm
  model_root = ''          # directory path (must exist)
  model_res = ''           # model resolution ('5.0' or '2.5')

  horires = ''
  vertres = ''
  zitop = ''
  mres = ''
  coupling = ''

  model_version = ''       # model version name
  machine = ''             # ch or linux
  script_default = ''      # job default script: file path (read)
  script_run = ''          # job script for run: file path (write)
  execdir = ''             # directory path (may or may not exist)
  bindir = ''             # directory path (may or may not exist)
  execute = ''             # execute model only if execute=TRUE
  stdin = ''               # namelist input file (read)
  stdout = ''              # standard out file (write)
  stdout_dir = ''          # directory to receive stdout files
  nprocs = ''              # number of processors (mpi tasks) (command line arg only)
  wallclock = ''           # wallclock limit for 5.0 or 2.5 deg res (ch only)
  project = ''             # authorized project number, e.g.: #PBS -P P28100036
  queue = ''               # LSF queue, e.g.: #PBS -q regular
  submitflag = ''          # Submit flag (not in job script)
  compiler = ''            # Compiler (intel, pgi, or gfort) (Linux platform only, not ch super)

  def make_jobscript(self,run):
    print('Openning job script',self.script_default)
    f = open(self.script_run,'w')
    found_input = 0
    found_output = 0
    found_execdir = 0
    found_bindir = 0
    found_runscript = 0
    for line in open(self.script_default):
#
# Set required parameters in job script (these have been set by default,
# or were received from the user by env var or prompts): 
#
# Model root directory:
      if 'set modeldir' in line:
        newline = 'set modeldir = '+self.model_root+'\n'
        f.write(newline)
#
# Execution directory:
      elif 'set execdir' in line and not found_execdir: # only first occurrence
        newline = 'set execdir = '+self.execdir+'\n'
        f.write(newline)
        found_execdir = 1
# Executable directory:
      elif 'set bindir' in line and not found_bindir: # only first occurrence      
        newline = 'set bindir = '+self.bindir+'\n'
        f.write(newline)
        found_bindir = 1        
#
# Namelist input file:
      elif 'set input' in line and not found_input: # only first occurrence
        newline = 'set input = '+self.stdin+'\n'
        f.write(newline)
        found_input = 1
#
# Stdout file:
      elif 'set output' in line and not found_output: # only first occurrence
        newline = 'set output = '+self.stdout+'\n'
        f.write(newline)
        found_output = 1
#
# Model resolution:
      elif 'set modelres' in line:
        newline = 'set modelres = '+self.model_res+'\n'
        f.write(newline)
#
# Horizontal resolution:
      elif 'set horires' in line:
        newline = 'set horires = '+self.horires+'\n'
        f.write(newline)
#
# Model resolution:
      elif 'set vertres' in line:
        newline = 'set vertres = '+self.vertres+'\n'
        f.write(newline)
#
# Model resolution:
      elif 'set zitop' in line:
        newline = 'set zitop = '+self.zitop+'\n'
        f.write(newline)
#
# Model resolution:
      elif 'set mres' in line:
        newline = 'set mres = '+self.mres+'\n'
        f.write(newline)
#
# Set some other options:
#
# Number of processors (mpi tasks):
# This depends on the default job script string "#PBS -l select" 
#   referring to nprocs for res2.5
# 
      elif '#PBS -l select' in line and '##PBS -l select' not in line: 
        print('Inside PBS select', line)
        if self.nprocs:
          newline = '#PBS -l select '+str(self.nprocs)+'\n'
          f.write(newline)
        else:
          if self.model_res == '5.0':
            newline = "#PBS -l select=2:ncpus=36:mpiprocs=36 \n"
            f.write(newline)
          else:
            f.write(line)
#
# Job name (pbs job only, e.g., ys):
      elif '#PBS -N' in line and not '##PBS -N' in line:
        newline = '#PBS -N '+run.fullname+'\n'
        f.write(newline)
#
# pbs script name (will be ignored if not ch):
      elif 'set runscript' in line and not found_runscript: # only first occurrence
        newline = 'set runscript = '+run.fullname+'.pbs\n'
        f.write(newline)
        found_runscript = 1
#
# Set wallclock limit e.g., for 30 minutes: #PBS -l walltime=00:30:00
# This will change for different resolutions.
# This will also change for different nprocs, but I will let that be
#   the responsibility of the user when setting command line options.
# Wallclock limit may be optionally set only by command line option.
# Default values are provided in set_run() method.
#
      elif '#PBS -l walltime' in line and '##PBS -l walltime' not in line:
        if self.wallclock:
          newline = '#PBS -l walltime='+self.wallclock
        f.write(newline+'\n')
#
# Project number (ys only) (command-line only), e.g.: #PBS -A P28100036
# Avoid changing any commented pbs commands (double '##')
# This will be validated by LSF when the job is submitted.
#
      elif '#PBS -A' in line and '##PBS -A' not in line:
        if self.project:             # on the command line
          newline = '#PBS -A '+self.project
          f.write(newline+'\n')
        else:                        # use default
          f.write(line)
#
# LSF queue, e.g.: #PBS -q regular (ys only)
# Avoid changing any commented pbs commands (double '##')
# Could test this for validity (premium, regular, etc), but
#   pbs will take care of this when the job is executed,
#   and queues may change. Also user's responsibility to
#   coordinate wallclock limit and queue.
#
      elif '#PBS -q' in line and '##PBS -q' not in line:
        if self.queue:             # on the command line
          newline = '#PBS -q '+self.queue
          f.write(newline+'\n')
        else:                        # use default
          f.write(line)
#
# Execution flag in job script:
#
      elif 'set execute' in line and 'set execdir' not in line:
        if self.execute:             # on the command line
          if self.execute == 'yes':
            newline = 'set execute = TRUE'
          else:
            newline = 'set execute = FALSE'
          f.write(newline+'\n')
        else:                        # use default
          f.write(line)
#
# Compiler (linux desktop platform only, not ch (not super)):
#
# if compiler=='intel', then makefile is Make.intel_hao64
# if compiler=='pgi', then makefile is Make.pgi_hao64
# if compiler=='gfort', then makefile is Make.gfort_hao64
#
      elif 'set make' in line:
        if self.machine != 'ch' and self.machine != 'de' and self.machine != 'pf':
          if self.compiler == 'intel':
            newline = 'set make = Make.intel_hao64'
          elif self.compiler == 'pgi':
            newline = 'set make = Make.pgi_hao64'
          elif self.compiler == 'gfort':
            newline = 'set make = Make.gfort_hao64'
          else:
            print('>>> Unknown compiler ',self.compiler)
            sys.exit()
        elif  self.machine == 'ch':
          newline = 'set make = Make.intel_ch'
        elif  self.machine == 'de':
          newline = 'set make = Make.intel_de'
        elif  self.machine == 'pf':
          newline = 'set make = Make.intel_pf'
        f.write(newline+'\n')
#
# Otherwise, no change to this line in default job script:
      else:               # no change
        f.write(line)

    f.close()
    os.popen('chmod u+x '+self.script_run)    # make it executable
#   print("Wrote job script ",self.script_run

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# def calc_wallclock

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  def submit(self,run_fullname):
    print("Submitting ",self.script_run," for run ",run_fullname)
    try:
      subprocess.call([self.script_run])
    except PermissionError:
      print("Setting execute permission for "+self.script_run)
      os.chmod(self.script_run, os.stat(self.script_run).st_mode | 0o111)
      subprocess.call([self.script_run])
  
#----------------------- Begin Class Run definition ---------------------
class Run(Job,Namelist):
  def __init__(self):
#
# There are 3 main categories of runs:
#   seasons (solstice,equinox, smin,smax)
#   storms (dec2006, nov2003, whi2008)
#   climatology (smin,smax)
#
# Short names and descriptions of available runs:
#
    self.names = [ \
      ['default_run','Default run'],
#
# Seasons:
      ['decsol_smax','December Solstice Solar Maximum'],
      ['decsol_smin','December Solstice Solar Minimum'],
      ['junsol_smax','June Solstice Solar Maximum'],
      ['junsol_smin','June Solstice Solar Minimum'],
      ['mareqx_smax','March Equinox Solar Maximum'],
      ['mareqx_smin','March Equinox Solar Minimum'],
      ['sepeqx_smax','September Equinox Solar Maximum'],
      ['sepeqx_smin','September Equinox Solar Minimum'],
#
# Storms:
      ['nov2003_heelis_gpi','November 2003 storm case, Heelis potential model, GPI data'],
      ['nov2003_weimer_imf','November 2003 storm case, Weimer potential model, IMF, GPI data'],
      ['dec2006_heelis_gpi','December 2006 "AGU storm", Heelis potential model, GPI data'],
      ['dec2006_weimer_imf','December 2006 "AGU storm", Weimer potential model, IMF and GPI data'],
      ['whi2008_heelis_gpi','2008 "Whole Heliosphere Interval", Heelis potential model, GPI data'],
      ['whi2008_weimer_imf','2008 "Whole Heliosphere Interval", Weimer potential model, IMF, GPI data'],
      ['jul2000_heelis_gpi','July 2000 "Bastille Day" storm, Heelis potential model, GPI data'],
      ['jul2000_weimer_imf','July 2000 "Bastille Day" storm, Weimer potential model, IMF, GPI data'],
#
# Climatology:
      ['climatology_smin','Climatology run with constant solar minimum conditions (Jan 1-5)'],
      ['climatology_smax','Climatology run with constant solar maximum conditions (Jan 1-5)']
      ]
    self.nruns = len(self.names)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  def get_number(self,run_name):
#
# Return list number given its short name.
#
    n=0
    for name in self.names:
      if name[0] in run_name: 
        self.number = n
        return self.number
      n=n+1
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  def get_name(self,run_number):
    if self.validate_number(run_number):
      self.name = self.names[int(run_number)][0]
      return self.name
    else:
      return ''
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  def validate_name(self,run_name):
#
# Return true (1) if run_name is valid (in the list of short names),
# return false (0) otherwise.
#
    for name in self.names:
      if name[0] in run_name: # "in" so user can modify (add to) the valid name
        return 1
    return 0
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  def validate_number(self,run_number):
#
# run.number must be an integer between 0 and run.nruns-1:
#   1. Try converting it to an integer, if this fails return 0
#   2. If it is an integer, and within range, return 1, otherwise return 0
# This is done silently, so calling function can print(its own errors.
#
    try:
      int(run_number)
      run_number = int(run_number)
      if run_number >= 0 and run_number <= self.nruns-1:
        return 1
      else:
        return 0
    except:
      return 0
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  def fullname(self,job,name):
#
# Run full names are model_name + model_res + short_name 
# (e.g. tiegcm_res5.0_decsol_smax)
#
    self.fullname = job.model_name+job.model_tag+'_res'+job.horires+'x'+job.vertres+'_z'+job.zitop+'_'+name 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# Print run numbers, short names, and descriptions:
# This will be used to prompt user for the desired run.
#
  def print_runs(self):
    print('\nThe following runs are available:')
    header = '''
NUMBER\tNAME\t\tDESCRIPTION
------\t----\t\t-----------'''
    print(header)
    n = 0
    while n < self.nruns:
      print(n,"\t",self.names[n][0],"\t",self.names[n][1])
      n = n+1
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  def set_run(self,job,tgcmdata):
    if job.machine == 'de' or job.machine == 'ch':
      prim_loc = '/glade/work/nikhilr/tiegcm3.0/prim'
    elif job.machine == 'pf':
      prim_loc = '/nobackup/nrao3/tiegcm/data/prim'

# Note: after going to tiegcm2.0, the model_res argument will no longer be necessary. 
#
# Set namelist modifications for each run:
# 
    n = int(self.number)
    self.name = self.names[n][0]
    self.desc = self.names[n][1]
    version = job.model_version
    res = job.model_res
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# Default run:
#
    if self.name == 'default_run':
      self.list_mods = [] # no changes to default namelist 
      self.list_rm   = []
      self.wc50_default = '00:30:00' # wallclock limit for 5.0-deg res (ch only)
      self.wc25_default = '01:15:00' # wallclock limit for 2.5-deg res (ch only)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# December Solstice, solar max:
#
    elif self.name == 'decsol_smax':

      source =  "'"+tgcmdata+"/"+self.fullname+"_prim.nc'"
      if not os.path.isfile(source):
        fullname = self.fullname.replace(job.model_tag,"")
        sec_source = job.stdout_dir+"/"+fullname+"_prim.nc"
        if not os.path.isfile(sec_source):
          fin = find_file('*'+self.name+'*', prim_loc)
          fout = job.stdout_dir+"/"+fullname+"_prim.nc"
          interpic(fin,float(job.horires),float(job.vertres),float(job.zitop),fout)
          source =  "'"+job.stdout_dir+"/"+fullname+"_prim.nc'"
        else:
          source =  "'"+job.stdout_dir+"/"+fullname+"_prim.nc'"
          

      self.list_mods = [ 
        ['LABEL'        , "'"+self.fullname+"'"],
        ['START_DAY'    , '355'],
        ['SOURCE'       , source],
        ['SOURCE_START' , '355 0 0 0'],
        ['PRISTART'        , '355 0 0 0'],
        ['PRISTOP'         , '360 0 0 0'],
        ['STEP'         , job.step],
        ['OUTPUT'       , "'"+job.hist_dir+self.fullname+"_prim_001.nc"+"'"],
        ['SECSTART'     , '355 1 0 0'],
        ['SECSTOP'      , '360 0 0 0'],
        ['SECOUT'       , "'"+job.hist_dir+self.fullname+"_sech_001.nc','to','"+job.hist_dir+self.fullname+"_sech_005.nc','by','1'"],
        ['POWER'        , '40.'],     
        ['CTPOTEN'      , '60.'],
        ['F107'         , '200.'],
        ['F107A'        , '200.']
        ]

      self.list_rm   = []
      self.wc50_default = '00:45:00' # wallclock limit for 5.0-deg res (ch only)
      self.wc25_default = '01:15:00' # wallclock limit for 2.5-deg res (ch only)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# December Solstice, solar min:
#
    elif self.name == 'decsol_smin':

    
      source =  "'"+tgcmdata+"/"+self.fullname+"_prim.nc'"
      if not os.path.isfile(source):
        fullname = self.fullname.replace(job.model_tag,"")
        sec_source = job.stdout_dir+"/"+fullname+"_prim.nc"
        if not os.path.isfile(sec_source):
          fin = find_file('*'+self.name+'*', prim_loc)
          fout = job.stdout_dir+"/"+fullname+"_prim.nc"
          interpic(fin,float(job.horires),float(job.vertres),float(job.zitop),fout)
          source =  "'"+job.stdout_dir+"/"+fullname+"_prim.nc'"
        else:
          source =  "'"+job.stdout_dir+"/"+fullname+"_prim.nc'"

      self.list_mods = [ 
        ['LABEL'        , "'"+self.fullname+"'"],
        ['START_DAY'    , '355'],
        ['SOURCE'       , source],
        ['SOURCE_START' , '355 0 0 0'],
        ['PRISTART'        , '355 0 0 0'],
        ['PRISTOP'         , '360 0 0 0'],
        ['STEP'         , job.step],
        ['OUTPUT'       , "'"+job.hist_dir+self.fullname+"_prim_001.nc"+"'"],
        ['SECSTART'     , '355 1 0 0'],
        ['SECSTOP'      , '360 0 0 0'],
        ['SECOUT'       , "'"+job.hist_dir+self.fullname+"_sech_001.nc','to','"+job.hist_dir+self.fullname+"_sech_005.nc','by','1'"],
        ['POWER'        , '18.'],     
        ['CTPOTEN'      , '30.'],
        ['F107'         , '70.'],
        ['F107A'        , '70.' ]
        ]

      self.list_rm = []
      self.wc50_default = '00:30:00' # wallclock limit for 5.0-deg res (ch only)
      self.wc25_default = '01:15:00' # wallclock limit for 2.5-deg res (ch only)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# June Solstice, solar max:
#
    elif self.name == 'junsol_smax':

      
      source =  "'"+tgcmdata+"/"+self.fullname+"_prim.nc'"
      if not os.path.isfile(source):
        fullname = self.fullname.replace(job.model_tag,"")
        sec_source = job.stdout_dir+"/"+fullname+"_prim.nc"
        if not os.path.isfile(sec_source):
          fin = find_file('*'+self.name+'*', prim_loc)
          fout = job.stdout_dir+"/"+fullname+"_prim.nc"
          interpic(fin,float(job.horires),float(job.vertres),float(job.zitop),fout)
          source =  "'"+job.stdout_dir+"/"+fullname+"_prim.nc'"
        else:
          source =  "'"+job.stdout_dir+"/"+fullname+"_prim.nc'"

      self.list_mods = [
        ['LABEL'        , "'"+self.fullname+"'"],
        ['START_DAY'    , '172'],
        ['SOURCE'       , source],
        ['SOURCE_START' , '172 0 0 0'],
        ['PRISTART'        , '172 0 0 0'],
        ['PRISTOP'         , '177 0 0 0'],
        ['STEP'         , job.step],
        ['OUTPUT'       , "'"+job.hist_dir+self.fullname+"_prim_001.nc"+"'"],
        ['SECSTART'     , '172 1 0 0'],
        ['SECSTOP'      , '177 0 0 0'],
        ['SECOUT'       , "'"+job.hist_dir+self.fullname+"_sech_001.nc','to','"+job.hist_dir+self.fullname+"_sech_005.nc','by','1'"],
        ['POWER'        , '40.'],     
        ['CTPOTEN'      , '60.'],
        ['F107'         , '200.'],
        ['F107A'        , '200.' ]
        ]

      self.list_rm   = []
      self.wc50_default = '00:30:00' # wallclock limit for 5.0-deg res (ch only)
      self.wc25_default = '01:15:00' # wallclock limit for 2.5-deg res (ch only)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# June Solstice, solar min:
#
    elif self.name == 'junsol_smin':

    
      source =  "'"+tgcmdata+"/"+self.fullname+"_prim.nc'"
      if not os.path.isfile(source):
        fullname = self.fullname.replace(job.model_tag,"")
        sec_source = job.stdout_dir+"/"+fullname+"_prim.nc"
        if not os.path.isfile(sec_source):
          fin = find_file('*'+self.name+'*', prim_loc)          
          fout = job.stdout_dir+"/"+fullname+"_prim.nc"
          interpic(fin,float(job.horires),float(job.vertres),float(job.zitop),fout)
          source =  "'"+job.stdout_dir+"/"+fullname+"_prim.nc'"
        else:
          source =  "'"+job.stdout_dir+"/"+fullname+"_prim.nc'"

      self.list_mods = [
        ['LABEL'        , "'"+self.fullname+"'"],
        ['START_DAY'    , '172'],
        ['SOURCE'       , source],
        ['SOURCE_START' , '172 0 0 0'],
        ['PRISTART'        , '172 0 0 0'],
        ['PRISTOP'         , '177 0 0 0'],
        ['STEP'         , job.step],
        ['OUTPUT'       , "'"+job.hist_dir+self.fullname+"_prim_001.nc"+"'"],
        ['SECSTART'     , '172 1 0 0'],
        ['SECSTOP'      , '177 0 0 0'],
        ['SECOUT'       , "'"+job.hist_dir+self.fullname+"_sech_001.nc','to','"+job.hist_dir+self.fullname+"_sech_005.nc','by','1'"],
        ['POWER'        , '18.'],     
        ['CTPOTEN'      , '30.'],
        ['F107'         , '70.'],
        ['F107A'        , '70.' ]
        ]

      self.list_rm   = []
      self.wc50_default = '00:30:00' # wallclock limit for 5.0-deg res (ch only)
      self.wc25_default = '01:15:00' # wallclock limit for 2.5-deg res (ch only)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# March Equinox, solar max:
#
    elif self.name == 'mareqx_smax':

    
      source =  "'"+tgcmdata+"/"+self.fullname+"_prim.nc'"
      if not os.path.isfile(source):
        fullname = self.fullname.replace(job.model_tag,"")
        sec_source = job.stdout_dir+"/"+fullname+"_prim.nc"
        if not os.path.isfile(sec_source):
          fin = find_file('*'+self.name+'*', prim_loc)          
          fout = job.stdout_dir+"/"+fullname+"_prim.nc"
          interpic(fin,float(job.horires),float(job.vertres),float(job.zitop),fout)
          source =  "'"+job.stdout_dir+"/"+fullname+"_prim.nc'"
        else:
          source =  "'"+job.stdout_dir+"/"+fullname+"_prim.nc'"

      self.list_mods = [ 
        ['LABEL'        , "'"+self.fullname+"'"],
        ['START_DAY'    , '80'],
        ['SOURCE'       , source],
        ['SOURCE_START' , '80 0 0 0'],
        ['PRISTART'        , '80 0 0 0'],
        ['PRISTOP'         , '85 0 0 0'],
        ['STEP'         , job.step],
        ['OUTPUT'       , "'"+job.hist_dir+self.fullname+"_prim_001.nc"+"'"],
        ['SECSTART'     , '80 1 0 0'],
        ['SECSTOP'      , '85 0 0 0'],
        ['SECOUT'       , "'"+job.hist_dir+self.fullname+"_sech_001.nc','to','"+job.hist_dir+self.fullname+"_sech_005.nc','by','1'"],
        ['POWER'        , '40.'],     
        ['CTPOTEN'      , '60.'],
        ['F107'         , '200.'],
        ['F107A'        , '200.' ]
        ]

      self.list_rm   = []
      self.wc50_default = '00:30:00' # wallclock limit for 5.0-deg res (ch only)
      self.wc25_default = '01:15:00' # wallclock limit for 2.5-deg res (ch only)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# March Equinox, solar min:
#
    elif self.name == 'mareqx_smin':

    
      source =  "'"+tgcmdata+"/"+self.fullname+"_prim.nc'"
      if not os.path.isfile(source):
        fullname = self.fullname.replace(job.model_tag,"")
        sec_source = job.stdout_dir+"/"+fullname+"_prim.nc"
        if not os.path.isfile(sec_source):
          fin = find_file('*'+self.name+'*', prim_loc)          
          fout = job.stdout_dir+"/"+fullname+"_prim.nc"
          interpic(fin,float(job.horires),float(job.vertres),float(job.zitop),fout)
          source =  "'"+job.stdout_dir+"/"+fullname+"_prim.nc'"
        else:
          source =  "'"+job.stdout_dir+"/"+fullname+"_prim.nc'"

      self.list_mods = [
        ['LABEL'        , "'"+self.fullname+"'"],
        ['START_DAY'    , '80'],
        ['SOURCE'       , source],
        ['SOURCE_START' , '80 0 0 0'],
        ['PRISTART'        , '80 0 0 0'],
        ['PRISTOP'         , '85 0 0 0'],
        ['STEP'         , job.step],
        ['OUTPUT'       , "'"+job.hist_dir+self.fullname+"_prim_001.nc"+"'"],
        ['SECSTART'     , '80 1 0 0'],
        ['SECSTOP'      , '85 0 0 0'],
        ['SECOUT'       , "'"+job.hist_dir+self.fullname+"_sech_001.nc','to','"+job.hist_dir+self.fullname+"_sech_005.nc','by','1'"],
        ['POWER'        , '18.'],     
        ['CTPOTEN'      , '30.'],
        ['F107'         , '70.'],
        ['F107A'        , '70.']
        ]

      self.list_rm   = []
      self.wc50_default = '00:30:00' # wallclock limit for 5.0-deg res (ch only)
      self.wc25_default = '01:15:00' # wallclock limit for 2.5-deg res (ch only)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# September Equinox, solar max:
#
    elif self.name == 'sepeqx_smax':

    
      source =  "'"+tgcmdata+"/"+self.fullname+"_prim.nc'"
      if not os.path.isfile(source):
        fullname = self.fullname.replace(job.model_tag,"")
        sec_source = job.stdout_dir+"/"+fullname+"_prim.nc"
        if not os.path.isfile(sec_source):
          fin = find_file('*'+self.name+'*', prim_loc)          
          fout = job.stdout_dir+"/"+fullname+"_prim.nc"
          interpic(fin,float(job.horires),float(job.vertres),float(job.zitop),fout)
          source =  "'"+job.stdout_dir+"/"+fullname+"_prim.nc'"
        else:
          source =  "'"+job.stdout_dir+"/"+fullname+"_prim.nc'"

      self.list_mods = [ 
        ['LABEL'        , "'"+self.fullname+"'"],
        ['START_DAY'    , '264'],
        ['SOURCE'       , source],
        ['SOURCE_START' , '264 0 0 0'],
        ['PRISTART'        , '264 0 0 0'],
        ['PRISTOP'         , '269 0 0 0'],
        ['STEP'         , job.step],
        ['OUTPUT'       , "'"+job.hist_dir+self.fullname+"_prim_001.nc"+"'"],
        ['SECSTART'     , '264 1 0 0'],
        ['SECSTOP'      , '269 0 0 0'],
        ['SECOUT'       , "'"+job.hist_dir+self.fullname+"_sech_001.nc','to','"+job.hist_dir+self.fullname+"_sech_005.nc','by','1'"],
        ['POWER'        , '40.'],     
        ['CTPOTEN'      , '60.'],
        ['F107'         , '200.'],
        ['F107A'        , '200.']
        ]

      self.list_rm   = []
      self.wc50_default = '00:30:00' # wallclock limit for 5.0-deg res (ch only)
      self.wc25_default = '01:15:00' # wallclock limit for 2.5-deg res (ch only)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# September Equinox, solar min:
#
    elif self.name == 'sepeqx_smin':

    
      source =  "'"+tgcmdata+"/"+self.fullname+"_prim.nc'"
      if not os.path.isfile(source):
        fullname = self.fullname.replace(job.model_tag,"")
        sec_source = job.stdout_dir+"/"+fullname+"_prim.nc"
        if not os.path.isfile(sec_source):
          fin = find_file('*'+self.name+'*', prim_loc)          
          fout = job.stdout_dir+"/"+fullname+"_prim.nc"
          interpic(fin,float(job.horires),float(job.vertres),float(job.zitop),fout)
          source =  "'"+job.stdout_dir+"/"+fullname+"_prim.nc'"
        else:
          source =  "'"+job.stdout_dir+"/"+fullname+"_prim.nc'"

      self.list_mods = [
        ['LABEL'        , "'"+self.fullname+"'"],
        ['START_DAY'    , '264'],
        ['SOURCE'       , source],
        ['SOURCE_START' , '264 0 0 0'],
        ['PRISTART'        , '264 0 0 0'],
        ['PRISTOP'         , '269 0 0 0'],
        ['STEP'         , job.step],
        ['OUTPUT'       , "'"+job.hist_dir+self.fullname+"_prim_001.nc"+"'"],
        ['SECSTART'     , '264 1 0 0'],
        ['SECSTOP'      , '269 0 0 0'],
        ['SECOUT'       , "'"+job.hist_dir+self.fullname+"_sech_001.nc','to','"+job.hist_dir+self.fullname+"_sech_005.nc','by','1'"],
        ['POWER'        , '18.'],     
        ['CTPOTEN'      , '30.'],
        ['F107'         , '70.'],
        ['F107A'        , '70.']
        ]

      self.list_rm   = []
      self.wc50_default = '00:30:00' # wallclock limit for 5.0-deg res (ch only)
      self.wc25_default = '01:15:00' # wallclock limit for 2.5-deg res (ch only)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# November 2003 storm case with Heelis and GPI:
#
    elif self.name == 'nov2003_heelis_gpi':

    
      source =  "'"+tgcmdata+"/"+self.fullname+"_prim.nc'"
      if not os.path.isfile(source):
        fullname = self.fullname.replace(job.model_tag,"")
        sec_source = job.stdout_dir+"/"+fullname+"_prim.nc"
        if not os.path.isfile(sec_source):
          fin = find_file('*'+self.name+'*', prim_loc)          
          fout = job.stdout_dir+"/"+fullname+"_prim.nc"
          interpic(fin,float(job.horires),float(job.vertres),float(job.zitop),fout)
          source =  "'"+job.stdout_dir+"/"+fullname+"_prim.nc'"
        else:
          source =  "'"+job.stdout_dir+"/"+fullname+"_prim.nc'"

      self.list_mods = [
        ['LABEL'        , "'"+self.fullname+"'"],
        ['START_YEAR'   , '2003'],
        ['START_DAY'    , '323'],
        ['SOURCE'       , source],
        ['SOURCE_START' , '323 0 0 0'],
        ['PRISTART'        , '323 0 0 0'],
        ['PRISTOP'         , '328 0 0 0'],
        ['STEP'         , job.step],
        ['OUTPUT'       , "'"+job.hist_dir+self.fullname+"_prim_001.nc"+"'"],
        ['SECSTART'     , '323 1 0 0'],
        ['SECSTOP'      , '328 0 0 0'],
        ['SECOUT'       , "'"+job.hist_dir+self.fullname+"_sech_001.nc','to','"+job.hist_dir+self.fullname+"_sech_005.nc','by','1'"],
        ['POTENTIAL_MODEL', "'HEELIS'"],
        ['GPI_NCFILE'     , "'"+tgcmdata+"/gpi_1960001-2015365.nc'"]
        ] 
#
# Pairs to be removed do not need values:
      self.list_rm = [
        ['POWER'        , '0.'],
        ['CTPOTEN'      , '0.'],
        ['F107'         , '0.'],
        ['F107A'        , '0.']]
      self.wc50_default = '00:30:00' # wallclock limit for 5.0-deg res (ch only)
      self.wc25_default = '01:15:00' # wallclock limit for 2.5-deg res (ch only)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# November 2003 storm case with Weimer and IMF, GPI:
#
    elif self.name == 'nov2003_weimer_imf':

     
      source =  "'"+tgcmdata+"/"+self.fullname+"_prim.nc'"
      if not os.path.isfile(source):
        fullname = self.fullname.replace(job.model_tag,"")
        sec_source = job.stdout_dir+"/"+fullname+"_prim.nc"
        if not os.path.isfile(sec_source):
          fin = find_file('*'+self.name+'*', prim_loc)          
          fout = job.stdout_dir+"/"+fullname+"_prim.nc"
          interpic(fin,float(job.horires),float(job.vertres),float(job.zitop),fout)
          source =  "'"+job.stdout_dir+"/"+fullname+"_prim.nc'"
        else:
          source =  "'"+job.stdout_dir+"/"+fullname+"_prim.nc'"
#
# As of 11/16/15, the current trunk code for tiegcm_res2.5 will crash in the 
#   first ~1.5 days if starting from tiegcm2.0 benchmark SOURCE history, and 
#   STEP is > 10 secs (e.g., 20, 30). It succeeds if timestep is reduced to 
#   10 secs, so force that here.
# (Interestingly, the current trunk code succeeds with STEP=30, if starting
#  from the old tiegcm1.95 benchmark SOURCE history)
#
        if job.model_res == '2.5': 
          if int(job.step) > 10: 
            print('NOTE ',self.name,': Changing timestep from ',job.step,' to 10 seconds')
            job.step = '10' # reduce timestep for res2.5 to 10 seconds

      self.list_mods = [
        ['LABEL'          , "'"+self.fullname+"'"],
        ['START_YEAR'     , '2003'],
        ['START_DAY'      , '323'],
        ['SOURCE'         , source],
        ['SOURCE_START'   , '323 0 0 0'],
        ['PRISTART'          , '323 0 0 0'],
        ['PRISTOP'           , '328 0 0 0'],
        ['STEP'           , job.step],
        ['OUTPUT'       , "'"+job.hist_dir+self.fullname+"_prim_001.nc"+"'"],
        ['SECSTART'       , '323 1 0 0'],
        ['SECSTOP'        , '328 0 0 0'],
        ['SECOUT'         , "'"+job.hist_dir+self.fullname+"_sech_001.nc','to','"+job.hist_dir+self.fullname+"_sech_005.nc','by','1'"],
        ['POTENTIAL_MODEL', "'WEIMER'"],
        ['IMF_NCFILE'     , "'"+tgcmdata+"/imf_OMNI_2003001-2003365.nc'"],
        ['GPI_NCFILE'     , "'"+tgcmdata+"/gpi_1960001-2015365.nc'"]
        ] 
#
# Pairs to be removed do not need values:
      self.list_rm = [
        ['POWER'        , '0.'],
        ['CTPOTEN'      , '0.'],
        ['F107'         , '0.'],
        ['F107A'        , '0.']]
      self.wc50_default = '00:30:00' # wallclock limit for 5.0-deg res (ch only)
      self.wc25_default = '02:30:00' # wallclock limit for 2.5-deg res (ch only) (step=10)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# December, 2006 AGU storm with Heelis and GPI:
#
    elif self.name == 'dec2006_heelis_gpi':

    
      source =  "'"+tgcmdata+"/"+self.fullname+"_prim.nc'"
      if not os.path.isfile(source):
        fullname = self.fullname.replace(job.model_tag,"")
        sec_source = job.stdout_dir+"/"+fullname+"_prim.nc"
        if not os.path.isfile(sec_source):
          fin = find_file('*'+self.name+'*', prim_loc)          
          fout = job.stdout_dir+"/"+fullname+"_prim.nc"
          interpic(fin,float(job.horires),float(job.vertres),float(job.zitop),fout)
          source =  "'"+job.stdout_dir+"/"+fullname+"_prim.nc'"
        else:
          source =  "'"+job.stdout_dir+"/"+fullname+"_prim.nc'"

      self.list_mods = [
        ['LABEL'          , "'"+self.fullname+"'"],
        ['START_YEAR'     , '2006'],
        ['START_DAY'      , '330'],
        ['SOURCE'         , source],
        ['SOURCE_START'   , '330 0 0 0'],
        ['PRISTART'          , '330 0 0 0'],
        ['PRISTOP'           , '360 0 0 0'],
        ['STEP'           , job.step],
        ['OUTPUT'         , "'"+job.hist_dir+self.fullname+"_prim_001.nc','to','"+job.hist_dir+self.fullname+"_prim_005.nc','by','1'"],
        ['SECSTART'       , '330 1 0 0'],
        ['SECSTOP'        , '360 0 0 0'],
        ['SECOUT'         , "'"+job.hist_dir+self.fullname+"_sech_001.nc','to','"+job.hist_dir+self.fullname+"_sech_032.nc','by','1'"],
        ['POTENTIAL_MODEL', "'HEELIS'"],
        ['GPI_NCFILE'     , "'"+tgcmdata+"/gpi_1960001-2015365.nc'"]
        ] 
#
# Pairs to be removed do not need values:
      self.list_rm = [
        ['POWER'        , '0.'],
        ['CTPOTEN'      , '0.'],
        ['F107'         , '0.'],
        ['F107A'        , '0.']]
      self.wc50_default = '12:00:00' # wallclock limit for 5.0-deg res (ch only)
      self.wc25_default = '12:00:00' # wallclock limit for 2.5-deg res (ch only)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# December, 2006 AGU storm with Weimer and IMF:
#
    elif self.name == 'dec2006_weimer_imf':

    
      source =  "'"+tgcmdata+"/"+self.fullname+"_prim.nc'"
      if not os.path.isfile(source):
        fullname = self.fullname.replace(job.model_tag,"")
        sec_source = job.stdout_dir+"/"+fullname+"_prim.nc"
        if not os.path.isfile(sec_source):
          fin = find_file('*'+self.name+'*', prim_loc)          
          fout = job.stdout_dir+"/"+fullname+"_prim.nc"
          interpic(fin,float(job.horires),float(job.vertres),float(job.zitop),fout)
          source =  "'"+job.stdout_dir+"/"+fullname+"_prim.nc'"
        else:
          source =  "'"+job.stdout_dir+"/"+fullname+"_prim.nc'"

      self.list_mods = [
        ['LABEL'          , "'"+self.fullname+"'"],
        ['START_YEAR'     , '2006'],
        ['START_DAY'      , '330'],
        ['SOURCE'         , source],
        ['SOURCE_START'   , '330 0 0 0'],
        ['PRISTART'          , '330 0 0 0'],
        ['PRISTOP'           , '360 0 0 0'],
        ['STEP'           , job.step],
        ['OUTPUT'         , "'"+job.hist_dir+self.fullname+"_prim_001.nc','to','"+job.hist_dir+self.fullname+"_prim_005.nc','by','1'"],
        ['SECSTART'       , '330 1 0 0'],
        ['SECSTOP'        , '360 0 0 0'],
        ['SECOUT'         , "'"+job.hist_dir+self.fullname+"_sech_001.nc','to','"+job.hist_dir+self.fullname+"_sech_032.nc','by','1'"],
        ['POTENTIAL_MODEL', "'WEIMER'"],
        ['IMF_NCFILE'     , "'"+tgcmdata+"/imf_OMNI_2006001-2006365.nc'"],
        ['GPI_NCFILE'     , "'"+tgcmdata+"/gpi_1960001-2015365.nc'"]
        ] 
#
# Pairs to be removed do not need values:
      self.list_rm = [
        ['POWER'        , '0.'],
        ['CTPOTEN'      , '0.'],
        ['F107'         , '0.'],
        ['F107A'        , '0.']]
      self.wc50_default = '12:00:00' # wallclock limit for 5.0-deg res (ch only)
      self.wc25_default = '12:00:00' # wallclock limit for 2.5-deg res (ch only)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# Whole Heliosphere Interval, 2008, with Heelis and GPI
#
    elif self.name == 'whi2008_heelis_gpi':

      source =  "'"+tgcmdata+"/"+self.fullname+"_prim.nc'"
      if not os.path.isfile(source):
        fullname = self.fullname.replace(job.model_tag,"")
        sec_source = job.stdout_dir+"/"+fullname+"_prim.nc"
        if not os.path.isfile(sec_source):
          fin = find_file('*'+self.name+'*', prim_loc)          
          fout = job.stdout_dir+"/"+fullname+"_prim.nc"
          interpic(fin,float(job.horires),float(job.vertres),float(job.zitop),fout)
          source =  "'"+job.stdout_dir+"/"+fullname+"_prim.nc'"
        else:
          source =  "'"+job.stdout_dir+"/"+fullname+"_prim.nc'"

      self.list_mods = [
        ['LABEL'          , "'"+self.fullname+"'"],
        ['START_YEAR'     , '2008'],
        ['START_DAY'      , '81'],
        ['SOURCE'         , source],
        ['SOURCE_START'   , '81 0 0 0'],
        ['PRISTART'          , '81 0 0 0'],
        ['PRISTOP'           , '106 0 0 0'],
        ['STEP'           , job.step],
        ['OUTPUT'         , "'"+job.hist_dir+self.fullname+"_prim_001.nc','to','"+job.hist_dir+self.fullname+"_prim_003.nc','by','1'"],
        ['SECSTART'       , '81 1 0 0'],
        ['SECSTOP'        , '106 0 0 0'],
        ['SECOUT'         , "'"+job.hist_dir+self.fullname+"_sech_001.nc','to','"+job.hist_dir+self.fullname+"_sech_025.nc','by','1'"],
        ['POTENTIAL_MODEL', "'HEELIS'"],
        ['GPI_NCFILE'     , "'"+tgcmdata+"/gpi_1960001-2015365.nc'"]
        ] 
#
# Lines (keys) to remove from default:
      self.list_rm = [
        ['POWER'        , '0.'],
        ['CTPOTEN'      , '0.'],
        ['F107'         , '0.'],
        ['F107A'        , '0.']]
#
# res5.0: 2 min/day * 25 days = 50 mins wc -> 1 hour
# res2.5: 8 min/day * 25 days = 200 mins -> 3.5 hours
#
      self.wc50_default = '12:00:00' # wallclock limit for 5.0-deg res (ch only)
      self.wc25_default = '12:00:00' # wallclock limit for 2.5-deg res (ch only)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# Whole Heliosphere Interval, 2008, with Weimer, IMF and GPI
#
    elif self.name == 'whi2008_weimer_imf':

    
      source =  "'"+tgcmdata+"/"+self.fullname+"_prim.nc'"
      if not os.path.isfile(source):
        fullname = self.fullname.replace(job.model_tag,"")
        sec_source = job.stdout_dir+"/"+fullname+"_prim.nc"
        if not os.path.isfile(sec_source):
          fin = find_file('*'+self.name+'*', prim_loc)          
          fout = job.stdout_dir+"/"+fullname+"_prim.nc"
          interpic(fin,float(job.horires),float(job.vertres),float(job.zitop),fout)
          source =  "'"+job.stdout_dir+"/"+fullname+"_prim.nc'"
        else:
          source =  "'"+job.stdout_dir+"/"+fullname+"_prim.nc'"

      self.list_mods = [
        ['LABEL'          , "'"+self.fullname+"'"],
        ['START_YEAR'     , '2008'],
        ['START_DAY'      , '81'],
        ['SOURCE'         , source],
        ['SOURCE_START'   , '81 0 0 0'],
        ['PRISTART'          , '81 0 0 0'],
        ['PRISTOP'           , '106 0 0 0'],
        ['STEP'           , job.step],
        ['OUTPUT'         , "'"+job.hist_dir+self.fullname+"_prim_001.nc','to','"+job.hist_dir+self.fullname+"_prim_003.nc','by','1'"],
        ['SECSTART'       , '81 1 0 0'],
        ['SECSTOP'        , '106 0 0 0'],
        ['SECOUT'         , "'"+job.hist_dir+self.fullname+"_sech_001.nc','to','"+self.fullname+"_sech_025.nc','by','1'"],
        ['POTENTIAL_MODEL', "'WEIMER'"],
        ['IMF_NCFILE'     , "'"+tgcmdata+"/imf_OMNI_2008001-2008366.nc'"],
        ['GPI_NCFILE'     , "'"+tgcmdata+"/gpi_1960001-2015365.nc'"]
        ] 
#
# Pairs to be removed do not need values:
      self.list_rm = [
        ['POWER'        , '0.'],
        ['CTPOTEN'      , '0.'],
        ['F107'         , '0.'],
        ['F107A'        , '0.']]
#
# res5.0: 2 min/day * 25 days = 50 mins wc -> 1 hour
# res2.5: 8 min/day * 25 days = 200 mins -> 3.5 hours
#
      self.wc50_default = '12:00:00' # wallclock limit for 5.0-deg res (ch only)
      self.wc25_default = '12:00:00' # wallclock limit for 2.5-deg res (ch only)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# July 2000 "Bastille Day" storm, with Heelis and GPI
#
    elif self.name == 'jul2000_heelis_gpi':

    
      source =  "'"+tgcmdata+"/"+self.fullname+"_prim.nc'"
      if not os.path.isfile(source):
        fullname = self.fullname.replace(job.model_tag,"")
        sec_source = job.stdout_dir+"/"+fullname+"_prim.nc"
        if not os.path.isfile(sec_source):
          fin = find_file('*'+self.name+'*', prim_loc)          
          fout = job.stdout_dir+"/"+fullname+"_prim.nc"
          interpic(fin,float(job.horires),float(job.vertres),float(job.zitop),fout)
          source =  "'"+job.stdout_dir+"/"+fullname+"_prim.nc'"
        else:
          source =  "'"+job.stdout_dir+"/"+fullname+"_prim.nc'"

      if job.model_res == '2.5': 
        if int(job.step) > 15: 
          print('NOTE ',self.name,': Changing timestep from ',job.step,' to 15 seconds')
          job.step = '15' # force timestep for res2.5 to 15 seconds

      self.list_mods = [
        ['LABEL'          , "'"+self.fullname+"'"],
        ['START_YEAR'     , '2000'],
        ['START_DAY'      , '192'],
        ['SOURCE'         , source],
        ['SOURCE_START'   , '192 0 0 0'],
        ['PRISTART'          , '192 0 0 0'],
        ['PRISTOP'           , '202 0 0 0'],
        ['STEP'           , job.step],
        ['OUTPUT'         , "'"+job.hist_dir+self.fullname+"_prim_001.nc','to','"+job.hist_dir+self.fullname+"_prim_003.nc','by','1'"],
        ['SECSTART'       , '192 1 0 0'],
        ['SECSTOP'        , '202 0 0 0'],
        ['SECOUT'         , "'"+job.hist_dir+self.fullname+"_sech_001.nc','to','"+job.hist_dir+self.fullname+"_sech_025.nc','by','1'"],
        ['POTENTIAL_MODEL', "'HEELIS'"],
        ['GPI_NCFILE'     , "'"+tgcmdata+"/gpi_1960001-2015365.nc'"]
        ] 
#
# Lines (keys) to remove from default:
      self.list_rm = [
        ['POWER'        , '0.'],
        ['CTPOTEN'      , '0.'],
        ['F107'         , '0.'],
        ['F107A'        , '0.']]
#
# res5.0: step=60 -> 2 min/day * 10 days = 20 mins wc -> 30 minutes
# res2.5: step=30 -> 17 min/day * 10 days = 170 mins -> 3.5 hours
#
      self.wc50_default = '00:30:00' # wallclock limit for 5.0-deg res (ch only)
      self.wc25_default = '01:00:00' # wallclock limit for 2.5-deg res (ch only)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# July 2000 "Bastille Day" storm, with Weimer, IMF, and GPI
#
    elif self.name == 'jul2000_weimer_imf':

    
      source =  "'"+tgcmdata+"/"+self.fullname+"_prim.nc'"
      if not os.path.isfile(source):
        fullname = self.fullname.replace(job.model_tag,"")
        sec_source = job.stdout_dir+"/"+fullname+"_prim.nc"
        if not os.path.isfile(sec_source):
          fin = find_file('*'+self.name+'*', prim_loc)          
          fout = job.stdout_dir+"/"+fullname+"_prim.nc"
          interpic(fin,float(job.horires),float(job.vertres),float(job.zitop),fout)
          source =  "'"+job.stdout_dir+"/"+fullname+"_prim.nc'"
        else:
          source =  "'"+job.stdout_dir+"/"+fullname+"_prim.nc'"
#
# Note: As of tiegcm/trunk -r1248 (1/31/16), STEP=10 and OPDIFFCAP=6e8 are
#       necessary for this run to succeed at 2.5-deg (resolution)
#
      opdiffcap = '0.'
      if job.model_res == '2.5': 
        if int(job.step) > 10: 
          print('NOTE ',self.name,': Changing timestep from ',job.step,' to 10 seconds')
          job.step = '10' # force timestep for res2.5 to 10 seconds
        opdiffcap = '6.e8'

      self.list_mods = [
        ['LABEL'          , "'"+self.fullname+"'"],
        ['START_YEAR'     , '2000'],
        ['START_DAY'      , '192'],
        ['SOURCE'         , source],
        ['SOURCE_START'   , '192 0 0 0'],
        ['PRISTART'          , '192 0 0 0'],
        ['PRISTOP'           , '202 0 0 0'],
        ['STEP'           , job.step],
        ['OUTPUT'         , "'"+job.hist_dir+self.fullname+"_prim_001.nc','to','"+job.hist_dir+self.fullname+"_prim_003.nc','by','1'"],
        ['SECSTART'       , '192 1 0 0'],
        ['SECSTOP'        , '202 0 0 0'],
        ['SECOUT'         , "'"+job.hist_dir+self.fullname+"_sech_001.nc','to','"+job.hist_dir+self.fullname+"_sech_025.nc','by','1'"],
        ['POTENTIAL_MODEL', "'WEIMER'"],
        ['IMF_NCFILE'     , "'"+tgcmdata+"/imf_OMNI_2000001-2000366.nc'"],
        ['GPI_NCFILE'     , "'"+tgcmdata+"/gpi_1960001-2015365.nc'"],
        ['OPDIFFCAP'      , opdiffcap],
        ] 
#
# Lines (keys) to remove from default:
      self.list_rm = [
        ['POWER'        , '0.'],
        ['CTPOTEN'      , '0.'],
        ['F107'         , '0.'],
        ['F107A'        , '0.']]
#
# res5.0: step=60 -> 2 min/day * 10 days = 20 mins wc -> 30 minutes
# res2.5: step=10 -> 25 min/day * 10 days = 250 mins -> 5.0 hours
#
      self.wc50_default = '00:30:00' # wallclock limit for 5.0-deg res (ch only)
      self.wc25_default = '12:00:00' # wallclock limit for 2.5-deg res (ch only)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# Climatology run with solar minimum:
# (Save daily primary histories (no secondaries)
# For now, just provide scripts for Jan 1-5. Then user can modify
#   dates, starting history, wallclock, etc for longer runs.
#
# ToDo: When full-year benchmark climatology runs are ready,
#       add comments about how to start on a different day, 
#       where to find source histories for arbitrary day of 
#       year, how to make longer runs, etc.
#
    elif self.name == 'climatology_smin':

    
      source =  "'"+tgcmdata+"/"+self.fullname+"_prim.nc'"
      if not os.path.isfile(source):
        fullname = self.fullname.replace(job.model_tag,"")
        sec_source = job.stdout_dir+"/"+fullname+"_prim.nc"
        if not os.path.isfile(sec_source):
          fin = find_file('*'+self.name+'*', prim_loc)          
          fout = job.stdout_dir+"/"+fullname+"_prim.nc"
          interpic(fin,float(job.horires),float(job.vertres),float(job.zitop),fout)
          source =  "'"+job.stdout_dir+"/"+fullname+"_prim.nc'"
        else:
          source =  "'"+job.stdout_dir+"/"+fullname+"_prim.nc'"
#
# 5-deg model at 60-sec timestep can simulate a full year in < 12 hours:
#   step=60 -> 1.8 min/day * 365 days = 657 mins / 60 = 11 hours
# 2.5-deg model at 30-sec timestep can complete 80 days in < 12 hours:
#   step=30 -> 8.2 min/day * 80 days = 656 mins / 60 = 11 hours
#
      stop = '366,0,0' # step=60 -> 1.8 min/day * 365 days = 657 mins / 60 = 11 hours
      if job.model_res == '2.5': 
        stop = '80,0,0'

      self.list_mods = [
        ['LABEL'          , "'"+self.fullname+"'"],
        ['START_YEAR'     , "2002  ! arbitrary for climatology run"], 
        ['START_DAY'      , '1'],
        ['SOURCE'         , source],
        ['SOURCE_START'   , '1 0 0 0'],
        ['PRISTART'          , '1 0 0 0'],
        ['PRISTOP'           , stop],
        ['STEP'           , job.step],
        ['OUTPUT'         , "'"+job.hist_dir+self.fullname+"_prim_001.nc','to','"+job.hist_dir+self.fullname+"_prim_020.nc','by','1'"],
        ['MXHIST_PRIM'    , '20'],     # default is 10
        ['POTENTIAL_MODEL', "'HEELIS'"],
        ['POWER'        , '18.'],     
        ['CTPOTEN'      , '30.'],
        ['F107'         , '70.'],
        ['F107A'        , '70.']]
#
# Pairs to be removed do not need values:
      self.list_rm = [
        ['SECSTART'    ,''],
        ['SECSTOP'     ,''],
        ['SECHIST'     ,''],
        ['SECOUT'      ,''],
        ['MXHIST_SECH' ,''],
        ['SECFLDS'     ,'']]
#
      self.wc50_default = '08:00:00' # wallclock limit for 5.0-deg res (ch only)
      self.wc25_default = '12:00:00' # wallclock limit for 2.5-deg res (ch only)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# Climatology run with solar maximum:
# Solar max run can crash around day 85 if STEP > 20 seconds.
#
    elif self.name == 'climatology_smax':

    
      source =  "'"+tgcmdata+"/"+self.fullname+"_prim.nc'"
      if not os.path.isfile(source):
        fullname = self.fullname.replace(job.model_tag,"")
        sec_source = job.stdout_dir+"/"+fullname+"_prim.nc"
        if not os.path.isfile(sec_source):
          fin = find_file('*'+self.name+'*', prim_loc)          
          fout = job.stdout_dir+"/"+fullname+"_prim.nc"
          interpic(fin,float(job.horires),float(job.vertres),float(job.zitop),fout)
          source =  "'"+job.stdout_dir+"/"+fullname+"_prim.nc'"
        else:
          source =  "'"+job.stdout_dir+"/"+fullname+"_prim.nc'"
#
# 5-deg model at 60-sec timestep can simulate a full year in < 12 hours:
#   step=60 -> 1.8 min/day * 365 days = 657 mins / 60 = 11 hours
#
# Climatology smax at 2.5-deg requires a 20-sec timestep:
# 2.5-deg model at 20-sec timestep can complete 55 days in < 12 hours:
#   step=20 -> 12.2 min/day * 55 days = 671 mins / 60 = 11 hours
#
      stop = '366,0,0' # step=60 -> 1.8 min/day * 365 days = 657 mins / 60 = 11 hours
      if job.model_res == '2.5': 
        if int(job.step) > 20: 
          print('NOTE ',self.name,': Changing timestep from ',job.step,' to 20 seconds')
          job.step = '20' # force timestep for res2.5 to 20 seconds
        stop = '55,0,0'

      self.list_mods = [
        ['LABEL'          , "'"+self.fullname+"'"],
        ['START_YEAR'     , "2002  ! arbitrary for climatology run"], 
        ['START_DAY'      , '1'],
        ['SOURCE'         , source],
        ['SOURCE_START'   , '1 0 0 0'],
        ['PRISTART'          , '1 0 0 0'],
        ['PRISTOP'           , stop],
        ['STEP'           , job.step],
        ['OUTPUT'         , "'"+job.hist_dir+self.fullname+"_prim_001.nc','to','"+job.hist_dir+self.fullname+"_prim_020.nc','by','1'"],
        ['MXHIST_PRIM'    , '20'],     # default is 10
        ['POTENTIAL_MODEL', "'HEELIS'"],
        ['POWER'        , '40.'],     
        ['CTPOTEN'      , '60.'],
        ['F107'         , '200.'],
        ['F107A'        , '200.']]
#
# Pairs to be removed do not need values:
      self.list_rm = [
        ['SECSTART'    ,''],
        ['SECSTOP'     ,''],
        ['SECHIST'     ,''],
        ['SECOUT'      ,''],
        ['MXHIST_SECH' ,''],
        ['SECFLDS'     ,'']]
#
      self.wc50_default = '08:00:00' # wallclock limit for 5.0-deg res (ch only)
      self.wc25_default = '12:00:00' # wallclock limit for 2.5-deg res (ch only)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
