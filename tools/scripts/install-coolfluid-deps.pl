#!/usr/bin/env perl

#==========================================================================
# TODO
#==========================================================================
# * user config file for overiding options

#==========================================================================
# Modules
#==========================================================================
use Term::ANSIColor;
use Getopt::Long;
use File::Path;
use File::Copy;

#==========================================================================
# Constants
#==========================================================================
my $ERRORCOLOR="bold red";
my $OKCOLOR="bold green";
my $HEADINGCOLOR = "bold";
my $DEBUGCOLOR = "yellow";
my $WARNCOLOR = "bold yellow";

#==========================================================================
# Global Variables
#==========================================================================

my $home = $ENV{HOME};
my $user = $ENV{USER};
my $arch = get_arch();

my $opt_help          = 0;
my $opt_dryrun        = 0;
my $opt_nocolor       = 0;
my $opt_envars        = 0;
my $opt_genconf       = 0;
my $opt_debug         = "1";
my $opt_int64         = "0";
my $opt_cuda_dir      = "";
my $opt_viennacl      = "0";
my $opt_nompi         = 0;
my $opt_mpi           = "openmpi";
my $opt_mpi_dir       = "";
my $opt_dependencies  = 0;
my $opt_sources       = 0;
my $opt_compile       = 0;
my $opt_fetchonly     = 0;
my $opt_many_mpi      = 0;
my $opt_with_cuda     = 0;
my $opt_static        = "0";
my $opt_install_dir   = "$home/local/coolfluid_deps/$arch";
#my $opt_install_dir   = "$home/COOLFluiD_DEPS";
my $opt_install_mpi_dir = "";
my $opt_petsc_dir = "";
my $opt_blaslapack_dir = "";
my $opt_parmetis_dir = "";
my $opt_cuda_bin  = "";
my $opt_cmake_dir     = "";
my $opt_confile       = "coolfluid.conf";
my $opt_tmp_dir       = "$home/coolfluid_tmp";
my $opt_build         = "optim";
#my $opt_dwnldsrc      = "http://coolfluidsrv.vki.ac.be/webfiles/coolfluid/packages";
my $opt_dwnldsrc      = "https://github.com/andrealani/COOLFluiD/tree/master/packages";
my $opt_wgetprog      = "wget -nc -nd";
my $opt_curlprog      = "curl -O -nc -nv --progress-bar";
my $opt_dwnldprog     = $opt_wgetprog;
my $opt_makeopts      = "-j12";
my $opt_svnrevision   = 0;
my $opt_cray          = 0;
my $opt_ibmgnu        = 0;
my @opt_install = ();

# list of packages, and their associated values
# [$vrs] : default version to install
# [$dft] : install by default ?
# [$ins] : should it be installed ?
# [$pri] : installation priority
# [$fnc] : function that implements the installation

my $vrs = 0;
my $dft = 1;
my $ins = 2;
my $pri = 3;
my $fnc = 4;

$priority = 0;

# these packages are listed by priority
my %packages = (  #  version   default install priority      function
    "coreutils"  => [ "6.7",    'off',  'off', $priority++,  sub { install_gnu("coreutils") } ],
    "make"       => [ "3.81",   'off',  'off', $priority++,  sub { install_gnu("make") } ],
#    "cmake"      => [ "2.8.11.2",  'on' ,  'off', $priority++,  \&install_cmake ],
#    "cmake"      => [ "3.10.0",  'on' ,  'off', $priority++,  \&install_cmake ],
     "cmake"      => [ "3.17.0",  'on' ,  'off', $priority++,  \&install_cmake ],
    "wget"       => [ "1.11.4", 'off',  'off', $priority++,  \&install_wgetprog],
    "binutils"   => [ "2.16.1", 'off',  'off', $priority++,  sub { install_gnu("binutils") } ],
    "m4"         => [ "1.4.4",  'off',  'off', $priority++,  sub { install_gnu("m4") } ],
    "tar"        => [ "1.15",   'off',  'off', $priority++,  sub { install_gnu("tar") } ],
    "gcc4"       => [ "4.1.1",  'off',  'off', $priority++,  \&install_gcc4 ],
    "gcc3"       => [ "3.4.5",  'off',  'off', $priority++,  \&install_gcc3 ],
    "autoconf"   => [ "2.59",   'off',  'off', $priority++,  sub { install_gnu("autoconf") } ],
    "automake"   => [ "1.9.5",  'off',  'off', $priority++,  sub { install_gnu("automake") } ],
    "libtool"    => [ "1.5.14", 'off',  'off', $priority++,  sub { install_gnu("libtool") } ],
    "openssl"    => [ "0.9.7i", 'off',  'off', $priority++,  \&install_openssl ],
    "blas"       => [ "3.0.3",  'off',  'off', $priority++,  \&install_blas ],
    "clapack"    => [ "",       'off',  'off', $priority++,  \&install_clapack ],
    "lapack"     => [ "3.0.3",  'off',  'off', $priority++,  \&install_lapack ],
    "log4cpp"    => [ "0.3.4b", 'off',  'off', $priority++,  sub { install_gnu("log4cpp") } ],
    "zlib"       => [ "1.2.3",  'off',  'off', $priority++,  sub { install_gnu("zlib") } ],
    "cppunit"    => [ "1.10.2", 'off',  'off', $priority++,  sub { install_gnu("cppunit") } ],
    "dateshift"  => [ "1.0",    'off',  'off', $priority++,  sub { install_gnu("dateshift") } ],
    "curl"       => [ "7.19.7", 'off' ,  'off', $priority++,  \&install_curl ],
    "libfaketime"=> [ "0.8",    'off',  'off', $priority++,  \&install_libfaketime ],
#     "boost"      => [ "1_42_0", 'on' ,  'off', $priority++,  \&install_boost ],
#    "boost"      => [ "1_59_0", 'on' ,  'off', $priority++,  \&install_boost ],
#    "boost"      => [ "1_54_0", 'on' ,  'off', $priority++,  \&install_boost ],
#    "boost"      => [ "1_66_0", 'on' ,  'off', $priority++,  \&install_boost ],
#    "boost"      => [ "1_70_0", 'on' ,  'off', $priority++,  \&install_boost ],
#    "boost"      => [ "1_79_0", 'on' ,  'off', $priority++,  \&install_boost ],
    "boost"      => [ "1_76_0", 'on' ,  'off', $priority++,  \&install_boost ],
#    "boost"      => [ "1_72_0", 'on' ,  'off', $priority++,  \&install_boost ],
#   "openmpi"    => [ "1.4.2",  'off',  'off', $priority++,  \&install_openmpi ],
#    "openmpi"    => [ "1.6.5",  'off',  'off', $priority++,  \&install_openmpi ],
#    "openmpi"    => [ "1.10.0",  'off',  'off', $priority++,  \&install_openmpi ],
#    "openmpi"    => [ "3.0.0",  'off',  'off', $priority++,  \&install_openmpi ],
     "openmpi"    => [ "4.0.3",  'off',  'off', $priority++,  \&install_openmpi ],
#    "mpich"      => [ "3.1b1",'off',  'off', $priority++,  \&install_mpich ],    	
    "mpich"      => [ "3.1.3",'off',  'off', $priority++,  \&install_mpich ], 
    "mpich2"     => [ "1.4.1p1",  'off',  'off', $priority++,  \&install_mpich2 ],
#    "mvapich2"   => [ "1.9",      'off',  'off', $priority++,  \&install_mvapich2 ],
     "mvapich2"  => [ "2.0.1", 'off',  'off', $priority++,  \&install_mvapich2 ], 
    "parmetis"   => [ "4.0.3",  'on' ,  'off', $priority++,  \&install_parmetis ],
    # "parmetis"   => [ "3.1.1",  'on' ,  'off', $priority++,  \&install_parmetis ],
    "hdf5"       => [ "1.6.4",  'off',  'off', $priority++,  \&install_hdf5 ],
    "subversion" => [ "1.4.3",  'off',  'off', $priority++,  \&install_subversion ],
    "trilinos"   => [ "10.10.2",'off',  'off', $priority++,  \&install_trilinos ],
#    "petsc"      => [ "3.2-p6",'on',  'off', $priority++,  \&install_petsc ], 
#    "petsc"      => [ "3.4.2",'on',  'off', $priority++,  \&install_petsc ], 
#    "petsc"      => [ "3.6.3-dev",'on',  'off', $priority++,  \&install_petsc ],
#    "petsc"      => [ "3.6.3",'on',  'off', $priority++,  \&install_petsc ], 
#     "petsc"      => [ "3.9.0",'on',  'off', $priority++,  \&install_petsc ],
#    "petsc"      => [ "3.9.3",'on',  'off', $priority++,  \&install_petsc ],
    "petsc"      => [ "3.12.3",'on',  'off', $priority++,  \&install_petsc ],
#    "petsc"      => [ "3.7.6",'on',  'off', $priority++,  \&install_petsc ],
#    "petsc"      => [ "3.7.3a",'on',  'off', $priority++,  \&install_petsc ],
#    "petsc"      => [ "3.7.3-next",'on',  'off', $priority++,  \&install_petsc ],
    "gmsh"       => [ "1.60.1", 'off',  'off', $priority++,  sub { install_gnu("gmsh") } ],
    "ccache"     => [ "2.4",    'off',  'off', $priority++,  sub { install_gnu("ccache") } ],
    "distcc"     => [ "2.18.3", 'off',  'off', $priority++,  sub { install_gnu("distcc") } ],
    "cgnslib"    => [ "2.5-4",  'off',  'off', $priority++,  \&install_cgnslib ],
    "cgnstools"  => [ "2-5-4",  'off',  'off', $priority++,  \&install_cgnstools ],
    "google-perftools" => [ "0.93",'off','off', $priority++, \&install_google_perftools ],
);

#==========================================================================
# Command Line
#==========================================================================

sub parse_commandline() # Parse command line
{
    $opt_help=1 unless GetOptions (
        'help'                  => \$opt_help,
        'nocolor'               => \$opt_nocolor,
        'envars'                => \$opt_envars,
        'genconf'               => \$opt_genconf,
        'debug=s'               => \$opt_debug,
        'int64=s'               => \$opt_int64,
        'static=s'              => \$opt_static,
        'viennacl=s'            => \$opt_viennacl,
	'nompi'                 => \$opt_nompi,
        'many-mpi'              => \$opt_many_mpi,
        'mpi=s'                 => \$opt_mpi,
        'mpi-dir=s'             => \$opt_mpi_dir,
        'fetchonly'             => \$opt_fetchonly,
        'dry-run'               => \$opt_dryrun,
        'install-dir=s'         => \$opt_install_dir,
        'install-mpi-dir=s'     => \$opt_install_mpi_dir,
        'cmake-dir=s'           => \$opt_cmake_dir,
        'cuda-dir=s'            => \$opt_cuda_dir,
	'tmp-dir=s'             => \$opt_tmp_dir,
        'dwnldsrc=s'            => \$opt_dwnldsrc,
        'branch=s'              => \$opt_branch,
        'build=s'               => \$opt_build,
        'makeopts=s'            => \$opt_makeopts,
        'install=s'             => \@opt_install,
        'install-petsc-dir=s'   => \$opt_petsc_dir,
	'install-parmetis-dir=s'   => \$opt_parmetis_dir,
        'cray'                  => \$opt_cray,
        'ibmgnu'                => \$opt_ibmgnu,
	'blaslapack-dir=s'      => \$opt_blaslapack_dir,
    );

    # show help if required
    if ($opt_help != 0)
    {
      print <<ZZZ;
install-coolfluid-deps.pl : Install software dependencies for COOLFluiD

usage: install-coolfluid-deps.pl [options]

By default will install a recomended set of dependencies: cmake,curl,boost,openmpi,parmetis,petsc

options:
        --help               Show this help.
        --nocolor            Don't color output
        --envars             Place enviromental variables into .bashrc

        --genconf            Generate a coolfluid.conf based on the installed dependencies

        --nompi              Don't compile with mpi support. This is only active for some packages.
        --mpi=               MPI compiler to use for compilations
                              Default: $opt_mpi.
        --many-mpi=          Install all mpi related packages in a separate directory
                             therefore allowing multiple mpi environments to coexist
                              Default: $opt_many_mpi.
       
        --cray               Compile on CRAY systems
        --ibmgnu             Compile on IBM systems with GNU compiler       
        --static=            Compile all libraries static
        --debug=             Compile dependencies and coolfluid with debug symbols (=1 is default)
        --int64=             If =1 compiles PETSc using long long int. 
        --fetchonly          Just download the sources. Do not install anything.
        --dry-run            Don't actually perform the configuration.
                             Just output what you would do.

        --blaslapack-dir     Location of blas and lapack libs
        --mpi-dir=           Location of existing MPI installation
        --install-dir=       Location of the software installation directory
                              Default: $opt_install_dir
        --install-mpi-dir=   Location where to install MPI
                              Default: $opt_install_mpi_dir
        --install-petsc-dir= Location of the directory where to install petsc (by default is $opt_install_mpi_dir/petsc )
        --install-parmetis-dir= Location of the directory where to install parmetis (by default is $opt_install_mpi_dir )
        --cmake-dir=         Location for the cmake installation
        --cuda-dir=          Location of the cuda installation
        --tmp-dir=           Location of the temporary directory for compilation
                              Default: $opt_tmp_dir

        --dwnldsrc=          URL of download server from where to download sources of dependencies
                              Default: $opt_dwnldsrc

        --install            Comma separated list of packages to install in the specified order.
                             Every test will be run for on the number of CPUs specified here.
                             Example: --install=all,hdf5
                              Default: all

Examples of usage:

# full installation of default packages (curl, cmake, boost, openmpi, petsc, parmetis) in default locations
./install-coolfluid-deps.pl

# full installation of default packages (curl, cmake, boost, openmpi, petsc, parmetis) in prescribed locations
# ./install-coolfluid-deps.pl --install-dir=\$opt_install_dir --install-mpi-dir=\$opt_install_mpi_dir --tmp-dir=\$opt_tmp_dir

# full installation of default packages plus CUDA bindings (thrust, cusp libraries needed by petsc) with some user predefined locations 
./install-coolfluid-deps.pl --install-dir=\$opt_install_dir --install-mpi-dir=\$opt_install_mpi_dir --tmp-dir=\$opt_tmp_dir --cuda-dir=\$opt_cuda_dir

#installation of individual (MPI-dependent) packages in \$opt_install_mpi_dir from \$opt_tmp_dir
./install-coolfluid-deps.pl --install=mpich2,petsc,parmetis --install-mpi-dir=\$opt_install_mpi_dir --tmp-dir=\$opt_tmp_dir

# installation of petsc in optimized mode and using long long int
./install-coolfluid-deps.pl --install=petsc --int64=1 --install-petsc-dir=\$opt_install_mpi_dir/petsc_optim_long --debug=0 --install-mpi-dir=\$opt_install_mpi_dir

ZZZ
    exit(0);
    }
    @opt_install = split(/,/,join(',',@opt_install));

   
   # gory fix to circumvent downloading from internal server
   run_command("mkdir $opt_tmp_dir"); 
   run_command_or_die("cp ../../packages/* $opt_tmp_dir");
   
   #run_command_or_die("mkdir $opt_tmp_dir");
   #my $status = get_command_status("cp ../../packages/* $opt_tmp_dir");
   #print my_colored("STATUS IS $status\n",$OKCOLOR);
   #if ($status eq 0)
   #{
   # run_command_or_die("svn co https://github.com/andrealani/COOLFluiD/trunk/packages $opt_tmp_dir");
   # run_command_or_die("mv $opt_tmp_dir/packages/* $opt_tmp_dir ; rm -fr packages");
   #}

}

#==========================================================================
# Helper funcions
#==========================================================================

sub my_colored ($$)
{
  return ($opt_nocolor ? shift : colored($_[0], $_[1]));
}

#==========================================================================

sub rm_file ($)
{
  my ($file) = @_;
  unlink($file) || warn "warn: not deleting $file: $!";
}

#==========================================================================

sub get_command_status($)
{
    my ($args)=@_;
    print my_colored("Executing   : $args\n",$OKCOLOR);
    unless ($opt_dryrun) {
        my $status = system($args);
        return $status;
    }
    return 0;
}

#==========================================================================

sub run_command_or_die($)
{
    my ($args)=@_;
    print my_colored("Executing   : $args\n",$OKCOLOR);
    unless ($opt_dryrun) {
        my $status = system($args);
        print my_colored("Exit Status : $status\n",$OKCOLOR);
        die "$args exited with error" unless $status == 0;
    }
}

#==========================================================================

sub run_command($)
{
    my ($args)=@_;
    my $output;
    # print my_colored("Executing : $args",$OKCOLOR);
    my $command = join("",$args,"|");
    my $pid=open READER, $command or die "Can't run the program: $args $!\n";
    while(<READER>){
       $output.=$_;
    }
    close READER;
    # print my_colored($output,$OK_COLOR);
    return $output;
}

#==========================================================================

sub safe_chdir($)
{
    my ($dir)=@_;
    print my_colored("Changing to dir $dir\n",$DEBUGCOLOR);
    chdir($dir) or die "Cannot chdir to $dir ($!)";
}

#==========================================================================

sub safe_copy($$)
{
    my ($orig,$targ)=@_;
    copy ($orig,$targ) or die "Cannot copy $orig to $targ ($!)";
}

#==========================================================================

sub safe_delete($)
{
    unlink("$_") or die "Failed to delete file $_\n";
}

#==========================================================================

sub get_arch() # returns the current architecture
{
    my $args="uname -m";
    my $arch = run_command($args);
    chomp($arch);
    return $arch;
}

#==========================================================================

sub is_mac()
{
    my $args="uname -s";
    my $arch = run_command($args);
    chomp($arch);
    if ($arch =~ Darwin) {
        return 1;
    } else {
        return 0;
    }
}

#==========================================================================

sub print_var($$) # create a recursive dir path
{
    my ($var,$value)=@_;
    print my_colored($var,$OKCOLOR); print " : ";
    print my_colored($value,$DEBUGCOLOR); print "\n";
}

#==========================================================================

sub parse_config_file($) # parse the config file to get the user overiding options
{
    my ($filename)=@_;
    open CONFIG, "<", $filename or die ("Error opening config file $filename!\n");

    while (<CONFIG>) {
        chomp;                  # no newline
        s/#.*//;                # no comments
        s/^\s+//;               # no leading white
        s/\s+$//;               # no trailing white
        next unless length;     # anything left?
        my ($var, $value) = split(/\s*=\s*/, $_, 2);
        $user_pref{$var} = $value;
    }

    close CONFIG;
}

#==========================================================================
# Local functions
#==========================================================================

sub domkdir ()
{
    mkpath "$opt_install_dir";
    mkpath "$opt_install_dir/bin";
    mkpath "$opt_install_dir/local";
    mkpath "$opt_install_dir/lib";
    mkpath "$opt_install_dir/include";
    mkpath "$opt_install_dir/share";
    mkpath "$opt_install_mpi_dir";
    mkpath "$opt_cmake_dir";
    mkpath "$opt_tmp_dir";
}

#==========================================================================

sub prepare ()
{
    # set the default mpi
    $packages{$opt_mpi}[$dft] = 'on';

    # set the mpi install dir if the user did not set
    if ($opt_install_mpi_dir eq "")
    {
      if ($opt_many_mpi)
      {
        $version = $packages{"$opt_mpi"}[$vrs];
        $opt_install_mpi_dir = "$opt_install_dir/mpi/$opt_mpi-$version";
      } else {
        $opt_install_mpi_dir = $opt_install_dir;
      }
    }

    # set the mpi dir if the user did not set
    if ($opt_mpi_dir eq "")
    {
        $opt_mpi_dir = $opt_install_mpi_dir;	  
    }

    # set the cmake dir if the user did not set
    if ($opt_cmake_dir eq "")
    {
      $opt_cmake_dir = $opt_install_dir;
    }

    # make directories for installation
    unless ($opt_dryrun) { domkdir(); }

    # normal paths
    $ENV{PATH} = "$opt_install_dir/bin:" . $ENV{PATH};
    $ENV{LD_LIBRARY_PATH} = "$opt_install_dir/lib:" . $ENV{LD_LIBRARY_PATH};

    # mpi specific paths
    $ENV{PATH} = "$opt_mpi_dir/bin:" . $ENV{PATH};
    $ENV{LD_LIBRARY_PATH} = "$opt_mpi_dir/lib:" . $ENV{LD_LIBRARY_PATH};

    $ENV{CFLAGS}   = "-O3 " . $ENV{CFLAGS};
    $ENV{CXXFLAGS} = "-O3 " . $ENV{CXXFLAGS};
    $ENV{FFLAGS}   = "-O3 " . $ENV{FFLAGS};
    $ENV{F77FLAGS} = $ENV{FFLAGS};
    $ENV{F90FLAGS} = $ENV{FFLAGS};

    if ($arch eq "x86_64" )
    {
        $ENV{CFLAGS}   = "-fPIC " . $ENV{CFLAGS};
        $ENV{CXXFLAGS} = "-fPIC " . $ENV{CXXFLAGS};
        $ENV{FFLAGS}   = "-fPIC " . $ENV{FFLAGS};
        $ENV{F77FLAGS}  = "-fPIC " . $ENV{F77FLAGS};
        $ENV{F90FLAGS}  = "-fPIC " . $ENV{F90FLAGS};
    }

    if ( !(exists $ENV{CC}) )
    {
      $ENV{CC} = "gcc";
      print "Setting C compiler to \"".$ENV{CC}."\". Override this with environment variable \"CC\"\n";
    }

    if ( !(exists $ENV{CXX}) )
    {
      $ENV{CXX} = "g++";
      print "Setting C++ compiler to \"".$ENV{CXX}."\". Override this with environment variable \"CXX\"\n";
    }

    if (!((exists $ENV{FC}) or (exists $ENV{F77})))
    {
      $ENV{FC} = "gfortran";
      print "Setting Fortran compiler to \"".$ENV{FC}."\". Override this with environment variable \"FC\".\n";
    }

    # makes sure the both compiler variable F77 and FC always exist
    if ( !(exists $ENV{FC}) )
    {
      print "Setting FC equal to F77\n";
      $ENV{FC} = $ENV{F77};
    }
    if ( !(exists $ENV{F77}) )
    {
      print "Setting F77 equal to FC\n";
      $ENV{F77} = $ENV{FC};
    }
}

#==========================================================================

sub generate_conf ()
{
   if (-e $opt_confile)
   {
	print my_colored("\nFile $opt_confile exists, not generating\n",$SECTIONCOLOR);
   }
   else
   {
   print my_colored("\nGenerating $opt_confile\n",$SECTIONCOLOR);

   unless ($opt_dryrun)
   {
   open CONF, ">", "$opt_confile" or die ("Error creating $opt_confile!\n");
   
   if ($opt_petsc_dir eq "") {
     $opt_petsc_dir = "$opt_install_mpi_dir";
   }
 
   if ($opt_parmetis_dir eq "") {
     $opt_parmetis_dir = "$opt_install_mpi_dir";
   }
 
 print CONF  <<ZZZ;

	cc     = $opt_mpi_dir/bin/mpicc
	cxx    = $opt_mpi_dir/bin/mpic++
	fc     = $ENV{FC}
        cudac  = $opt_cuda_bin

	fflags    = -O3
	cflags    = -O3
	cxxflags  = -O3

	mpi_dir       = $opt_mpi_dir
	boost_dir     = $opt_install_dir
	petsc_dir     = $opt_petsc_dir
	parmetis_dir  = $opt_parmetis_dir
        cuda_dir      = $opt_cuda_dir

    withcuda = $opt_with_cuda

	mods-getall = 0

	mod_AnalyticalModel    = on
	mod_CFmesh2THOR        = on
	mod_CFmeshCellSplitter = on
	mod_CFmeshExtruder     = on
	mod_CFmeshFileReader   = on
	mod_CFmeshFileWriter   = on
	mod_CGNS2CFmesh        = on
	mod_Gambit2CFmesh      = on
	mod_Solver             = on
	mod_THOR2CFmesh        = on
	mod_TecplotWriter      = on
	mod_Petsc              = on
	mod_XCFcaseConverter   = on

ZZZ
   close CONF;
   }
 }
}

#==========================================================================

sub post () {
  if ($opt_envars)
  {
    run_command_or_die("echo export ARCH=\\\"$arch\\\" >> $home/.bash_profile");
    run_command_or_die("echo export PATH=\\\"$opt_install_dir/bin:\\\$PATH\\\" >> $home/.bash_profile");
    run_command_or_die("echo export LD_LIBRARY_PATH=\\\"$opt_install_dir/lib:\\\$LD_LIBRARY_PATH\\\" >> $home/.bash_profile");
  }
}

#==========================================================================

sub download_file ($) {
  my ($url)=@_;
  return get_command_status("$opt_dwnldprog $url");
}

#==========================================================================

sub remote_file_exists($) {
  my ($file)=@_;
  my $status = "";

  if ($opt_dwnldprog eq $opt_curlprog) {
    $status = run_command("curl -sl $opt_dwnldsrc/ | grep --quiet '$file' && echo 1");
  } elsif ($opt_dwnldprog eq $opt_wgetprog) {
    $status = run_command("wget -q --spider $opt_dwnldsrc/$file && echo 1");
  } else {
    print my_colored("could not check for file.\n",$DEBUGCOLOR);
  }
  if ($status eq "") {
    return 0;
  } else {
    return 1;
  }
}

#==========================================================================

sub download_src ($$) {
  my ($lib,$version)=@_;

  my $file1 = "$lib-$version.tar.gz";
  my $file2 = "$lib-$version.tar.bz2";
  my $status = 0;

  if ( not -e $file1 and not -e $file2 )
  {
    if (remote_file_exists($file1)) {
      $status = download_file("$opt_dwnldsrc/$file1");
    } elsif (remote_file_exists($file2)) {
      $status = download_file("$opt_dwnldsrc/$file2");
    } else {
      print my_colored("File $file1 or $file2 does not exist on server. \n",$OKCOLOR);
      $status = 1;
    }
    print my_colored("Exit Status : $status\n",$OKCOLOR);

    if ( $status )
    {
      die "$args exited with error" unless $status == 0;
    }
  }
  else { print my_colored("file already exists, not retrieving.\n",$OK_COLOR); }
}

#==========================================================================

sub check_curlprog() {
  my $status = run_command("which curl");
  if ($status eq "") {
    return 0;
  } else {
    return 1;
  }
}

#==========================================================================

sub check_wgetprog() {
  my $status = run_command("which wget");
  if ($status eq "") {
    print my_colored("wget is not installed, checking for curl...\n",$DEBUGCOLOR);
    if(check_curlprog()) {
      $opt_dwnldprog = $opt_curlprog;
      print my_colored("curl found, using curl instead of wget\n",$DEBUGCOLOR);
      return 1;
    } else{
      print my_colored("curl and wget not found... install wget manually\n",$DEBUGCOLOR);
      return 0;
    }
  } else {
    return 1;
  }
}

#==========================================================================

sub install_wgetprog() {
  my $lib = "wget";
  my $version = $packages{$lib}[0];
  safe_chdir($opt_tmp_dir);
  download_src($lib,$version);
  untar_src($lib,$version);
  safe_chdir("$opt_tmp_dir/$lib-$version/");
  run_command_or_die("./configure --prefix=$opt_install_dir");
  run_command_or_die("make");
  run_command_or_die("make install");
}

#==========================================================================

sub untar_src ($$) {
  my ($lib,$version)=@_;
  my  $status = get_command_status("tar zxf $lib-$version.tar.gz");
  if ($status) {
    $status = get_command_status("tar jxf $lib-$version.tar.bz2");
    print my_colored("Exit Status : $status\n",$OKCOLOR);
    die "$args exited with error" unless $status == 0;
  }
}

#==========================================================================

sub install_google_perftools ()
{
  my $lib="google-perftools";
  my $version = $packages{$lib}[0];

  print my_colored("Installing $lib-$version\n",$HEADINGCOLOR);

  safe_chdir($opt_tmp_dir);
  download_src($lib,$version);
  unless ($opt_fetchonly) {
    rmtree "$opt_tmp_dir/$lib-$version";
    untar_src($lib,$version);
    safe_chdir("$opt_tmp_dir/$lib-$version/");
    run_command_or_die("./configure --enable-frame-pointers  --prefix=$opt_install_dir");
    run_command_or_die("make");
    run_command_or_die("make install");
  }
}

#==========================================================================

sub install_gnu ($)
{
  my ($lib)=@_;
  my $version = $packages{$lib}[0];

  print my_colored("Installing $lib-$version\n",$HEADINGCOLOR);

  safe_chdir($opt_tmp_dir);
  download_src($lib,$version);
  unless ($opt_fetchonly) {
    rmtree "$opt_tmp_dir/$lib-$version";
    untar_src($lib,$version);
    safe_chdir("$opt_tmp_dir/$lib-$version/");
    run_command_or_die("./configure --prefix=$opt_install_dir");
    run_command_or_die("make");
    run_command_or_die("make install");
  }
}

#==========================================================================

sub install_gcc4() {
  my $lib = "gcc";
  my $version = $packages{"gcc4"}[$vrs];
  print my_colored("Installing $lib-$version\n",$HEADINGCOLOR);

  my $objdir = "$opt_tmp_dir/$lib-$version-obj";

  safe_chdir($opt_tmp_dir);
  download_src($lib,$version);
  unless ($opt_fetchonly) {
    rmtree "$opt_tmp_dir/$lib-$version";
    rmtree "$objdir";
    untar_src($lib,$version);
    safe_chdir("$opt_tmp_dir/$lib-$version/");

    mkpath $objdir;
    safe_chdir($objdir);
    run_command_or_die("$opt_tmp_dir/$lib-$version/configure --prefix=$opt_install_dir --enable-languages=c,c++,fortran");
    run_command_or_die("make");
    run_command_or_die("make install");
  }
}

#==========================================================================

sub install_gcc3() {
  my $lib = "gcc";
  my $version = $packages{"gcc3"}[$vrs];
  print my_colored("Installing $lib-$version\n",$HEADINGCOLOR);

  my $objdir = "$opt_tmp_dir/$lib-$version-obj";

  safe_chdir($opt_tmp_dir);
  download_src($lib,$version);
  unless ($opt_fetchonly) {
    rmtree "$opt_tmp_dir/$lib-$version";
    rmtree "$objdir";
    untar_src($lib,$version);
    safe_chdir("$opt_tmp_dir/$lib-$version/");

    mkpath $objdir;
    safe_chdir($objdir);
    run_command_or_die("$opt_tmp_dir/$lib-$version/configure --prefix=$opt_install_dir --enable-languages=c,c++,f77");
    run_command_or_die("make");
    run_command_or_die("make install");
  }
}

#==========================================================================

sub install_openssl ()
{
  my ($lib)= "openssl";
  my $version = $packages{$lib}[0];

  print my_colored("Installing $lib-$version\n",$HEADINGCOLOR);

  safe_chdir($opt_tmp_dir);
  download_src($lib,$version);
  unless ($opt_fetchonly) {
    rmtree "$opt_tmp_dir/$lib-$version";
    untar_src($lib,$version);
    safe_chdir("$opt_tmp_dir/$lib-$version/");
    run_command_or_die("./config --prefix=$opt_install_dir --shared --openssldir=$opt_install_dir");
    run_command_or_die("make");
    run_command_or_die("make install");
  }
}

#==========================================================================

sub install_curl ()
{
  my ($lib)= "curl";
  my $version = $packages{$lib}[0];

  print my_colored("Installing $lib-$version\n",$HEADINGCOLOR);

  safe_chdir($opt_tmp_dir);
  download_src($lib,$version);
  unless ($opt_fetchonly) {
    rmtree "$opt_tmp_dir/$lib-$version";
    untar_src($lib,$version);
    safe_chdir("$opt_tmp_dir/$lib-$version/");
    run_command_or_die("./configure --prefix=$opt_install_dir --without-ssl --without-libidn --without-gnutls --disable-ipv6 ");
    run_command_or_die("make");
    run_command_or_die("make install");
  }
}

#==========================================================================

sub install_blas()
{
  if (is_mac()) {
    print "Skipping because MacOSX is smarter and already has it ;) \n"
  } else {
    my $lib = "blas";
    my $version = $packages{$lib}[$vrs];
    print my_colored("Installing $lib-$version\n",$HEADINGCOLOR);

    safe_chdir($opt_tmp_dir);
    download_src($lib,$version);
    unless ($opt_fetchonly) {
      rmtree "$opt_tmp_dir/$lib-$version";
      untar_src($lib,$version);
      safe_chdir("$opt_tmp_dir/$lib-$version/");

      # fix Makefile
      my $filename = 'Makefile';
      safe_copy($filename,"$filename.orig");
      open(OUT, ">$filename") or die ("Error opening config file $filename !\n");
      open(IN,  "<$filename.orig") or die ("Error opening config file $filename.orig !\n");
      while (<IN>) {
      chomp;
          s/FFLAGS=(([\w]*)(\S*))*/FFLAGS=$ENV{FFLAGS}/g;
          s/cc -shared/gcc -shared $ENV{LDFLAGS}/g;
          print "$_\n";
          print OUT "$_\n";
      }
      close IN;
      close OUT;


      run_command_or_die("make all");
      safe_copy("libblas.so.$version","$opt_install_dir/lib/libblas.so.$version") or die;
      safe_copy("libblas.a","$opt_install_dir/lib/libblas.a");

      # fix some links
      safe_chdir("$opt_install_dir/lib");
      rm_file("libblas.so");   # if it fails is OK
      rm_file("libblas.so.3"); # if it fails is OK
      run_command_or_die("ln -sf libblas.so.$version libblas.so.3");
      run_command_or_die("ln -sf libblas.so.$version libblas.so");
    }
  }
}

#==========================================================================

sub install_lapack() {
  if (is_mac()) {
    print "Skipping because MacOSX is smarter and already has it ;) \n"
  } else {
    my $lib = "lapack";
    my $version = $packages{$lib}[$vrs];
    print my_colored("Installing $lib-$version\n",$HEADINGCOLOR);

    safe_chdir($opt_tmp_dir);
    download_src($lib,$version);

    unless ($opt_fetchonly) {
      rmtree "$opt_tmp_dir/$lib-$version";
      untar_src($lib,$version);
      safe_chdir("$opt_tmp_dir/$lib-$version/SRC");

      # fix Makefile
      my $filename = 'Makefile';
      safe_copy($filename,"$filename.orig");
      open(OUT, ">$filename") or die ("Error opening config file $filename !\n");
      open(IN,  "<$filename.orig") or die ("Error opening config file $filename.orig !\n");
      while (<IN>) {
      chomp;
          s/FFLAGS=(([\w]*)(\S*))*/FFLAGS=$1 $ENV{FFLAGS}/g;
          s/BLAS_PATH=/BLAS_PATH=$opt_install_dir\/lib/g;
          s/INSTALL_PATH=/INSTALL_PATH=$opt_install_dir\/lib/g;
          s/cc -shared/gcc -shared $ENV{LDFLAGS}/g;
          print "$_\n";
          print OUT "$_\n";
      }
      close IN;
      close OUT;

      run_command_or_die("make all");

      safe_copy("liblapack.so.$version","$opt_install_dir/lib/liblapack.so.$version") or die;
      safe_copy("liblapack.a","$opt_install_dir/lib/liblapack.a");

      # fix some links
      safe_chdir("$opt_install_dir/lib");
      rm_file("liblapack.so");   # if it fails is OK
      rm_file("liblapack.so.3"); # if it fails is OK
      run_command_or_die("ln -sf liblapack.so.$version liblapack.so.3");
      run_command_or_die("ln -sf liblapack.so.$version liblapack.so");
    }
  }
}

#==========================================================================

sub install_openmpi() {
  my $lib = "openmpi";
  my $version = $packages{$lib}[$vrs];
  print my_colored("Installing $lib-$version\n",$HEADINGCOLOR);

  safe_chdir($opt_tmp_dir);
  download_src($lib,$version);
  unless ($opt_fetchonly)
  {
    rmtree "$opt_tmp_dir/$lib-$version";
    untar_src($lib,$version);
    safe_chdir("$opt_tmp_dir/$lib-$version/");
    run_command_or_die("./configure --enable-shared --enable-static --with-threads=posix --disable-mpi-f90  --prefix=$opt_mpi_dir");
    run_command_or_die("make");
    run_command_or_die("make install");
  }
}

#==========================================================================

sub install_mpich2() {
  my $lib = "mpich2";
  my $version = $packages{$lib}[$vrs];
  print my_colored("Installing $lib-$version\n",$HEADINGCOLOR);

  safe_chdir($opt_tmp_dir);
  download_src($lib,$version);
  unless ($opt_fetchonly)
  {
      rmtree "$opt_tmp_dir/$lib-$version";
      untar_src($lib,$version);
      safe_chdir("$opt_tmp_dir/$lib-$version/");

      if ((exists $ENV{F90}) or (exists $ENV{F90FLAGS}))
      {
        $ENV{F90} = "";
        $ENV{F90FLAGS} = "";
        print "Unsetting F90 and F90FLAGS";
      }

      run_command_or_die("./configure --disable-fc --enable-fast=all --enable-shared --prefix=$opt_mpi_dir");
      run_command_or_die("make");
      run_command_or_die("make install");
  }
}

#==========================================================================

sub install_mvapich2() {
  my $lib = "mvapich2";
  my $version = $packages{$lib}[$vrs];
  print my_colored("Installing $lib-$version\n",$HEADINGCOLOR);

  safe_chdir($opt_tmp_dir);
  download_src($lib,$version);
  unless ($opt_fetchonly)
  {
      rmtree "$opt_tmp_dir/$lib-$version";
      untar_src($lib,$version);
      safe_chdir("$opt_tmp_dir/$lib-$version/");

      if ((exists $ENV{F90}) or (exists $ENV{F90FLAGS}))
      {
        $ENV{F90} = "";
        $ENV{F90FLAGS} = "";
        print "Unsetting F90 and F90FLAGS";
      }

#      run_command_or_die("./configure --with-device=ch3:psm --enable-shared --prefix=$opt_mpi_dir");
      run_command_or_die("./configure --prefix=$opt_mpi_dir --with-pmi=slurm --with-pm=none --enable-romio --enable-f77 --enable-fc --enable-cxx --enable-shared");
      run_command_or_die("make");
      run_command_or_die("make install");
  }
}

#==========================================================================

sub install_cgnslib() {
  my $lib = "cgnslib";
  my $version = $packages{$lib}[$vrs];
  print my_colored("Installing $lib-$version\n",$HEADINGCOLOR);

  safe_chdir($opt_tmp_dir);
  download_src($lib,$version);
  unless ($opt_fetchonly)
  {
    rmtree "$opt_tmp_dir/$lib-$version";
    untar_src($lib,$version);
    safe_chdir("$opt_tmp_dir/$lib-$version/");

    my $loptions = "--enable-gcc --enable-lfs --without-fortran --enable-shared=all --with-zlib --prefix=$opt_install_dir";
    if ($arch eq "x86_64" )
    {
        $loptions  = $loptions . " --enable-64bit";
    }
    run_command_or_die("./configure " . $loptions);
    run_command_or_die("make");
    run_command_or_die("make install");
  }
}

#==========================================================================

sub install_cgnstools() {
  my $lib = "cgnstools";
  my $version = $packages{$lib}[$vrs];
  my $cgnslib = "cgnslib";
  my $cgnsversion = $packages{$cgnslib}[$vrs];
  print my_colored("Installing $lib-$version\n",$HEADINGCOLOR);

  # first install a cgnlslib distribution dir
  safe_chdir($opt_tmp_dir);
  download_src($lib,$version);
  unless ($opt_fetchonly)
  {
    run_command_or_die("tar -zxf  $cgnslib-$cgnsversion.tar.gz -C $opt_install_dir/local");
    safe_chdir("$opt_install_dir/local/$cgnslib-$cgnsversion/");

    my $loptions = "--enable-gcc --enable-lfs --without-fortran --enable-shared=all --with-zlib";
    if ($arch eq "x86_64" )
    {
        $loptions  = $loptions . " --enable-64bit";
    }
    run_command_or_die("./configure " . $loptions);
    run_command_or_die("make");
  }

  safe_chdir($opt_tmp_dir);
  download_src($lib,$version);
  unless ($opt_fetchonly)
  {
    rmtree "$opt_tmp_dir/$lib-$version";
    untar_src($lib,$version);
    safe_chdir("$opt_tmp_dir/$lib-$version/");

    my $loptions = "--enable-gcc --with-cgns=$opt_install_dir/local/$cgnslib-$cgnsversion/ --prefix=$opt_install_dir";
    if ($arch eq "x86_64" )
    {
        $loptions  = $loptions . " --enable-64bit";
    }
    run_command_or_die("./configure " . $loptions);
    run_command_or_die("make");
    run_command_or_die("make install");
}
}

#==========================================================================

sub install_mpich() {
  my $lib = "mpich";
  my $version = $packages{$lib}[$vrs];
  print my_colored("Installing $lib-$version\n",$HEADINGCOLOR);

  safe_chdir($opt_tmp_dir);
  download_src($lib,$version);
  unless ($opt_fetchonly)
  {
      rmtree "$opt_tmp_dir/$lib-$version";
      untar_src($lib,$version);
      safe_chdir("$opt_tmp_dir/$lib-$version/");
      
      my $optim_options = "--enable-g=dbg,mem,log";
      if ($opt_debug eq "0") {  
        $optim_options = "--enable-fast=all,O3";
      }

      if ((exists $ENV{F90}) or (exists $ENV{F90FLAGS}))
      {
        $ENV{F90} = "";
        $ENV{F90FLAGS} = "";
        print "Unsetting F90 and F90FLAGS";
      }
 
      # --enable-shared
      run_command_or_die("./configure --prefix=$opt_mpi_dir --disable-f77 --disable-fortran --disable-fc $optim_options --with-device=ch3:nemesis"); 
      run_command_or_die("make");
      run_command_or_die("make install");
  }
}

#==========================================================================

sub install_parmetis () {
  my $lib = "parmetis";
  my $version = $packages{$lib}[$vrs];
 
  if ( $version eq "3.1.1") {
    install_parmetis3 ();
  }
  else {
   print my_colored("Installing $lib-$version\n",$HEADINGCOLOR);

  if ($opt_parmetis_dir eq "") {
     $opt_parmetis_dir = "$opt_install_mpi_dir";
  }

  my $mpicc_comp  = "mpicc";
  my $mpicxx_comp = "mpicxx";
  if ($opt_cray) {
     $mpicc_comp  = "cc";
     $mpicxx_comp = "CC";
  }

  #
  if ($opt_ibmgnu) {
     $mpicc_comp  = "mpigcc";
     $mpicxx_comp = "mpig++";
  }

  my $sharedflag = "shared=1";
  if ( $opt_static eq "1" )
  {
    $sharedflag = "";
  }
  safe_chdir($opt_tmp_dir);
  download_src("parmetis",$version);
  unless ($opt_fetchonly) {
    rmtree "$opt_tmp_dir/parmetis-$version";
    untar_src("parmetis",$version);
    # build metis
    if ( $version eq "4.0.3") {
       safe_chdir("$opt_tmp_dir/parmetis-$version/metis");
       run_command_or_die("make config $sharedflag prefix=$opt_parmetis_dir cc=$opt_mpi_dir/bin/$mpicc_comp cxx=$opt_mpi_dir/bin/$mpicxx_comp");
       run_command_or_die("make install");
    }
    #build parmetis
    safe_chdir("$opt_tmp_dir/parmetis-$version/");
    run_command_or_die("make config $sharedflag prefix=$opt_parmetis_dir cc=$opt_mpi_dir/bin/$mpicc_comp cxx=$opt_mpi_dir/bin/$mpicxx_comp");
    run_command_or_die("make install");
  }  
 }
}

#==========================================================================

sub install_parmetis3 () {
  my $lib = "parmetis";
  my $version = $packages{$lib}[$vrs];
  print my_colored("Installing $lib-$version\n",$HEADINGCOLOR);

  if ($opt_parmetis_dir eq "") {
     $opt_parmetis_dir = "$opt_install_mpi_dir";
  }

  my $include_dir = "$opt_parmetis_dir/include/";
  my $lib_dir = "$opt_parmetis_dir/lib/";

  mkpath $include_dir;
  mkpath $lib_dir;

  safe_chdir($opt_tmp_dir);
  download_src("ParMetis",$version);
  unless ($opt_fetchonly) {

    rmtree "$opt_tmp_dir/ParMetis-$version";
    untar_src("ParMetis",$version);
    safe_chdir("$opt_tmp_dir/ParMetis-$version/");

    if (is_mac()) { # add include for malloc.h
       my $filename = 'Makefile.in';
        safe_copy($filename,"$filename.orig");
        open(OUT, ">$filename") or die ("Error opening config file $filename !\n");
        open(IN,  "<$filename.orig") or die ("Error opening config file $filename.orig !\n");
        while (<IN>) {
        chomp;
        s/(^INCDIR\s=\s?$)/INCDIR = -I\/usr\/include\/malloc\//g;
        #print "$_\n";
        print OUT "$_\n";
        }
        print my_colored("Modified Makefile.in to include malloc for MacOSX\n",$DEBUGCOLOR);
        close IN;
        close OUT;
    }

    if ($opt_nompi) { # substitute mpicc for gcc
        my $filename = 'Makefile.in';
        safe_copy($filename,"$filename.orig");
        open(OUT, ">$filename") or die ("Error opening config file $filename !\n");
        open(IN,  "<$filename.orig") or die ("Error opening config file $filename.orig !\n");
        while (<IN>) {
        chomp;
        s/mpicc/gcc/g;
        print "$_\n";
        print OUT "$_\n";
        }
        close IN;
        close OUT;
    }

    safe_chdir("METISLib");
    run_command_or_die("make $opt_makeopts");
    safe_chdir("..");

    safe_chdir("ParMETISLib");
    run_command_or_die("make $opt_makeopts");
    safe_chdir("..");

    safe_copy("parmetis.h","$include_dir/parmetis.h");
    safe_copy("libmetis.a","$lib_dir/libmetis.a");
    safe_copy("libparmetis.a","$lib_dir/libparmetis.a");
  }
}

#==========================================================================

sub install_petsc ()
{
  my $lib = "petsc";
  my $version = $packages{"$lib"}[$vrs];
  my $source_file = "$lib-$version.tar.gz";
  my $fblas_name = "fblaslapack-3.1.1.tar.gz";
  my $fblas_file = "$opt_tmp_dir/$fblas_name";
  
  print my_colored("Installing $lib-$version\n",$HEADINGCOLOR);

  safe_chdir($opt_tmp_dir);

  if ( not -e $source_file ) { download_file("$opt_dwnldsrc/$source_file") };
  if ( not -e $fblas_file  ) { download_file("$opt_dwnldsrc/$fblas_name") };

  unless ($opt_fetchonly)
  {
    my $build_dir   = "$opt_tmp_dir/$lib-$version";

    my $install_dir = "$opt_install_mpi_dir/petsc";
    if ($opt_petsc_dir) {
      # set the installation directory to the one specified by the user if any 
      $install_dir = "$opt_petsc_dir";
    }

    my $petsc_arch  = "arch-$arch";
    if (is_mac()) { $petsc_arch = "arch-darwin"; };

    $ENV{PETSC_DIR}  = "$build_dir";
    $ENV{PETSC_ARCH} = $petsc_arch;

    # extract sources to build dir
    rmtree $build_dir;
    untar_src($lib,$version);

    safe_chdir("$build_dir");

    my $wdebug = "--with-debugging=1";
    if ($opt_debug eq "0") { $wdebug = "--with-debugging=0" };

    my $use_int64 = "";
    if ($opt_int64 eq 1) {
       $use_int64 = "--with-64-bit-indices";
    } 
   
    my $wviennacl = "";
    my $cuda_support = "";
    if ($opt_cuda_dir) {
      print my_colored("Installing PETSc with CUDA support\n", $HEADINGCOLOR);
      # set paths for CUDA bin and libraries
      $install_dir = "$opt_install_mpi_dir/petsc_cuda";
       if ($opt_petsc_dir) {
         # set the installation directory to the one specified by the user if any 
         $install_dir = "$opt_petsc_dir";
       }

      $opt_cuda_bin = "$opt_cuda_dir/bin/nvcc";
      $opt_with_cuda = 1;
      run_command_or_die("export PATH=$opt_cuda_dir/bin:\\\$PATH");
      run_command_or_die("export LD_LIBRARY_PATH=/usr/lib64:$opt_cuda_dir/lib64:/usr/lib64/nvidia:\\\$LD_LIBRARY_PATH");
     
      if ($version eq "3.6.3" or $version eq "3.7.3a" or $version eq "3.7.3-next" or $version eq "3.7.6" or $version eq "3.9.0" or $version eq "3.9.3" or $version eq "3.12.3") {
         my $cuda_arch="--with-cuda-arch=sm_30";
	 my $cusp_dir="cusp-v0.4.0";  
	 if ($version eq "3.9.0" or $version eq "3.9.3" or $version eq "3.12.3") {
           $cuda_arch="--CUDAFLAGS=-arch=sm_30";
           $cusp_dir="cusplibrary-0.5.0";
         }
	 
	 if ($opt_viennacl eq "0") {
           # download and unpack compatible cusp source files
           my $cusp_name = "$cusp_dir.zip";
           my $cusp_file = "$opt_tmp_dir/$cusp_name";
           if ( not -e $cusp_file ) { download_file("$opt_dwnldsrc/$cusp_name") }; 
           run_command_or_die("cp $cusp_file .");
           run_command_or_die("cp $cusp_name $opt_install_dir/ ; cd $opt_install_dir ; unzip $cusp_name ; rm -f $cusp_name ; cd -");
        
	   my $cuda_arch="--with-cuda-arch=sm_30";
	   if ($version eq "3.9.0" or $version eq "3.9.3" or $version eq "3.12.3") {$cuda_arch="--CUDAFLAGS=-arch=sm_30";} 

         # gory fixes to fix inclusions paths for cups and thrust     
	 run_command_or_die("cp -R $opt_install_dir/$cusp_dir/cusp $opt_install_dir/include");
         run_command_or_die("cp -R $opt_cuda_dir/include/thrust $opt_install_dir/include/cusp");
         $cuda_support = "--with-cudac=$opt_cuda_bin --with-cuda-dir=$opt_cuda_dir --with-cuda=1 --with-cusp=1 --with-cusp-dir=$opt_install_dir/$cusp_dir --with-thrust=1 --with-thrust-dir=$opt_cuda_dir --with-precision=double --with-clanguage=c $cuda_arch";
        }
       else  {
        $wviennacl = "--download-viennacl --with-opencl-include=/usr/local/cuda-6.5/targets/x86_64-linux/include --with-opencl-lib=/usr/local/cuda-6.5/targets/x86_64-linux/lib/libOpenCL.so --with-openmp=1";
        $cuda_support = "--with-cudac=$opt_cuda_bin --with-cuda-dir=$opt_cuda_dir --with-cuda=1 --with-precision=double --with-clanguage=c $cuda_arch";
       } 
     } 
      else {

      # download and unpack compatible thrust source files
      my $thrust_name = "thrust-1.4.0.zip";
      my $thrust_file = "$opt_tmp_dir/$thrust_name";
      if ( not -e $thrust_file ) { download_file("$opt_dwnldsrc/$thrust_name") };
      run_command_or_die("cp $thrust_name $opt_install_dir/ ; cd $opt_install_dir ; unzip $thrust_name ; rm -f $thrust_name ; cd -");
      
      # download and unpack compatible cusp source files
      my $cusp_name = "cusp-v0.2.0.zip";
      my $cusp_file = "$opt_tmp_dir/$cusp_name";
      if ( not -e $cusp_file ) { download_file("$opt_dwnldsrc/$cusp_name") };
      run_command_or_die("cp $cusp_name $opt_install_dir/ ; cd $opt_install_dir ; unzip $cusp_name ; rm -f $cusp_name ; cd -");
   
      # configuration options for enabling GPU support 
      $cuda_support = "--with-cudac=$opt_cuda_bin --with-cuda-dir=$opt_cuda_dir --with-cuda=1 --with-cusp=1 --with-thrust=1 --with-thrust-dir=$opt_install_dir --with-cusp-dir=$opt_install_dir"
      }
    }
    
    # update the petsc installation directory
    $opt_petsc_dir = "$install_dir";
  
    my $dynamicload = "--with-dynamic-loading=1"; 
    my $fcflag = "0"; 
    my $wblaslib = "";
    #if (is_mac()) { 
      # use built-in optimized blas-lapack lib
    #  $wblaslib = "--with-blas-lapack-lib=\"-framework vecLib\"";
    #} else { 
      # use the downloaded blas sources
    #  $wblaslib = "--download-f-blas-lapack=1";
    if ( is_mac() and $opt_static eq "1") {
     #$wblaslib = "--with-blas-lapack-dir=$opt_blaslapack_dir";
     #$wblaslib  = "--with-blas-lib=$opt_blaslapack_dir/libblas.a --with-lapack-lib=$opt_blaslapack_dir/liblapack.a"
     # $wblaslib = "--download-f2blaslapack=1";
     #$fcflag = "1";
     $wblaslib = "--download-f2cblaslapack=1";
     # removed duplicated symbols that could create linking problems for static building 
     run_command_or_die("ar d $opt_petsc_dir/lib/libf2clapack.a xerbla_array.o");
     run_command_or_die("ar d $opt_petsc_dir/lib/libf2clapack.a xerbla.o");	
    }   
    else {
     $wblaslib = "--download-f2cblaslapack=1";
    }
    
  if ($version eq "3.2-p6") {
   run_command_or_die("./configure --prefix=$install_dir $wdebug COPTFLAGS='-O3 ' FOPTFLAGS='-O3 ' --with-mpi-dir=$opt_mpi_dir $wblaslib --with-fortran=$fcflag --with-shared-libraries=1 $dynamicload --with-c++-support $cuda_support --PETSC_ARCH=$petsc_arch");
   run_command_or_die("make $opt_makeopts");
  }
  elsif ($version eq "3.6.3" or $version eq "3.7.3a" or $version eq "3.7.3-next" or $version eq "3.7.6" or $version eq "3.9.0" or $version eq "3.9.3" or $version eq "3.12.3") {
    my $FFLAGS = "-O3 ";
    my $F90FLAGS ="-O3 ";
    my $CFLAGS ="-O3 "; 
    my $CXXFLAGS ="-O3 ";
    my $FOPTFLAGS ="-O3 ";
    my $COPTFLAGS ="-O3 ";
    my $CXXOPTFLAGS ="-O3 ";

    # IBM BG/Q
    if (not $opt_ibmgnu and $petsc_arch eq "arch-ppc64" ) {
       $CXXFLAGS="-O3 -qstrict -qarch=qp -qtune=qp";
       $CFLAGS="-O3 -qstrict -qarch=qp -qtune=qp";
       $FFLAGS="-O3 -qstrict -qarch=qp -qtune=qp";
       $F90FLAGS="-O3 -qstrict -qarch=qp -qtune=qp";
       $CXXOPTFLAGS="-O3 -qstrict -qarch=qp -qtune=qp"; 
       $COPTFLAGS="-O3 -qstrict -qarch=qp -qtune=qp";
       $FOPTFLAGS="-O3 -qstrict -qarch=qp -qtune=qp";
    }

    if ($version eq "3.7.6" or $version eq "3.9.0" or $version eq "3.9.3" or $version eq "3.12.3") {
      $dynamicload ="";
    } 
    my $link_libraries = "--with-shared-libraries=1"; 
    if ( $opt_static eq "1" ) 
    {
      $link_libraries = "--with-shared-libraries=0";
    } 

   if ($wdebug eq "--with-debugging=0") { 
    run_command_or_die("./configure --prefix=$install_dir $wdebug $use_int64 FFLAGS=$FFLAGS F90FLAGS=$F90FLAGS CFLAGS=$CFLAGS CXXFLAGS=$CXXFLAGS CXXOPTFLAGS=$CXXOPTFLAGS COPTFLAGS=$COPTFLAGS FOPTFLAGS=$FOPTFLAGS --with-mpi-dir=$opt_mpi_dir $wblaslib --with-fortran=$fcflag $link_libraries $wviennacl $cuda_support --PETSC_ARCH=$petsc_arch");
   }
   else {
     run_command_or_die("./configure --prefix=$install_dir $wdebug $use_int64 --with-mpi-dir=$opt_mpi_dir $wblaslib --with-fortran=$fcflag $link_libraries $cuda_support --PETSC_ARCH=$petsc_arch");
   }
   run_command_or_die("make");
  }
  else {
   run_command_or_die("./configure --prefix=$install_dir $wdebug $use_int64 COPTFLAGS='-O3 ' FOPTFLAGS='-O3 ' --with-mpi-dir=$opt_mpi_dir $wblaslib --with-fortran=$fcflag $link_libraries $dynamicload $cuda_support --PETSC_ARCH=$petsc_arch");
   run_command_or_die("make");
  }
    run_command_or_die("make install PETSC_DIR=$build_dir");
  }
}

#==========================================================================

sub install_trilinos() {
  my $lib = "trilinos";
  my $version = $packages{"trilinos"}[$vrs];
  print my_colored("Installing $lib-$version\n",$HEADINGCOLOR);
  print my_colored("In trilinos subroutine..\n",$HEADINGCOLOR);


  safe_chdir($opt_tmp_dir);
  download_src($lib,$version);
  mkdir "Tribuild";

  my $mpiopt;
  unless ($opt_nompi) {
      $mpiopt = " -D TPL_ENABLE_MPI:BOOL=ON \\
-D MPI_BASE_DIR_PATH:PATH=$opt_mpi_dir \\
-D CMAKE_C_COMPILER:FILEPATH=$opt_mpi_dir/bin/mpicc \\
-D CMAKE_CXX_COMPILER:FILEPATH=$opt_mpi_dir/bin/mpic++ \\
-D CMAKE_Fortran_COMPILER:FILEPATH=$opt_mpi_dir/bin/mpif77 " 
 }


  unless ($opt_fetchonly) {
    rmtree "$opt_tmp_dir/$lib-$version";
    untar_src($lib,$version);
    safe_chdir("$opt_tmp_dir/Tribuild/");


    run_command_or_die("$opt_cmake_dir/bin/cmake \\
-D CMAKE_BUILD_TYPE:STRING=RELEASE \\
-D Trilinos_ENABLE_DEFAULT_PACKAGES:BOOL=OFF \\
-D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=ON \\
-D Trilinos_ENABLE_TESTS:BOOL=OFF \\
-D Trilinos_ENABLE_Amesos:BOOL=ON \\
-D Trilinos_ENABLE_AztecOO:BOOL=ON \\
-D Trilinos_ENABLE_Belos:BOOL=ON \\
-D Trilinos_ENABLE_Didasko:BOOL=OFF \\
-D Didasko_ENABLE_TESTS:BOOL=OFF \\
-D Didasko_ENABLE_EXAMPLES:BOOL=OFF \\
-D Trilinos_ENABLE_Epetra:BOOL=ON \\
-D Trilinos_ENABLE_EpetraExt:BOOL=ON \\
-D Trilinos_ENABLE_Ifpack:BOOL=ON \\
-D Trilinos_ENABLE_Meros:BOOL=ON \\
-D Trilinos_ENABLE_ML:BOOL=ON \\
-D Trilinos_ENABLE_RTOp:BOOL=ON \\
-D Trilinos_ENABLE_Teuchos:BOOL=ON \\
-D Trilinos_ENABLE_Thyra:BOOL=ON \\
-D Trilinos_ENABLE_Triutils:BOOL=ON \\
-D Trilinos_ENABLE_Stratimikos:BOOL=ON \\
-D Trilinos_ENABLE_Zoltan:BOOL=OFF \\
-D TPL_ENABLE_BLAS:BOOL=ON \\
-D TPL_ENABLE_LAPACK:BOOL=ON \\
$mpiopt \\
-D CMAKE_VERBOSE_MAKEFILE:BOOL=FALSE \\
-D BUILD_SHARED_LIBS:BOOL=ON\\
-D Trilinos_VERBOSE_CONFIGURE:BOOL=FALSE \\
-D CMAKE_INSTALL_PREFIX:PATH=$opt_install_mpi_dir \\
$opt_tmp_dir/$lib-$version-Source"
);

    run_command_or_die("make");
    run_command_or_die("make install");
  }
}

#==========================================================================

sub install_thrust ()
{
  my ($lib)= "thrust";
  my $version = $packages{$lib}[0];

  print my_colored("Installing $lib-$version\n",$HEADINGCOLOR);

  safe_chdir($opt_tmp_dir);
  download_src($lib,$version);
  unless ($opt_fetchonly) {
    rmtree "$opt_tmp_dir/$lib-$version";
    untar_src($lib,$version);
    safe_chdir("$opt_tmp_dir/$lib-$version/");
    run_command_or_die("cp -R $opt_tmp_dir/$lib-$version/ $opt_install_dir/thrust");
  }
}

#==========================================================================

sub install_cusp ()
{
  my ($lib)= "cusp";
  my $version = $packages{$lib}[0];

  print my_colored("Installing $lib-$version\n",$HEADINGCOLOR);

  safe_chdir($opt_tmp_dir);
  download_src($lib,$version);
  unless ($opt_fetchonly) {
    rmtree "$opt_tmp_dir/$lib-$version";
    untar_src($lib,$version);
    safe_chdir("$opt_tmp_dir/$lib-$version/");
    run_command_or_die("cp -R $opt_tmp_dir/$lib-$version/ $opt_install_dir/cusp");
  }
}

#==========================================================================

sub install_subversion() {
  my $lib = "subversion";
  my $version = $packages{$lib}[$vrs];
  print my_colored("Installing $lib-$version\n",$HEADINGCOLOR);

  safe_chdir($opt_tmp_dir);
  download_file("$opt_dwnldsrc/$lib-$version-linux-x86.tar.gz");
  unless ($opt_fetchonly) {
  run_command_or_die("tar zxvf $lib-$version-linux-x86.tar.gz -C $opt_install_dir");
  }
}

#==========================================================================

sub install_libfaketime() {
  my $lib = "libfaketime";
  my $version = $packages{$lib}[$vrs];
  my $pack = "$lib-$version";
  print my_colored("Installing $pack\n",$HEADINGCOLOR);

  safe_chdir($opt_tmp_dir);
  if ( not -e "$pack.tar.gz" ) { download_file("$opt_dwnldsrc/$pack.tar.gz"); }

  unless ($opt_fetchonly)
  {
    rmtree "$opt_tmp_dir/$pack";
    run_command_or_die("tar zxf $pack.tar.gz");
    safe_chdir("$opt_tmp_dir/$pack/");

    run_command_or_die("make");

    safe_copy("libfaketime.so.1","$opt_install_dir/lib/libfaketime.so.1") or die;
    safe_copy("libfaketimeMT.so.1","$opt_install_dir/lib/libfaketimeMT.so.1") or die;
  }
}
#==========================================================================

sub install_boost()
{
  my $lib = "boost";
  my $version = $packages{$lib}[$vrs];
  my $pack = "$lib\_$version";
  print my_colored("Installing $pack\n",$HEADINGCOLOR);

  safe_chdir($opt_tmp_dir);
  # if ( not -e "$pack.tar.bz2" ) { download_file("$opt_dwnldsrc/$pack.tar.bz2"); }

  if ( not -e "$pack.tar.bz2" ) {
    my $mversion = $version;
    $mversion =~ s/_/./g;
    print my_colored("Downloading $pack\n with version $mversion",$HEADINGCOLOR);
    run_command_or_die("wget https://sourceforge.net/projects/boost/files/boost/$mversion/$pack.tar.bz2");
  }

  unless ($opt_fetchonly)
  {
    rmtree "$opt_tmp_dir/$pack";
    run_command_or_die("tar jxf $lib\_$version.tar.bz2");
    safe_chdir("$opt_tmp_dir/$pack/");

    if ($version  eq "1_42_0")
      {
	# build toolset
	safe_chdir("tools/jam/src");
      }
    
    my $toolset = "gcc";
    if($ENV{CC} eq "icc") { $toolset = "intel-linux"; }
    
    # in case g++ is speaical path
    $ENV{GCC} = $ENV{CC};
    $ENV{GXX} = $ENV{CXX};


    my $boost_arch;
    if($arch eq "x86_64") { $boost_arch = "linuxx86_64" ;  }
    if($arch eq "i686")   { $boost_arch = "linuxx86" ;  }

    if(is_mac())         
    { 
	  $toolset = "darwin";
      $boost_arch = "macosxx86"; 
      
      # If Snow Leopard
      my $capable64 = run_command("sysctl hw | grep 'hw.cpu64bit_capable: [0-9]'");
      my $OSversion = run_command("sw_vers | grep 'ProductVersion:'");
      
      # AL: change here, removed check on version
      # if ($capable64 =~ /hw.cpu64bit_capable:\s1/ && $OSversion =~ /10\.6\.*/) 
      if ($capable64 =~ /hw.cpu64bit_capable:\s1/)
      {
          $boost_arch = "macosxx86_64";    
      }
    }

    if ($version  eq "1_42_0")
      {
	# disable compression filters in boost because some systems like ubuntu
	# dont have the zlib-dev installed by default
	$ENV{NO_COMPRESSION} = "1";
	
	run_command_or_die("sh build.sh $toolset");
	
	# build boost libs
	safe_chdir("../../..");
      }    
    
    my $boostmpiopt=" --without-mpi ";
    unless ($opt_nompi or $version  eq "1_59_0" or $version eq "1_62_0" or $version eq "1_66_0" or $version eq "1_70_0" or $version eq "1_72_0" or $version eq "1_76_0" or $version eq "1_79_0")  {
      $boostmpiopt=" --with-mpi cxxflags=-DBOOST_MPI_HOMOGENEOUS ";
      open  (USERCONFIGJAM, ">>./tools/build/v2/user-config.jam") || die("Cannot Open File ./tools/build/v2/user-config.jam") ;
      print  USERCONFIGJAM <<ZZZ;


# ----------------------
# mpi configuration.
# ----------------------
using mpi : $opt_mpi_dir/bin/mpicxx ;

ZZZ
      close (USERCONFIGJAM); 
    }
    
    if ($version  eq "1_42_0") 
      {
	run_command_or_die("./tools/jam/src/bin.$boost_arch/bjam --prefix=$opt_install_dir $opt_makeopts --with-test --with-thread --with-iostreams --with-filesystem --with-system --with-regex --with-date_time --with-program_options $boostmpiopt toolset=$toolset threading=multi variant=release stage install");
      }
    
    my $static_link = "";
    if ( $opt_static eq "1" ) 
    {
       $static_link = " link=static";
    }
    if ($version  eq "1_54_0" or $version  eq "1_59_0" or $version  eq "1_62_0" or $version eq "1_66_0" or $version eq "1_70_0" or $version eq "1_72_0" or $version eq "1_76_0") 
      {
	run_command_or_die("./bootstrap.sh --prefix=$opt_install_dir -with-libraries=test,thread,iostreams,filesystem,system,regex,date_time toolset=$toolset threading=multi variant=release stage");
	run_command_or_die("./b2 $static_link install");
      }
   if ($version eq "1_79_0")
      {
        run_command_or_die("./bootstrap.sh --prefix=$opt_install_dir -with-libraries=test,thread,iostreams,filesystem,system,regex,date_time,atomic toolset=$toolset threading=multi variant=release stage");
        run_command_or_die("./b2 $static_link install");
      } 
  }
}

#==========================================================================

sub install_cmake() {
  my $lib = "cmake";
  my $version = $packages{$lib}[$vrs];
  my $pack = "$lib-$version";
  print my_colored("Installing $pack\n",$HEADINGCOLOR);

  safe_chdir($opt_tmp_dir);
  download_src($lib,$version);

  unless ($opt_fetchonly) {

    rmtree "$opt_tmp_dir/$pack";
    run_command_or_die("tar zxf $pack.tar.gz");
    safe_chdir("$opt_tmp_dir/$pack/");

    run_command_or_die("./bootstrap --prefix=$opt_cmake_dir");
    run_command_or_die("make");
    run_command_or_die("make install");

  }
}

#==========================================================================

sub install_hdf5() {
  my $lib = "hdf5";
  my $version = $packages{$lib}[$vrs];
  print my_colored("Installing $lib-$version\n",$HEADINGCOLOR);

  safe_chdir($opt_tmp_dir);
  download_src($lib,$version);
  unless ($opt_fetchonly) {
    rmtree "$opt_tmp_dir/$lib-$version";
    untar_src($lib,$version);
    safe_chdir("$opt_tmp_dir/$lib-$version/");

    my $old_cc  = $ENV{CC};
    my $old_cxx = $ENV{CXX};
    my $mpiopt;
    unless ($opt_nompi) {
        $ENV{CC}   = "mpicc";
        $ENV{CXX}  = "mpic++";
        $mpiopt = "--enable-parallel";
    }

    run_command_or_die("./configure --prefix=$opt_install_mpi_dir $mpiopt --enable-cxx --enable-zlib --enable-linux-lfs --with-gnu-ld --enable-hl");
    run_command_or_die("make");
    run_command_or_die("make install");

    $ENV{CC}   = $old_cc;
    $ENV{CXX}  = $old_cxx;
  }
}

#==========================================================================

sub print_info() # print information about the
{
    print my_colored("Installing COOLFluiD and its dependencies\n",$HEADINGCOLOR);

    print_var("Install     dir ","$opt_install_dir");
    print_var("Install MPI dir ","$opt_install_mpi_dir");
    print_var("MPI         dir ","$opt_mpi_dir");
    print_var("Temporary   dir ","$opt_tmp_dir");

# Env vars
    print_var(PATH,$ENV{PATH});
    print_var(LD_LIBRARY_PATH,$ENV{LD_LIBRARY_PATH});
    print_var(CC,$ENV{CC});
    print_var(CXX,$ENV{CXX});
    print_var(FC,$ENV{FC});
    print_var(CFLAGS,$ENV{CFLAGS});
    print_var(CXXFLAGS,$ENV{CXXFLAGS});
    print_var(FFLAGS,$ENV{FFLAGS});
    print_var(F77FLAGS,$ENV{F77FLAGS});
    print_var(F90FLAGS,$ENV{F90FLAGS});

# Options
#     while ( my ($key, $value) = each(%options) ) {
#         print_var($key,get_option($key));
#     }

# User prefs
#     while ( my ($key, $value) = each(%user_pref) ) {
#         print_var($key,$value);
#     }
}

#==========================================================================

sub set_install_all()
{
  foreach $pname (keys %packages) {
    $packages{$pname}[$ins] = $packages{$pname}[$dft];
  }
}

#==========================================================================

sub install_packages()
{
  print_info();
  check_wgetprog();

    # if 'all' exists, copy the [$dft] to [$ins]
    for ($i=0; $i < scalar @opt_install; $i++)
    {
        if ($opt_install[$i] eq 'all') { set_install_all(); }
    }

    # if there is no package selected, then also copy the [$dft] to [$ins]
    if (scalar @opt_install == 0) { set_install_all(); }

    # turn on the manually selected packages
    for ($i=0; $i < scalar @opt_install; $i++)
    {
        my $opt = $opt_install[$i];
        if (exists $packages{$opt})
        {
            $packages{$opt}[$ins] = 'on';
        }
        elsif (!($opt eq 'all')) {
            print my_colored("Package does not exist: $opt\n",$ERRORCOLOR);
        }
    }

    my %install_packages = ();

    # sort the packages to install by priority
    foreach $pname (keys %packages) {
    #       print "$pname\n";
        if ($packages{$pname}[$ins] eq 'on') {
            $install_packages{$packages{$pname}[$pri]} = $pname;
        }
    }

    $actually_installed = "";

    # install the packages by priority
    foreach $p (sort {$a <=> $b} keys %install_packages) {
        my $pname = $install_packages{$p};
        my $pversion = $packages{$pname}[$vrs];
        print my_colored("Package marked for installation: $pname\t[$pversion]\n",$WARNCOLOR);
        unless ($opt_dryrun)
        {
          $packages{$pname}[$fnc]->();
          $actually_installed .= "$pname ";
        }
    }

    unless ($opt_dryrun)
    {
      print my_colored("\n\nInstalled sucessfully: $actually_installed\n",$OKCOLOR);
      print my_colored("\n!!! FINISHED INSTALLING ALL SELECTED DEPENDENCIES !!!\n\n",$OKCOLOR);
    }
}

#==========================================================================
# Main execution
#==========================================================================

parse_commandline();

prepare();

install_packages();

if ($opt_genconf)      { generate_conf(); }

post();
