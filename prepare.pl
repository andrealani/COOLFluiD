#!/usr/bin/env perl
#
# TODO:
#  - check modules with name colision

# modules
use Term::ANSIColor;
use Getopt::Long;
use File::Path;
use File::Basename;
use Cwd;
use Data::Dumper;
use strict;
# use warnings;

# constant names of files
my $preparepl  = "prepare.pl";
my $cmakelist  = "CMakeLists.txt";

# svn server structure
my $svn_kernel_dir  = "src";
my $svn_plugins_dir = "plugins";

# dependency libraries
my @libraries = ( "mpi",
                  "curl",
                  "trilinos",
                  "parmetis",
                  "metis",
                  "mitrem",
                  "cgns",
                  "get",
                  "gsl",
                  "blitz",
                  "flens",
                  "tvmet",
                  "cblas",
                  "google_perftools",
                  "mutationpp",
                  "cuda");

# constants
my $ERRORCOLOR="bold red";
my $OKCOLOR="bold green";
my $HEADINGCOLOR = "bold";
my $SUBSECTIONCOLOR = "yellow";
my $SECTIONCOLOR = "bold yellow";
my $DEBUGCOLOR = "yellow";

# options
my %default_options = (
    'help'         => 0,
    'nocolor'      => 0,
    'mods-list'    => 0,
    'mods-update'  => 0,
    'mods-getall'  => 0,
    'mods-revision'=> 0,
    'coolfluid-version' => "",
    'dry-run'      => 0,
    'debug'        => 0,
    'verbose'      => 0,
    'allstatic'    => "",
    'allactive'    => 1,
    'find'         => 0,
    'config'       => 0,
    'nodeps'       => 0,
    'with_unit_tests' => "",
    'search_dirs'  => "src plugins",
    'install_api'  => "",
    'extra_search_dirs'  => "",
    'extra_mods_url'     => "",
    'extra_mods_list'    => "",
    'svnserver'    => "https://github.com/andrealani/COOLFluiD/trunk",
    'config_file'  => "coolfluid.conf",
    'build'        => "RelWithDebInfo",
    'buildtype'    => "",
    'basebuild_dir'=> getcwd(),
    'plugins_dir'  => "plugins",
    'coolfluid_dir'=> getcwd(),
    'install_dir'  => "",
    'cmake_generator' => "make",
    'cc'              => "",
    'cxx'             => "",
    'fc'              => "",
    'cudac'            => "",
    'withbuildtype'   => 1,
    'withmpi'         => "",
    'withcuda'        => "",
    'withcurl'        => 1,
    'with_mutationpp' => 0,
    'withdocs'        => "",
    'mpi_extra_libs'  => "",
    'nofortran'       => "",
    'explicit_templates' => "",
    'cflags'          => "",
    'cxxflags'        => "",
    'fflags'          => "",
    'cudacflags'      => "",
    'ldflags'         => "",
    'with_testcases'  => "",
    'internal_deps'   => "",
    'warnings'        => "",
    'assertions'      => "",
    'trace'           => "",
    'logall'          => "",
    'logdebug'        => "",
    'debug_macros'    => "",
    'profiling'       => "",
    'profiler_tool'   => "gprof",
    'blas_dir'        => "",
    'lapack_dir'      => "",
    'lapack_libraries' => "",
    'curl_dir'        => "",
    'boost_dir'      => "",
    'boost_includedir' => "",
    'boost_librarydir' => "",
    'get_dir'         => "",
    'cgns_dir'        => "",
    'trilinos_dir'    => "",
    'metis_dir'       => "",
    'parmetis_dir'    => "",
    'google_perftools_dir' => "",
    'plas_dir'        => "",
    'lesmodels_dir'   => "",
    'mitrem_dir'      => "",
    'pardiso_include_dir'  => "", # /opt/intel/mkl/10.x.x.xxx/include
    'pardiso_dir'          => "", # /opt/intel/mkl/10.x.x.xxx/lib/32
    'pardiso_gfortran_dir' => "", # /usr/lib/gcc/x/x
    'petsc_dir'            => "",
    'samg_dir'             => "",
    'samg_options'         => "",
    'superlu_dir'          => "",
    'mpi_dir'              => "",
    'flens_dir'            => "",
    'blitz_dir'            => "",
    'tvmet_dir'            => "",
    'cblas_dir'            => "",
    'cuda_dir'             => "",
    'gsl_dir'              => "",
    'gsl_librarydir'       => "/usr/lib64",
    'gsl_includedir'       => "/usr/include",
    'mutationpp_dir'       => "",
    'mutationpp_librarydir' => "/usr/lib64",
    'mutationpp_includedir' => "/usr/include",
    'command'              => "",
    'single_precision'     => 0
);

# add skip default skip options
foreach my $libname ( @libraries )
{
  $default_options{"$libname\_skip"} = 0;
}

# variables
my %options = %default_options;   # for command line options
my %user_pref;                    # for config file options

my $extrasearchdirs  = "";        # record the extra search dirs
my $extramodsurl  = "";          # record the extra mods dirs
my $install_api = "";             # which libraries should have their API installed

my %moddir;                       # in which dir the module resides
my %modlibs;                      # which libs are in each module
my %modapps;                      # which apps are in each module
my %modtests;                     # which tests are in each module

my @all_mods;                     # all the modules
my @all_libs;                     # all the plugin libs
my @all_apps;                     # all the plugin apps
my @all_tests;                    # all the plugin unit tests

# list of options to choose the compiler
my $comp_options = "";
# list of other options that affect configuration
my $other_options = "";
# list of options with which plugins to compile or skip
my $plugin_options = "";

#==========================================================================
# Helper funcions
#==========================================================================

sub parse_command_line_options()
{
   # parse command line
   $options{'help'}=1 unless GetOptions (
       'help'           => \$options{'help'},
       'debug'          => \$options{'debug'},
       'nocolor'        => \$options{'nocolor'},
       'mods-list'      => \$options{'mods-list'},
       'mods-update'    => \$options{'mods-update'},
       'mods-select=s'  => \$options{'mods-select'},
       'mods-getall'    => \$options{'mods-getall'},
       'mods-revision=i'=> \$options{'mods-revision'},
       'coolfluid-version=s' => \$options{'coolfluid-version'},
       'config'         => \$options{'config'},
       'nodeps'         => \$options{'nodeps'},
       'find'           => \$options{'find'},
       'dry-run'        => \$options{'dry-run'},
       'with_unit_tests'    => \$options{'with_unit_tests'},
       'verbose'        => \$options{'verbose'},
       'allstatic'      => \$options{'allstatic'},
       'allactive'      => \$options{'allactive'},
       'search-dirs=s'  => \$options{'search_dirs'},
       'extra-search-dirs=s'  => \$options{'extra_search_dirs'},
       'extra-mods-url=s'  => \$options{'extra_mods_url'},
       'extra-mods-list=s'  => \$options{'extra_mods_list'},
       'install-api=s'  => \$options{'install_api'},
       'config-file=s'  => \$options{'config_file'},
       'build=s'        => \$options{'build'},
       'buildtype=s'    => \$options{'buildtype'},
       'coolfluid_dir=s'=> \$options{'coolfluid_dir'},
       'install_dir=s'  => \$options{'install_dir'},
       'cmake_generator=s' => \$options{'cmake_generator'},
       );
   
   # remove duplicated entries in options
   while ( my ($key, $value) = each(%default_options) ) {
       if ( $value eq $options{$key}) {
	   delete $options{$key};
       };
   }
   
   # show help if required
   if ($options{'help'} != 0)
   {
       print <<ZZZ;
       prepare.pl : this script prepares the COOLFluiD Build Environment
	   
	   OPTIONS
	   
	 behavior:
	   --help              Show this help.
	   --nocolor           Don't color output
           --dry-run           Don't actually execute the system calls.
	 Just output what you would do.
            --nodeps            Don't do check for action dependencies.
            --verbose           Print extra information.

  configuration:
         --allactive         Makes all libraries active by default.
                               This is the default.
         --allstatic         Makes all libraries static by default, instead of dynamic.
                               They still can be desactivated in the configuration file.
         --with_unit_tests   Create the unit tests. [Default: $default_options{'with_unit_tests'}]
         --config-file=      User config file to overide default configuration options
                                Default: $default_options{'config_file'}
         --basebuild_dir=    Location where to build the sources
                               Default: Current directory
         --coolfluid_dir=    Location of the COOLFluiD sources directory
                               Default: Current directory
         --install-dir=      Location of the COOLFluiD installation directory
                               Default: $default_options{'install_dir'}
         --cmake_generator=  Chooses which type of build system will cmake generate
                               'make'  for 'Unix Makefiles'
                               'kdev'  for 'KDevelop3'
                               'xcode' for 'XCode'
                               'cb'    for 'Code::Blocks'
                               'vs8'   for 'Visual Studio 8 2005'
                               'vs9'   for 'Visual Studio 9 2008'
                               'eclipse-make' for 'Eclipse CDT4 - Unix Makefiles'
                               Default: $default_options{'cmake_generator'}
         --mods-getall       Gets all the modules found in the plugins repository
         --mods-revision     Chechout or update the modules to a specific svn revision
         --mods-select       Comma separated list of modules to apply the svn actions

         --coolfluid-version Select a tagged version of coolfluid to check out. Ex: "2009.0"

         --extra-search-dirs=   Comma separated list of directories where to search for modules
                                  Default: $default_options{'extra_search_dirs'}

         --extra-mods-url=   Comma separated list of svn directories where to search for modules
                                  Default: $default_options{'extra_mods_url'}
  
         --extra-mods-list=   Comma separated list of extra modules
                                  Default: $default_options{'extra_mods_list'}

        


  ACTIONS

  build actions:
         --find            Just report which modules and plugins inside each module you find.
         --config          Performs the configuration.
         --build=          Build configuration for the compilation. Implies --config
                              Default: $default_options{'build'}

  modules actions:   ( needs network access )
         --mods-list       Lists the available modules in the svn server
         --mods-update     Checkout or Updates the selected modules from the svn server

ZZZ
   exit(0);
   }
}

#==========================================================================

sub wordcaps {
    my $line = shift;
    $line =~ tr/a-z/A-Z/; # all letters
  # $line =~ s/\b(\w)/\U$1/g; # only first letter
    return $line;
}

#==========================================================================

sub my_colored($$)
{
  return (get_option('nocolor') ? shift : colored($_[0], $_[1]));
}

#==========================================================================

sub get_build() # get the build option
{
   my $build;
    if ( exists($options{'build'}) ) {
      $build = $options{'build'};
    } else {
      if ( exists($user_pref{'build'}) ) {
        $build = $user_pref{'build'};
      } else {
        $build = $default_options{'build'};
      }
    }
    $build =~ s/\/$//;
    return $build;
}

#==========================================================================

sub get_build_path() # gets the build path
{
   my $build = get_build();
   my $basebuild_dir = get_option('basebuild_dir');
   return "$basebuild_dir/$build";
}

#==========================================================================

sub get_search_paths() # gets the paths where to search for sources
{
   my $coolfluid_dir = get_option('coolfluid_dir');
   my $search_dirs = get_option('search_dirs');
   my @sdirs = split (" " , $search_dirs );
   my $search_paths;

   foreach my $dir ( @sdirs )
   {
     $search_paths .= "$coolfluid_dir/$dir ";
   }
   return "$search_paths";
}

#==========================================================================

sub get_option($) # get the option, if the case, overidden by the user
{
    my ($key) = @_;
    my $build = get_build();

    # per build options in the config file have top priority
    if ( exists($user_pref{"$build\_$key"}) ) {
      return $user_pref{"$build\_$key"};
    }

    # next come the command line options
    if ( exists($options{"$key"}) ) {
      return $options{"$key"};
    }

    # followed by the normal config file options
    if ( exists($user_pref{"$key"}) ) {
      return $user_pref{"$key"};
    }

    # followed by the default options
    if ( exists($default_options{"$key"}) ) {
      return $default_options{$key};
    }

    # die if option doesnt exist
    die("Unknown option $key required");
}

#==========================================================================

sub get_arch() # returns the current architecture
{
    my $args="uname -m";

    my $output;
    # print my_colored("Executing : $args",$OKCOLOR);
    my $command = join("",$args,"|");
    my $pid=open READER, $command or die "Can't run the program: $args $!\n";
    while(<READER>){ $output.=$_; }
    close READER;
    # print my_colored($output,$OK_COLOR);
    my $arch = $output;

    chomp($arch);
    return $arch;
}

#==========================================================================

sub parse_config_file($) # parse the config file to get the user overiding options
{
    my ($filename) = @_;
    if ( -e "$filename" ) {

    print "Reading config file : $filename\n";

    open CONFIG, "<", $filename or die ("Error opening config file $filename!\n");

    my $arch = get_arch();

    while (<CONFIG>) {
        chomp;                  # no newline
        s/#.*//;                # no comments
        s/^\s+//;               # no leading white
        s/\s+$//;               # no trailing white
        next unless length;     # anything left?
        s/\$ARCH/$arch/;        # subst $ARCH
        s/\$(\w+)/$ENV{$1}/g;   # subst all $X for environmental variables
        my ($key, $value) = split(/\s*=\s*/, $_, 2);
        $user_pref{$key} = $value;
    }
    close CONFIG;
   }
}

#==========================================================================

sub find_modules()
{
   print my_colored("\nSearching for module directories containing $cmakelist\n",$SECTIONCOLOR);

   # check that we are in the toplevel directory
   die (my_colored("Cannot find $preparepl; You need to run this command in the toplevel COOLFluiD directory!\n",$ERRORCOLOR))  unless -f "$preparepl";

   if (get_option('extra_search_dirs'))
   {
     # fix the search dir option which comes comma separated
     my @array_extradirs = split(',',get_option('extra_search_dirs'));
     $extrasearchdirs = join(" ",@array_extradirs);
   }
   
   # check out additional user-defined modules into the "plugins" folder
   if (get_option('extra_mods_list') and get_option('extra_mods_url'))
   {
     # fix the mods dir option which comes comma separated
     my @array_extramlist = split(',',get_option('extra_mods_list'));
     my @array_extramurl  = split(',',get_option('extra_mods_url'));
     my $coolfluid_dir = get_option('coolfluid_dir');
     for (my $idir = 0; $idir <= $#array_extramlist; $idir++)
     {
       my $mdir = "$coolfluid_dir/plugins/$array_extramlist[$idir]";
       print my_colored("\nChecking whether dir $mdir exists\n",$SECTIONCOLOR);
       if (!(-d $mdir) ) 
       {	 
        print my_colored("\nChecking out extra module $array_extramlist[$idir] from $array_extramurl[$idir]\n",$SECTIONCOLOR);
        run_command("svn co $array_extramurl[$idir] $mdir");  
       }
     }
   }
   
   my $searchdirs = get_search_paths() . " " . $extrasearchdirs;
   my $search_unittests = get_option('with_unit_tests');
   
   open FILE, "-|", "find -L $searchdirs -name $cmakelist" or die ("Error searching for  $cmakelist - $!");

   while (my $file=<FILE>)
   {
      my $dir     = dirname($file);
      my $modname = $dir;
      $modname =~ s/\//_/g;

      #print my_colored("\n##### FILE   => $file\n",$SECTIONCOLOR);
      #print my_colored("\n##### DIR    => $dir\n",$SECTIONCOLOR);
      #print my_colored("\n##### MODULE => $modname\n",$SECTIONCOLOR);

      push(@all_mods, $modname);

      $moddir{$modname} = $dir;

      # array handle for local libs, apps and tests
      $modlibs{"$modname"} = [];
      $modapps{"$modname"} = [];
      $modtests{"$modname"} = [];
      my $local_libs = $modlibs{"$modname"};
      my $local_apps = $modapps{"$modname"};
      my $local_tests = $modtests{"$modname"};

      # parse module file to find libs, apps or unit tests
      open MODMK, "<$file" or die "can't open $file $!";
      while (<MODMK>) {
        my $line = $_;

        # found a plugin library
        if ($line =~ m/(\s*)CF_ADD_PLUGIN_LIBRARY(\s*)\((\s*)(.*)(\s*)\)/ )
        {
          my $libs = $4;
          $libs =~ s/(\s+)/\:/g;
          push(@all_libs,split('\:',$libs));
          push(@$local_libs,split('\:',$libs));
        }

        # found a kernel library
        if ($line =~ m/(\s*)CF_ADD_KERNEL_LIBRARY(\s*)\((\s*)(.*)(\s*)\)/ )
        {
          my $libs = $4;
          $libs =~ s/(\s+)/\:/g;
          push(@all_libs,split('\:',$libs));
          push(@$local_libs,split('\:',$libs));
        }

        if ($line =~ m/(\s*)CF_ADD_PLUGIN_APP(\s*)\((\s*)(.*)(\s*)\)/ )
        {
          my $apps = $4;
          $apps =~ s/(\s+)/\:/g;
          push(@all_apps,split('\:',$apps));
          push(@$local_apps,split('\:',$apps));
        }

        if ( $search_unittests and $line =~ m/(\s*)CF_ADD_TEST(\s*)\((\s*)(.*)(\s*)\)/ )
        {
          my $tests = $4;
          $tests =~ s/(\s+)/\:/g;
          push(@all_tests,split('\:',$tests));
          push(@$local_tests,split('\:',$tests));
        }
      }
      close MODMK;

   # print info about this config file
   if (get_option('verbose')) {
       my $nblibs = scalar (@$local_libs);
       my $nbapps = scalar (@$local_apps);
       my $nbtests = scalar (@$local_tests);
       print "Module file " . $file . " defines:\n\t$nblibs lib(s) : @$local_libs\n\t$nbapps app(s) : @$local_apps\n\t$nbtests test(s) : @$local_tests\n";
   }

   }
   close FILE;

   @all_libs  = sort @all_libs;
   @all_apps  = sort @all_apps;
   @all_tests = sort @all_tests;
   @all_mods  = sort @all_mods;

   my $nblibs  = scalar @all_libs;
   my $nbapps  = scalar @all_apps;
   my $nbtests = scalar @all_tests;
   my $nbmods  = scalar @all_mods;

   print my_colored("Searched in directories $searchdirs\n",$OKCOLOR);
   print my_colored("Found $nbmods MODULE directories with $cmakelist defining a total of $nblibs lib(s), $nbapps app(s) and $nbtests unit tests(s)\n",$OKCOLOR);

   # print summary of module info
   if (get_option('verbose')) {
      print my_colored("\nList of Modules\n",$SECTIONCOLOR);

      foreach my $modname (@all_mods)
      {
         print my_colored("Module $modname : ",$SUBSECTIONCOLOR) if get_option('verbose');

         my $libs_handle  = $modlibs{$modname};
         my $apps_handle  = $modapps{$modname};
         my $tests_handle = $modtests{$modname};

         if (scalar(@$libs_handle))  { print " " . scalar(@$libs_handle)  ." libs (". join(',',@$libs_handle) .") ";  }
         if (scalar(@$apps_handle))  { print " " . scalar(@$apps_handle)  ." apps (". join(',',@$apps_handle) .") ";  }
         if (scalar(@$tests_handle)) { print " " . scalar(@$tests_handle) ." tests (".join(',',@$tests_handle).") ";  }
         print "\n";
      }
   }

   print "\n";

   print my_colored("All libs  :\n  ",$OKCOLOR) . " @all_libs\n";
   print my_colored("All apps  :\n  ",$OKCOLOR) . " @all_apps\n";
   print my_colored("All tests :\n  ",$OKCOLOR) . " @all_tests\n";
}

#==========================================================================

sub set_default_plugins_on() # put the plugins into the default options
{
   # set the default value for the libs
   # should be the default value of allactive
   foreach (@all_libs)
   {
      if (get_option('allactive'))
      { $default_options{"lib_$_"} = 'on';  }
      else
      { $default_options{"lib_$_"} = 'off'; }
   }

   # by default apps are active
   foreach (@all_apps)
   {
      $default_options{"app_$_"} = 'on';
   }

   # by default tests are active
   foreach (@all_tests)
   {
      $default_options{"test_$_"} = 'on';
   }
}

#==========================================================================

sub print_hash($$)
{
   my ($ref_hash, $ordered) = @_;

# unordered
  if ($ordered)
  {
   my %hash = %{$ref_hash};
   my @keys = sort keys %{$ref_hash};
   foreach my $key (@keys)
   {
         my $value = $hash{$key};
         print "$key => $value\n";
   }
  }
  else
  {
      while ( my ($key, $value) = each(%{$ref_hash}) ) {
         print "$key => $value\n";
      }
  }
}

#==========================================================================

sub print_debug_options()
{
   print my_colored("DefaultOptions\n",$OKCOLOR);
   print_hash(\%default_options,1);

   print my_colored("Options\n",$OKCOLOR);
   print_hash(\%options,1);

   print my_colored("UserOptions\n",$OKCOLOR);
   print_hash(\%user_pref,1);
}

#==========================================================================

sub run_command($)
{
    my ($args)=@_;
    my $output;
    if (get_option('dry-run'))
    {
       print my_colored("[dry-run] would have executed : ",$OKCOLOR).$args."\n";
    }
    else {
      print my_colored("[verbose] executing : ",$OKCOLOR).$args."\n" if get_option('verbose');
      my $command = join("",$args,"|");
      my $pid=open READER, $command or die "Can't run the program: $args $!\n";
      while(<READER>){
         $output.=$_;
      }
      close READER;
    }
    return $output;
}

#==========================================================================

sub print_var($$) # prints a variable and its value in color
{
    my ($var,$value)=@_;
    print my_colored($var,$OKCOLOR); print " : ";
    print my_colored($value,$DEBUGCOLOR); print "\n";
}

#==========================================================================

sub suredir($) # make sure a directory exists, if not, create it
{
    my ($dir)=@_;
    if (-d $dir)
    {
        print my_colored($dir,$OKCOLOR); print "\n";
    }
    else
    {
        print my_colored($dir,$ERRORCOLOR); print "\n";
        mkdir($dir, 0755) || die "Cannot mkdir $dir: $!";
    }
}

#==========================================================================

sub cfg_lib($$) # put the disabled libraries in the list
{
  my ($lib, $list) = @_;
  if ( ! ( get_option("lib_$lib") eq 'on' or get_option("lib_$lib") eq 'ON' ) )
  { push(@$list,"$lib") }
}

#==========================================================================

sub cfg_app($$) # put the disabled application in the list
{
  my ($app, $list) = @_;
  if ( ! ( get_option("app_$app") eq 'on' or get_option("app_$app") eq 'ON' ) )
  { push(@$list,"$app") }
}

#==========================================================================

sub cfg_test($$) # put the disabled unit test in the list
{
  my ($test, $list) = @_;
  if ( ! ( get_option("test_$test") eq 'on' or get_option("test_$test") eq 'ON' ) )
  { push(@$list,"$test") }
}


#==========================================================================

sub add_cmake_option($)
{ 
  my ($cmopt) = @_;

  my $cmopt_caps = wordcaps ( $cmopt );
  my $cmopt_data = get_option( $cmopt );
  
  if ( !($cmopt_data eq '') )
  {  $plugin_options .= " -D$cmopt_caps=\"$cmopt_data\"";  }
}

#==========================================================================

sub setup_deps()
{ 
  if(get_option('boost_dir')) 
    { $plugin_options .= " -DBOOST_ROOT=\"".get_option('boost_dir')."\""; }

  # dependency variables
  my @dep_variables = ( "boost_includedir",
			"boost_librarydir",
			"petsc_dir",
			"plas_dir",
			"pardiso_include_dir",
			"pardiso_dir",
			"pardiso_gfortran_dir",
			"lesmodels_dir",
			"samg_dir",
			"samg_options",
			"blas_dir",
			"lapack_dir",
			"gsl_includedir",
			"gsl_librarydir",
                        "mutationpp_librarydir",
                        "mutationpp_includedir",
			"lapack_libraries" );
  
  foreach (@dep_variables) {  add_cmake_option($_); }
  
  foreach (@libraries) {  setup_library($_); }
}

#==========================================================================

sub setup_library($) # setup a dependecy library
{
    my ($lib)=@_;

    my $lkey = "$lib\_dir";
    my $lib_dir  = get_option($lkey);
    my $lkey_s = "$lib\_skip";
    my $lib_skip = get_option($lkey_s);

    my $libcaps = wordcaps($lib);

    if    ($lib_skip) {  $plugin_options .= " -DCF_SKIP_$libcaps=ON";   }
    else              
    {
	if ( !($lib_dir eq '') )
	{ $plugin_options .= " -D$libcaps\_HOME=\"$lib_dir\""; }  
    }
}

#==========================================================================

sub setup_option($$)
{
  my ($opt_name, $opt_cmake_name) = @_;
  my $lopt = get_option("$opt_name");

  if ( ($lopt eq 'on') or ($lopt eq 'off') or ($lopt eq '1') or ($lopt eq '0') )
  {
      $other_options .= " -D$opt_cmake_name=$lopt";
  }
  else
  {
    unless ( $lopt eq "" ) { die "Option \'$opt_name\' must have value 'on' of 'off'"; }
  }
}

#==========================================================================

sub setup_cfgoptions()
{
  setup_deps();

  my $enablelibs   = [];
  my $disablelibs  = [];
  my $disableapps  = [];
  my $disabletests = [];

  foreach (@all_libs)  { cfg_lib($_,$disablelibs);   }
  foreach (@all_apps)  { cfg_app($_,$disableapps);   }
  foreach (@all_tests) { cfg_test($_,$disabletests); }

  # enable libs are all minus the disabled
  # compute the intersection of the arrays
  my @isect = ();
  my %count = ();
  foreach my $e (@all_libs, @$disablelibs) { $count{$e}++ }
  foreach my $e (keys %count)
  {
    push @{ $count{$e} == 2 ? \@isect : \@$enablelibs }, $e;
  }

  # plugins are ON by default, no need to add them to the configuration

  foreach (@$disablelibs)  
  {  
   $plugin_options .= " -DCF_BUILD_$_=OFF"; 
   $plugin_options .= " -DCF_COMPILES_$_=OFF";
  }

  foreach (@$disableapps)  {  $plugin_options .= " -DCF_BUILD_$_=OFF"; }

  foreach (@$disabletests) {  $plugin_options .= " -DCF_BUILD_$_=OFF"; }

  print my_colored("List of enabled libs:\n", $OKCOLOR);
  foreach (@$enablelibs) 
  {  
   print "$_ "; 
   $plugin_options .= " -DCF_COMPILES_$_=ON";
  }
  print "\n";

  if ( get_option('install_api') )
  {
    # fix the option which comes comma separated
    my @array_install_api = split(',',get_option('install_api'));
    foreach my $lib_api (@array_install_api)
    {
	  $plugin_options .= " -DCF_BUILD_$lib_api\_API=ON";
    }
  }

  if (get_option('extra_search_dirs'))
  {
    # fix the search dir option which comes comma separated
    my @array_extradirs = split(',',get_option('extra_search_dirs'));
    $extrasearchdirs = join(";",@array_extradirs);

    $other_options .= " -DCF_EXTRA_SEARCH_DIRS=$extrasearchdirs";
  }

  if (get_option('mpi_extra_libs'))
  {
    # fix the MPI extra libs option which comes with comma separated list
    my @array_mpiextralibs = split(',',get_option('mpi_extra_libs'));
    my $mpiextralibs = join(";",@array_mpiextralibs);

    $other_options .= " -DMPI_EXTRA_LIBRARY_NAMES=$mpiextralibs";
  }

  if (get_option('profiling'))
  {
    my $opt_profiling = get_option('profiling');
    my $opt_profiler_tool = get_option('profiler_tool');
    unless ( ($opt_profiling eq "on") or ($opt_profiling eq "off") ) { die "Option 'profiling' must have value 'on' of 'off'"; }
    $other_options .= " -DCF_ENABLE_PROFILING=$opt_profiling -DCF_PROFILER_TOOL=$opt_profiler_tool";
  }

  setup_option('nofortran',           'CF_SKIP_FORTRAN');
  setup_option('withmpi',             'CF_ENABLE_MPI');
  setup_option('withcuda',            'CF_ENABLE_CUDA');
  setup_option('withcurl',            'CF_ENABLE_CURL');
  setup_option('with_mutationpp',     'CF_ENABLE_MUTATIONPP');
  setup_option('withdocs',            'CF_ENABLE_DOCS');
  setup_option('explicit_templates',  'CF_ENABLE_EXPLICIT_TEMPLATES');
  setup_option('with_testcases',      'CF_ENABLE_TESTCASES');
  setup_option('with_unit_tests',     'CF_ENABLE_UNITTESTS');
  setup_option('internal_deps',       'CF_ENABLE_INTERNAL_DEPS');
  setup_option('assertions',          'CF_ENABLE_ASSERTIONS');
  setup_option('trace',               'CF_ENABLE_TRACE');
  setup_option('logall',              'CF_ENABLE_LOGALL');
  setup_option('logdebug',            'CF_ENABLE_LOGDEBUG');
  setup_option('debug_macros',        'CF_ENABLE_DEBUG_MACROS');
  setup_option('allstatic',           'CF_ENABLE_STATIC');
  setup_option('warnings',            'CF_ENABLE_WARNINGS');
  setup_option('single_precision',    'CF_PRECISION_SINGLE');
}

#==========================================================================

sub run_configuration() # run configuration
{
  my $build = get_build();
  my $coolfluid_dir = get_option('coolfluid_dir');
  my $buildtype = get_option('buildtype');
  my $caps_build;
  if ( $buildtype eq "" ) { $caps_build = uc ($build); }
  else { $caps_build = uc ($buildtype) }

  print my_colored("\n--- CMAKE Configuration ---\n",$HEADINGCOLOR);
  print my_colored("Build profile : $build\n",$OKCOLOR);
  print my_colored("Build type    : $caps_build\n",$OKCOLOR);

  my $install_dir = get_option('install_dir');
  my $path = get_build_path();

  mkpath $path; # make sure dir exists
  chdir($path) || die "Cannot chdir to $path ($!)";

  $comp_options = "";

  my $cc   = get_option('cc');
  my $cxx  = get_option('cxx');
  my $fc   = get_option('fc');
  my $cudac = get_option('cudac'); 

  my $cflags   = get_option('cflags');
  my $cxxflags = get_option('cxxflags');
  my $fflags   = get_option('fflags');
  my $cudacflags = get_option('cudacflags');

  if (get_option('cflags'))    { $comp_options .= "-DCF_C_FLAGS:STRING=\"$cflags\" ";   }
  if (get_option('cxxflags'))  { $comp_options .= "-DCF_CXX_FLAGS:STRING=\"$cxxflags\" "; }
  if (get_option('fflags'))    { $comp_options .= "-DCF_Fortran_FLAGS:STRING=\"$fflags\" ";   }
  if (get_option('cudacflags'))    { $comp_options .= "-DCF_CUDAC_FLAGS:STRING=\"$cudacflags\" ";   }
 
  my $withcuda = get_option('withcuda');
  if ($withcuda eq '1')    { $comp_options .= "-DCMAKE_CUDA_COMPILER_ENV_VAR=\"CUDACC\" ";   }

  if (get_option('install_dir')) { $other_options .= " -DCMAKE_INSTALL_PREFIX=$install_dir"; }

  my $gen = '';
  my $genopt = get_option('cmake_generator');
  if ($genopt eq 'make')  { $gen = 'Unix Makefiles'; }
  if ($genopt eq 'kdev')  { $gen = 'KDevelop3'; }
  if ($genopt eq 'xcode') { $gen = 'Xcode'; }
  if ($genopt eq 'cb')    { $gen = 'CodeBlocks - Unix Makefiles'; }
  if ($genopt eq 'vs8')   { $gen = 'Visual Studio 8 2005'; }
  if ($genopt eq 'vs9')   { $gen = 'Visual Studio 9 2008'; }
  if ($genopt eq 'vs8win64')     { $gen = 'Visual Studio 8 2005 Win64'; }
  if ($genopt eq 'vs8win64')     { $gen = 'Visual Studio 9 2008 Win64'; }
  if ($genopt eq 'eclipse-make') { $gen = 'Eclipse CDT4 - Unix Makefiles'; }
  if ( $gen eq '' ) { die "CMake generator option has value \'$genopt\' which does not map to a known generator" }

  if ( get_option('withbuildtype') ) 
  {
	unless ( ( $caps_build eq "RELWITHDEBINFO" )  )
	{ $other_options .= " -DCMAKE_BUILD_TYPE=$caps_build";  }
  }

  my $args = "cmake $comp_options -DCMAKE_C_COMPILER=$cc -DCMAKE_CXX_COMPILER=$cxx -DCMAKE_Fortran_COMPILER=$fc -DCMAKE_CUDA_COMPILER=$cudac -G\"$gen\" $other_options $plugin_options $coolfluid_dir";

  print run_command($args);
}

#==========================================================================

sub print_configinfo() # print information about the
{
    my $build = get_build();
    my $coolfluid_dir = get_option('coolfluid_dir');
    my $install_dir = get_option('install_dir');
    my $extra_search_dirs = get_option('extra_search_dirs'); 
    my $extra_mods_url = get_option('extra_mods_url');
    my $buildpath = get_build_path();

    print my_colored("\nDirectories\n",$SECTIONCOLOR);

    print_var("COOLFluiD    dir  ", $coolfluid_dir);
    print_var("Build        dir  ", $buildpath);
    print_var("Extra        dirs ", "$extra_search_dirs");
    print_var("Extra module dirs ", "$extra_mods_url");
    print_var("Installation dir ", $install_dir);

  my $cc   = get_option('cc');
  my $cxx  = get_option('cxx');
  my $fc   = get_option('fc');
  my $cudac   = get_option('cudac');
  print my_colored("\nCompilers\n",$SECTIONCOLOR);
  print_var('CC ',$cc);
  print_var('CXX',$cxx);
  print_var('FC ',$fc);
  print_var('CUDAC',$cudac);


  # Options
  if ( get_option('verbose') )
  {
     # User prefs
     while ( my ($key, $value) = each(%options) ) 
     {
         print_var($key,get_option($key));
     }

     # User prefs
     while ( my ($key, $value) = each(%user_pref) ) 
     {
         print_var($key,$value);
     }
  }
}

#==========================================================================

sub find_phase() # find phase
{
   print my_colored("\n--- Module Discovery Phase ---\n",$HEADINGCOLOR);

   # change into the coolfluid dir
   my $dir = get_option('coolfluid_dir');
   chdir ($dir) || die "Cannot change into COOLFluiD dir - $dir : $!" ;

   find_modules();
   set_default_plugins_on();
}

#==========================================================================

sub config_phase() # configuration phase
{
   my $build = get_build();
   print my_colored("\n--- Configuration Phase ---\n",$HEADINGCOLOR);
   print my_colored("Build profile : $build\n",$OKCOLOR);

   setup_cfgoptions();

   print_configinfo();

   run_configuration();
}

#==========================================================================

sub mods_list() # lists the available modules in the subversion server
{
   print my_colored("\nListing available modules in subversion server\n",$SECTIONCOLOR);

   my $svnserver = get_option('svnserver');
   my $args = "svn list $svnserver/$svn_plugins_dir";
   print run_command($args);
}

#==========================================================================

sub mods_svnupdate ()
{
   print my_colored("\nChecking out modules from subversion server\n",$SECTIONCOLOR);

   my $svnserver = get_option('svnserver');
   my $coolfluid_dir = get_option('coolfluid_dir');
   my $plugins_dir = "$coolfluid_dir/".get_option('plugins_dir');
   my $cf_version = get_option("coolfluid-version");
   my $rev = get_option('mods-revision');

   # find last revision
   my $revout = run_command("svn info $svnserver");
   $revout    =~ m/Revision\:\s*(\d+)/;
   my $last_rev = $1;
   print "Last revision is $last_rev\n";

   # if there is a not a revision, find out which in the current one
   if ($rev eq 0) { $rev = $last_rev  }
   else { print "User selected revision $rev\n" }

   # if there is a version, find out which in which revision was modified
   #if (!($cf_version eq ''))
   #{
   #    my $versout = run_command("svn info $svnserver/../tags/$cf_version");
   #    $versout    =~ m/Last Changed Rev\:\s*(\d+)/;
   #    $rev = $1;
   #    print "User selected version $cf_version which has revision $rev\n";
   #}
   print "\n";

  # run on the main directory
  my $args;
  # add specific revision if specified
  #if (!($cf_version eq ''))
  #	{ $args = "svn switch -r$rev $svnserver/../tags/$cf_version"  }
  #else
  # 	{
          $args = " svn switch -r$rev $svnserver";
  #      }
  my $output = run_command($args);
  unless ($output eq '') { print "$output\n"; }

   # gets the modules list
   my $args = "svn list $svnserver/$svn_plugins_dir";
   my $modules  = run_command($args);
   # remove trailing slash
   $modules =~ s/\/\n/\n/g;
   my @mods = split("\n",$modules);

   # by default the modules are not retrieved unless getall option is active
   my $defvalue;
   if (get_option('mods-getall')) { $defvalue = 'on'; } else { $defvalue = 'off'; }
   foreach (@mods)
   {
     $default_options{"mod_$_"} = $defvalue;
   }

   # check if a module has been removed from the server
   for my $each_file ( glob($plugins_dir.'/*') )
   {
    ## if the $each_file is a directory
    if( -d $each_file)
    {
      my $moddir = $each_file;
      $moddir =~ s/$plugins_dir\///;
      my $found = 0;
      foreach my $mod (@mods)
      {
        if ( $mod eq $moddir ) { $found = 1; }
      }
      unless ($found)
      {
          my $ans;
          while ( ! ( $ans eq 'y' or $ans eq 'Y' or $ans eq 'n' or $ans eq 'N' )  )
          {
          print "Module \'$moddir\' does not exist anymore in the svn server\n";
          print "Do you want to remove the directory? [y/n]\n";
          $ans=<STDIN>;
          chomp($ans);
#          print "Answer was [$ans]\n";
          if ( $ans eq 'y' or $ans eq 'Y' )
          {
            print "removing directory $plugins_dir/$moddir\n";
            rmtree( "$plugins_dir/$moddir" );
          }
          elsif ( $ans eq 'n' or $ans eq 'N' )
          {
            print "keeping directory $plugins_dir/$moddir\n";
          }
          else { print "invalid answer [$ans]\n"; }
          print "\n";
        }
      }
    }
  }

   # try to checkout each module
   foreach my $mod (@mods)
   {
     print "Module $mod ... ";
     if ( get_option("mod_$mod") eq 'on' or get_option("mod_$mod") eq 'ON' )
     {
       my $output;
       my $ddir;

       # if exists then just update it
       if ( -e "$plugins_dir/$mod" )
       {
           #  die "$plugins_dir/$mod exists but is not a directory " unless ( -d "$plugins_dir/$mod" );
           # if is a directory update it
           print "updating ";
           $ddir = "$plugins_dir/$mod";
	         if ($cf_version eq '')
		         { $args = "svn update -r$rev" }
#           else
#		         { $args = "svn switch -r$rev $svnserver/../tags/$cf_version/$svn_plugins_dir/$mod " }
       }
       else # check it out
       {
          print "checking out ";
          $ddir = "$coolfluid_dir";
          if ( $cf_version  eq '' )
          { $args = "svn co -r$rev $svnserver/$svn_plugins_dir/$mod $plugins_dir/$mod"; }
#          else
#          { $args = "svn co -r$rev $svnserver/../tags/$cf_version/$svn_plugins_dir/$mod $plugins_dir/$mod"; }
       }

       # run the command
       chdir $ddir;
       print "[$args]\n";
       $output = run_command($args);
       unless ($output eq '') { print "$output\n"; }
       chdir $coolfluid_dir;
     }
     else
     {
       print "skipping\n";
     }
   }
}

#==========================================================================

# fix dependencies in the actions 'build' implies 'config' implies 'find'
sub fix_action_dependencies()
{
  # if user passes only --build, interpert it has to do all actions
  if ( !(get_build() eq $default_options{'build'})
    && !(get_option('config'))
    && !(get_option('find'))
    ) 
  { 
	$options{'config'} = 1 ;  
  };

  if (get_option('config'))        { $options{'find'} = 1};
}

#==========================================================================
# Main execution
#==========================================================================

parse_command_line_options();

print my_colored("\n === COOLFluiD Build Environment ===\n\n",$HEADINGCOLOR);

parse_config_file(get_option('config_file'));

fix_action_dependencies() unless (get_option('nodeps'));

if (get_option('debug')) { print_debug_options(); }

if (get_option('mods-list'))     { mods_list()     };
if (get_option('mods-update'))   { mods_svnupdate() };
if (get_option('find') )         { find_phase();   };
if (get_option('config'))        { config_phase();  }

exit 0;
