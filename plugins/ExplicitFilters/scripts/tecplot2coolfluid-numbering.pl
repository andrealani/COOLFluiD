#!/usr/bin/env perl

#==========================================================================
# Modules
#==========================================================================
use Term::ANSIColor;
use Getopt::Long;
use File::Path;
use File::Copy;
use Switch;
use File::Find;
use File::Basename;



#==========================================================================
# Global Variables
#==========================================================================
my $opt_help             = 0;
my $opt_fileSearch       = "bad_filters";
my $opt_zone             = 0;
my $opt_element          = 0;

my @files = ();
my $filename;
my $CF_element;

#==========================================================================
# Command Line
#==========================================================================
sub parse_commandline() # Parse command line
{
    $opt_help=1 unless GetOptions (
        'help'                  => \$opt_help,
        'zone=s'                => \$opt_zone,
        'element=s'             => \$opt_element,
        'filename=s'            => \$opt_fileSearch,
    );

    # show help if required
    if ($opt_help != 0)
    {
      print <<ZZZ;
plot-transferfunction.pl : Plot transfer function from explicit filtering
options:
        --help            Show this help.
        --zone=           The zone in tecplot
        --element=        The element in tecplot
        --filename=       Search files with this filename
ZZZ
    exit(0);
    }
}




#==========================================================================
# Helper funcions
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

sub is_mac()
{
    $arch = run_command("uname -s");
    if ($arch =~ m/Darwin/i) {
        return 1;
    } else {
        return 0;
    }
}

#==========================================================================

sub check_in_path($) {
    my ($prog) = @_;
    my $status = run_command("which $prog");
    if ($status eq "") {
        print "$prog not found\n";
        return 0;
    } else {
        return 1;
    }
}

#==========================================================================

sub find_file()
{
    # http://perldoc.perl.org/File/Find.html
    find(\&mySub,"."); #custom subroutine find, parse $dir
    $size = scalar @files;
    if (scalar $size > 1) {
        print "ERROR: Too many datafiles found\n";
        # print @files;
        exit;
    }
    elsif ($size == 0) {
        print "ERROR: No file datafile found\n";
        exit;
    }
    else {
        # print "@files\n";
        $filename = @files[0];
    }
}
sub mySub()
{
    print "$_\n" if(/$opt_fileSearch.*\.log$/i);
    push @files, $_ if(/$opt_fileSearch.*\.log$/i);
}

#==========================================================================

sub find_coolfluid_element() 
{
    print "Finding coolfluid element in $filename\n";
    open(IN,"<$filename") or die ("Error opening $filename\n");
    while (<IN>)
    {
        chomp;
        $line = $_;
        if ($line =~ m/^\s*zone\s*$opt_zone\s*tec_element\s*$opt_element\s*CF_element\s*([0-9]+)/) { $CF_element = $1;}
    }
    close(IN);

    print "CF_element = $CF_element\n";
  
}

#==========================================================================
# Main execution
#==========================================================================

print "Script to convert tecplot element numbering to coolfluid numbering\n";

parse_commandline();
find_file();
find_coolfluid_element();
