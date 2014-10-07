#!/usr/bin/env perl

# modules
use Term::ANSIColor;
use Getopt::Long;
use Tie::File;
use File::Basename;
use Shell qw(ls cat);

# constants
my $ERRORCOLOR="bold red";
my $OKCOLOR="bold green";
my $HEADINGCOLOR = "bold";
my $DEBUGCOLOR = "yellow";
my $WARNCOLOR = "bold yellow";

# command line options
my $opt_help = 0;
my $opt_dryrun = 0;
my $opt_file = "lolo.hh";
my $opt_mod = "LoloModule";

# Parse command line
$opt_help=1 unless GetOptions (
        'help'                  => \$opt_help,
        'dryrun'                => \$opt_dryrun,
        'file=s'                => \$opt_file,
        'mod=s'                 => \$opt_mod,
      );

# show help if required
if ($opt_help != 0)
{
  print <<ZZZ;
convert-CFcase-CFL.pl : move the CFL from the subsystem to the convergence method
options:
        --help            Show this help.
        --file=           The cxx file to convert
                            Default: $opt_file
        --mod=            The module to use in the conversion
                            Default: $opt_module
ZZZ
  exit(0);
}


print $opt_file."\n";

# put all file into $file
open FILE, "<$opt_file" or die "can't open $opt_file $!";
my $file = do { local $/; <FILE> };
close FILE;

# Detect what is the convergence method
$file =~ m/ConvergenceMethod(\s*)=(\s*)(\S+)/;
$conv = $3;
print "Convergence Method: ";

# Substitute Subsys.CFL... by SubSys.ConvergenceMethod.Data.CFL
$file =~ s/\.CFL\.(\S+)/\.$conv\.Data\.CFL\.$1/g ;

print $file."\n";

unless ($opt_dryrun) {
    open FILE, ">$opt_file" or die "can't open $opt_file $!";
    print FILE $file;
    close FILE;
}

#Reopen the file as an array and move lines around
my @lines;
tie @lines, 'Tie::File', $opt_file or die ("Error opening file $opt_file - $!\n");

my $index = 0;
my $convLine = 100000;
my $cflLine = 100000;
my @newLines;

my $newLineIdx = 0;
foreach $line(@lines)
{
  if ($lines[$index] =~ m/ConvergenceMethod(\s*)=(\s*)(\S+)/)
  {
    $convLine = $index;
    @newLines[0] = $lines[$index];
  }
  $index = $index+1;
}

$index = 0;
foreach $line(@lines)
{
  if ($lines[$index] =~ m/\.Data\.CFL\.(\S+)/)
  {
    @cflLines[$newLineIdx] = $index;
    @newLines[$newLineIdx + 1] = $lines[$index];
    $newLineIdx++;
  }
  $index = $index+1;
}


my $count = @newLines;

#Remove the lines with the CFL
my $iLine = 0;
foreach $newLine(@newLines)
{
  if ($iLine<$count-1)
  {
#    print "Removing line: " , @cflLines[$iLine]-$iLine,"\n";
    splice(@lines,@cflLines[$iLine]-$iLine,1);
  }
  $iLine++;
}

my $nbLines= @lines;
for ($i=$nbLines; $i>0; $i--)
{
  if($i > $convLine-$count+1){
    $lines[$i + $count-1] = $lines[$i];
  }
}

#Add the removed lines at the correct position
splice(@lines,$convLine - $iLine+1,$count,@newLines);


# {
#
# #Add the removed lines at the correct position
# splice(@lines,$convLine - $iLine + $i,2,@newLines[$i],@newLines[$i+1]);
#
# }


