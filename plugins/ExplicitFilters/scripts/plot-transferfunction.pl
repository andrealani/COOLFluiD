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
my $opt_fileSearch       = "transferFunction";
my $opt_cell             = 0;
my $opt_plotter          = "dgraph";
my $opt_gnuplot_terminal = "aqua dashed";
my $opt_yrange        = "[:]";

my @files = ();
my $filename;
my $tmpfile = "tmp.dat";
my $FGR_target;
my $FGR_xy;
my $FGR_x;
my $dgraph_script = "$ENV{HOME}/workspace/coolfluid/plugins/ExplicitFilters/scripts/transferfunction-template.dgraph";
my $gnuplot_script = "script.gnuplot";
my $ymin = "-∞";
my $ymax = "∞";




#==========================================================================
# Command Line
#==========================================================================
sub parse_commandline() # Parse command line
{
    $opt_help=1 unless GetOptions (
        'help'                  => \$opt_help,
        'cell=s'                => \$opt_cell,
        'filename=s'            => \$opt_fileSearch,
        'plotter=s'             => \$opt_plotter,
        'term=s'                => \$opt_gnuplot_terminal,
        'yrange=s'                => \$opt_yrange,
    );

    # show help if required
    if ($opt_help != 0)
    {
      print <<ZZZ;
plot-transferfunction.pl : Plot transfer function from explicit filtering
options:
        --help            Show this help.
        --cell=           Plot transfer function of this number
        --filename=       Search files with this filename
        --plotter=        Plot with DataGraph (datagraph, dgraph) or gnuplot
        --term=           Gnuplot terminal
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
    print "$_\n" if(/$opt_fileSearch.*-$opt_cell\.dat$/i);
    push @files, $_ if(/$opt_fileSearch.*-$opt_cell\.dat$/i);
}

#==========================================================================

sub find_labels() 
{
    print "Finding labels in $filename\n";
    open(IN,"<$filename") or die ("Error opening $filename\n");
    while (<IN>)
    {
        chomp;
        $line = $_;
        if ($line =~ m/^\#\s*FGR_target\s*=\s*([0-9]+(\.[0-9]+)?)/) { $FGR_target = $1;}
        elsif ($line =~ m/^\#\s*FGR_xy\s*=\s*([0-9]+(\.[0-9]+)?)/) { $FGR_xy = $1;}
        elsif ($line =~ m/^\#\s*FGR_x\s*=\s*([0-9]+(\.[0-9]+)?)/) { $FGR_x = $1;}
    }
    close(IN);

    $FGR_target = sprintf("%.2f", $FGR_target);
    $FGR_xy = sprintf("%.2f", $FGR_xy);
    $FGR_x = sprintf("%.2f", $FGR_x);
    
    print "FGR_target = $FGR_target\n";
    print "FGR_xy = $FGR_xy\n";
    print "FGR_x = $FGR_x\n";    
}

#==========================================================================

sub parse_yrange()
{
    print "yrange = $opt_yrange\n";
    if ($opt_yrange =~ m/\[([-+]?[0-9]+(\.[0-9]+)?)?:([-+]?[0-9]+(\.[0-9]+)?)?\]/i)
    {
        $ymin = $1;
        $ymax = $3;
    }
    if($ymin eq "") {$ymin="-∞";}
    if($ymax eq "") {$ymax="∞";}
    print "ymin = $ymin\n";
    print "ymax = $ymax\n";
}

#==========================================================================

sub write_tmp_datafile_without_comments() 
{
    open(IN,"<$filename") or die ("Error opening $filename\n");
    open(TMP,">$tmpfile") or die ("Error opening $tmpfile\n");
    while (<IN>)
    {
        chomp;
        if ($_ =~ m/^\#/i) {
            # do nothing
        }
        else {
            print TMP "$_\n";
        }
    }
    close(IN);
    close(TMP);
}

#==========================================================================

sub plot_with_datagraph()
{
    print "Plot with DataGraph\n";
    if (check_in_path("dgraph")) 
    {
        parse_yrange();
        write_tmp_datafile_without_comments();
        run_command("dgraph -s $dgraph_script -v cell=$opt_cell -v FGR_target=$FGR_target -v FGR_xy=$FGR_xy -v FGR_x=$FGR_x -v ymin=$ymin -v ymax=$ymax $tmpfile");
        sleep(2);
        # run_command("rm $tmpfile");
    }
    else
    {
        plot_with_gnuplot();
    }
}

#==========================================================================

sub plot_with_gnuplot()
{
    print "Plot with gnuplot\n";
    if(check_in_path("gnuplot"))
    {
        # Make gnuplot script
        open(OUT,">$gnuplot_script") or die ("Error opening $tmpfile\n");
        if (is_mac()) {
            print OUT "set terminal $opt_gnuplot_terminal\n";
        }
        print OUT "set grid\n";
        print OUT "set yrange $opt_yrange\n";
        print OUT "set style line 1 lt 5 lc 2 lw 1\n";
        print OUT "set style line 2 lt 1 lc 3 lw 1\n";
        print OUT "set style line 3 lt 1 lc 7 lw 1\n";
        print OUT "plot '$filename' using 1:2 w l ls 1 title 'target       FGR=$FGR_target',\\\n";
        print OUT "'$filename' using 1:3 w l ls 2      title 'diagonal xy  FGR=$FGR_xy',\\\n";
        print OUT "'$filename' using 4:5 w l ls 3      title 'diagonal x   FGR=$FGR_x'\n";

        if ($opt_gnuplot_terminal =~ m/aqua/i) {
            # do nothing
        }
        else {
            print OUT "pause -1 'Hit enter to continue...'\n";
        }

        close(OUT);

        # execute script
        run_command("gnuplot $gnuplot_script");
    }
    else
    {
        print "ERROR: Could not find plotting program to plot\n";
        exit(1);
    }
    
}

#==========================================================================

sub plot()
{
    if ($opt_plotter =~ m/d(ata)?graph/i) 
    {
        plot_with_datagraph();
    }
    elsif ($opt_plotter =~ m/gnuplot/i)
    {
        plot_with_gnuplot();
    }
    else
    {
        print "ERROR: plotting program not supported. Choose gnuplot or datagraph (for Mac)\n";
    }
}

#==========================================================================
# Main execution
#==========================================================================

print "Plot transferfunction script\n";

parse_commandline();
find_file();
find_labels();
plot();
