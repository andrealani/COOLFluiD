#!/usr/bin/perl
#TODO: 
#add the use of a SRC_DIR to write the files to SRC_DIR/namespace1/.../
#consider option of putting them in the cwd or not.

use strict;
#use warnings;

sub usage() {
  print <<'END_OF_USAGE';
usage:
cpphead.pl [-hucs] namespace1::namespace2::class1 [namespace3::baseClass]
           
END_OF_USAGE
}

sub help() {
  &usage();
  print <<'END_OF_HELP';

Alternatively, when your perl installation is not located in /usr/bin, do
perl test ...

This script writes the basic stuff in a header and implementation file
for C++. 
The following things are set:
header file:
 + ifdef statement
 + include statement for base class if appropriate
 + namespace calls, for this class and the base class if appropriate
 + class name
 + public default constructor, with documentation
 + public default destructor, with documentation
 + private keyword
implementation file:
 + include statement for the header file
 + namespace calls
 + empty implementation for default constructor
 + empty implementation for default destructor

Pieter De Ceuninck, 2002

Usage:
1. A simple class, no inheritance:
 cpphead.pl [-hucs] [namespace1::[namespace2::]]class1
2. A class deriving from a base class (the latter doesn't have to exist):
 cpphead.pl [-hucs] [ns1::[ns2::]]derivedClass [ns3::[ns4::]]baseClass

The following files will be created (if one of them already exists, 
the script will come up with a question whether to overwrite or not):
 classname.hh
 classname.cxx

Mind the fact that you should _not_ type "cpphead.pl classname.hh".

Options:
  -u, --usage                   Show usage.
  -h, --help                    Show this help message.
  -c, --current-dir             Write the files in the current directory. 
                                -c and -s are mutually exclusive.
  -s, --source-dir srcDir       Write the files in a directory, starting from
                                srcDir following the path, given the 
                                namespaces. 
                                -c and -s are mutually exclusive.


END_OF_HELP
}

#function to write a string to the .hh file.
# just use    h   to print an empty line,
# use    h "my string"  to print   my string   on a new line. 
# don't bother about the use of   *   it's ok thx to "$1" instead of $1
sub h() {
    print HFILE "$_[0]\n";
}

#insert beginning of block comment to .hh file
sub bh() {
    &h("  /**");
}

#insert a line with comments to .hh file
sub ch() {
    my $str = "   * ";
    $str .= $_[0];
    &h($str);
}

#insert end of block comment to .hh file
sub eh() {
    &h("   */");
}

#function to write a separation line to the .hh file
sub hsep() {
	&h;
	&h("//////////////////////////////////////////////////////////////////////////////");
	&h;
}

#function to write a string to the .cxx file
# use: @see h()
sub x() {
    print CFILE "$_[0]\n";
}

#function to write a separation line to the .cxx file
sub xsep() {
    &x;
	&x("//////////////////////////////////////////////////////////////////////////////");
    &x;
}

##############################################################################

my $author = `which ypmatch 2>/dev/zero > /dev/zero && ypmatch \`whoami\` passwd || grep \`whoami\` /etc/passwd`;
$author =~ s/([^:]*):([^:]*):([^:]*):([^:]*):([^,:]*).*/\5/;
# | sed -e "s/\([^:]*\):\([^:]*\):\([^:]*\):\([^:]*\):\([^,^:]*\).*/\5/"`;
chop $author;

my $curDirYN = 1; #use current directory for output or not. Default: current.
my $srcDirYN = 0; #don't use a src dir to start from, use current dir by def
my $curDir = "";
my $srcDir = "";

my $target1 = "";
my $target2 = "";

my $c1 = ""; # first class, must exist
my $c2 = ""; # second class, optional. 
            # If it exists, it's interpreted as the base class $c1 derives from.
my $namespace = ""; 

my $ifdefstat = ""; # to be put after #ifndef and #define   
					# no leading or trailing _ or __ 
my $nsCounter = 1;

if ($#ARGV < 0) {
  &usage();
  exit(0);
}

while (my $arg = shift(@ARGV)) {
    if ($arg =~ /^-u$|^--usage$/) {
        &usage();
        exit(0);
    }
    elsif ($arg =~ /^-h$|^--help$/) {
        &help();
        exit(0);
    }
    elsif ($arg =~ /^-s$|^--source-dir$/) {
		$curDirYN = 0;
		$srcDirYN = 1;
        $srcDir = shift(@ARGV);
        if (! $srcDir =~ /[a-zA-Z]+/) {
            die "wrong parameter for the option -s or --source-dir
use -h for help\n";
        }
    }
    elsif ($arg =~ /^-c$|^--current-dir$/) {
        $curDirYN = 1;
        $curDir = `pwd`;
        chop $curDir;
    }
    else {
        $target1 = $arg;
        if ($arg = shift(@ARGV)) {
            $target2 = $arg;
        }
    }
}

if ($curDirYN && $srcDirYN) {
	die "don't combine current dir and src dir!\n";
}

my @target1Array = split("::", $target1);
$c1 = $target1Array[$#target1Array]; # the first class
pop(@target1Array);
#from now on, target1Array contains all namespaces
# get the name of the hfile
my $hfile = "";
my $cfile = "";

if ($curDirYN) {
	$hfile = $c1.".hh";
	$cfile = $c1.".cxx";
} else {
	my $pathToFiles = $srcDir."/";
	$pathToFiles .= join("/", @target1Array);
	$hfile = $pathToFiles."/".$c1.".hh";
	$cfile = $pathToFiles."/".$c1.".cxx";
}

$ifdefstat .= join("_", @target1Array);
$ifdefstat .= "_".$c1."_hh"; # no trailing _ or __ 
#$ifdefstat = uc($ifdefstat); # keep correct caps


my @target2Array = split("::", $target2);
$c2 = $target2Array[$#target2Array];
pop(@target2Array);
#from now on, target2Array contains all namespaces of the baseclass

#############################################################
# Start
#############################################################
my $moveon = "";
if ( -e $hfile or -e $cfile ) {
    print "Overwrite existing files [y,n]";
    $moveon = <STDIN>;
    chop $moveon; 
    if ($moveon !~ /y/) {
        die "Quit without writing\n" ;
    }
}
if ( $c2 =~ /[a-zA-Z]+/ ) {
    print "Generating basic header and cpp file for class $target1, 
deriving from class $target2\n";
}
else {
    print "Generating basic header and cpp file for class $target1\n"
}

##################################################################
# adding lines to .hh
##################################################################
open(HFILE, ">$hfile");
&h("#ifndef ".$ifdefstat);
&h("#define ".$ifdefstat);
#include statement
if ( $c2 =~ /[a-zA-Z]+/ ) {
    &hsep();
	my $incFile = "";
	if ( $target2 =~ /::/ ) {
		( $incFile = $target2 ) =~ s/::/\//g ;
	} else { 
		$incFile = $c2; 
	}
    &h("#include \"".$incFile.".hh\"");
}
&hsep();

&h("// include statements for header files here");

&hsep();

#following block is incorrect, read about namespaces on the web page.

##check if the namespaces of derived class and base class differ
#if (@target1Array != @target2Array) {
##name the namespaces for the base class
#    foreach $namespace (@target2Array) {
#        &h("namespace $namespace { ");
#    }
#    foreach $namespace (@target2Array) {
#        &h("} // end of namespace $namespace ");
#    }
#    if ( $c2 =~ /[a-zA-Z]+/ ) {
#        &h();
#    }
#    foreach $namespace (@target2Array) {
#        &h("using namespace $namespace;");
#    }
#    &hsep();
#}

# open the namespaces for the derived (=this) class, 
# even if the two namespaces are equal.
foreach $namespace (@target1Array) {
    &h("namespace $namespace { ");
}
if ($#target1Array >= 0) { 
    &hsep();
}
#class comment
&bh();
&ch("This class represents");
&ch();
&ch("\@author $author");
if ( $c2 =~ /[a-zA-Z]+/ ) {
  &ch();
  &ch("\@see $c2" );
}
&ch();
&ch("\$Revision: \$");
&ch("\$Source: \$");
&ch("\$Date: \$");
&eh();
if ( $c2 =~ /[a-zA-Z]+/ ) {
	if ( $target2 =~ /::/ ) {
	    &h("class $c1 : public $target2 { " );
	} else {
		&h("class $c1 : public $c2 { " );
	}
}
else {
    &h("class $c1 { "); 
}
&h();  
&h("public: " );
&h();
#Constructors
&bh();
&ch("Default constructor without arguments");
if ( $c2 =~ /[a-zA-Z]+/ ) {
    &ch("\@see $target2()");
}
&eh();
&h("  $c1();");
&h();
&bh();
&ch("Default destructor");
&eh();
&h("  ~$c1();");
&h();
&h("private: ");
&h();
&h("}; // end of class $c1");
#namespaces
&hsep();
foreach $namespace (@target1Array) {
    &h("} // end of namespace $namespace ");
}
if ($#target1Array >= 0) { 
    &hsep();
}
&h("// #include \"$c1.ci\"");
&hsep();
&h("#endif // $ifdefstat");
close(HFILE);
##################################################################
# adding lines 2-end to file.cxx
##################################################################
open(CFILE, ">$cfile");
&x("#include \"$c1.hh\"");
&xsep();
#open namespaces
foreach $namespace (@target1Array) {
    &x("namespace $namespace { ");
}
if ($#target1Array >= 0) { 
    &xsep();
}

if ($target2 =~ /::/ ) {
	my $usingStatement = "using namespace ";
	foreach $namespace  (@target2Array) {
		$usingStatement .= $namespace."::"
	}
	$usingStatement =~ s/::$//;
	&x("$usingStatement;");
	&xsep();
}
#default constructor: empty implementation
if ( $c2 =~ /[a-zA-Z]+/ ) {
    &x($c1 . "::$c1() : $c2()");
}
else {
    &x($c1 . "::$c1()");
}
&x("{");
&x("}");
&xsep();
#default destructor: empty implementation
&x($c1."::~$c1()");
&x("{");
&x("}");
&xsep();
#namespaces
# TQTQTQ: there is a BUG here although not important: the namespaces should be printed in reverse order.
foreach $namespace (@target1Array) {
    &x("} // end of namespace $namespace ");
}
if ($#target1Array >= 0) { 
    &xsep();
}
close(CFILE);
##################################################################

print "done \n";
