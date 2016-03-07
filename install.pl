#! /usr/bin/perl -w
use strict;
use LWP::Simple;

my $bioperl_home = undef;
my $path = `pwd`;
chomp $path;

#Check if BioPerl is installed

my $perl_config = `perl -MBio::Root::Version -e 'print \$Bio::Root::Version::VERSION'`;
if($perl_config =~ /\d.+/){
	print "BioPerl is installed. Setting library path...\n";
	set_path();
}
else{
	print("BioPerl not found, installing BioPerl now.\n");
	#Install BioPerl
	getstore('http://bioperl.org/DIST/BioPerl-1.6.1.tar.gz', "BioPerl-1.6.1.tar.gz");
	system("tar xvfz BioPerl-1.6.1.tar.gz");
	chdir "BioPerl-1.6.1";
	system("perl Build.PL");
	system("./Build install");
	print "BioPerl is installed. Setting library path...\n";
	set_path();
}



#Set Bioperl modules path in PhyloCore script
sub set_path {
	($bioperl_home) = `perl -MBio::Seq -e'print \$INC{"Bio/Seq.pm"}'` =~ /(^\/.*\/)Bio\/Seq/;
	my $text = undef;
	open (FILE,"<",$path."/PhyloCore.pl") || die "Can't open $path.\/PhyloCore.pl. Make sure PhyloCore.pl is in current directory.\n";
	while (<FILE>) {
		$text .= $_;
	}
	close FILE;

	$text =~ s/#use lib "bioperl_home_dir";/use lib "$bioperl_home";/g;

	open (FILE, ">$path/PhyloCore.pl") || die "Can't write $path/PhyloCore.pl\n";
	print FILE $text;
	close FILE;
	print "Done!\n"
}


