#!/usr/bin/perl -w
#
# Example code that separates the individual solutions in the output file
# from adaf_family.pl.

$input=$ARGV[0];


# Opens the output file from the dynamics code
open (INFILE, $input) || 
	die "Can't open $input !";
	
	
	
$i=0;
$stop=0;
while ($stop==0) {

    while (<INFILE>) {    
      print "$. ";
      if (eof(INFILE)) { $stop=1; }  
      if ($_ =~ /eigenvalue/) { $i++; }  
      if ($_ =~ /eigenvalue/ && $i>1) { 
        print "Section $i \n";
        last; 
      }
    }

}

close(INFILE);