#!/usr/bin/env perl
#this is replacement for ufasta executable
use strict;
use warnings;

sub process_args {
    my %opts = (
        cmd   => undef,
        fstr   => undef,
        Sflag  => 0,
        Vflag  => 0,
        N    => undef,
        Hflag  => 0,
    );

    # First argument must be a word
    die "Usage: $0 <cmd> <options>\nfasta file must be piped in from STDIN\n"
        unless @ARGV;

    $opts{cmd} = shift @ARGV;

    # Process remaining arguments manually
    while (@ARGV) {
        my $arg = shift @ARGV;

        if ($arg eq '-f') {
            $opts{fstr} = shift @ARGV;
            die "Error: file $opts{fstr} does not exist \n" if(not(-e $opts{fstr}));
        } elsif ($arg eq '-N') {
            $opts{N} = shift @ARGV;
            $opts{N} =~ /^\d+$/ or die "Error: -N requires a number\n";
        } elsif ($arg eq '-H') {
            $opts{Hflag} = 1;
        } elsif ($arg eq '-S') {
            $opts{Sflag} = 1;
        } elsif ($arg eq '-v') {
            $opts{Vflag} = 1;  
        } else {
            die "Unknown argument: $arg\n";
        }
    }

    return \%opts;
}

# Example usage:
my $opts = process_args();

if($opts->{cmd} eq "one"){
  my $contline=0;
  while(my $line=<STDIN>){
    chomp($line);
    if($line =~ /^>/){
      print "\n" if($contline);
      $contline=1;
      print $line,"\n";
    }else{
      print $line;
    }
  }
  print "\n" if($contline);
}elsif($opts->{cmd} eq "sizes"){  
  my $seq="";
  my $name="";
  while(my $line=<STDIN>){
    chomp($line);
    next if($line eq '');
    if($line =~ /^>/){
      if(length($name)>0){
        print "$name " if($opts->{Hflag} == 1);
        print length($seq),"\n";
      }
      $name=(split(/\s+/,substr($line,1)))[0];
      $seq="";
    }else{
      $seq.=$line;
    }
  }
  if(length($name)>0){
    print "$name " if($opts->{Hflag} == 1);
    print length($seq),"\n";
  }
}elsif($opts->{cmd} eq "extract"){
  my $flag=0;
  my %names=();
  open(FILE,$opts->{fstr});
  while(my $line=<FILE>){
    chomp($line);
    $line=(split(/\s+/,$line))[0];
    $names{$line}=1;
  }
  while(my $line=<STDIN>){
    chomp($line);
    next if($line eq '');
    if($line =~ /^>/){
      my $name=(split(/\s+/,substr($line,1)))[0];
      $flag=0;
      $flag=1 if(($opts->{Vflag}==0 && defined($names{$name})) || ($opts->{Vflag}==1 && not(defined($names{$name}))));
      print ">$name\n" if($flag);
    }else{
      print $line,"\n" if($flag);
    }
  }
}elsif($opts->{cmd} eq "n50"){

  my %len;
  my $current;

# Read FASTA and accumulate lengths
  while (<>) {
    chomp;
    next if $_ eq '';

    if (/^>(\S+)/) {
      $current = $1;
      $len{$current} = 0;
    } else {
      s/\s+//g;
      $len{$current} += length($_);
    }
  }
  my @lengths = sort { $b <=> $a } values %len;
  my $total   = 0;
  $total += $_ for @lengths;

# Compute N50
  my $frac=0;
  if(not(defined($opts->{N}))){
    $frac=0;
  }else{
    if($opts->{N}>0 && $opts->{N}<=100){
      $frac=$opts->{N}/100;
    }else{
      die "Error: argument for N must be  >0 and <=100";
    }
  }
  my $half = $total * $frac;
  my $running = 0;
  my $n50;

  for my $L (@lengths) {
    $running += $L;
    if ($running >= $half) {
      $n50 = $L;
      last;
    }
  }
  print "Sequence $total\n" if($opts->{Sflag} ==1);
  print "N$opts->{N} $n50\n" if(defined($opts->{N}));
}else{
  die "Error: unknown command\n";
}




