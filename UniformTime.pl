#! /usr/bin/perl

$NumArgs = $#ARGV + 1;
$BBHFile = $ARGV[0];
$dt = "0.25";
if($NumArgs>1) { $dt = $ARGV[1]; }

system("cp ${BBHFile} ${BBHFile}.original");

open(MetaFileIn,  "<${BBHFile}.original");
$BBHFile =~ s/.minimal//;
open(MetaFileOut, ">${BBHFile}");

while(<MetaFileIn>) {
    my($line) = $_;
    chomp($line);
    if($line =~ m/^[^\#]*\[ht-ampphi-data\]/) { # Look for the [ht-ampphi-data] section
        print MetaFileOut "[ht-data]\n";
        while(<MetaFileIn>) {
            my($line) = $_;
            chomp($line);
            if($line =~ m/(.*?,.*?)= *([^\#]*)(.*)/) {
                $Mode = $1;
                $File = $2;
                $Comment = $3;
                $File =~ s/\.gz//;
                $NewFile = $File;
                if($File =~ m/minimal/) {
                    $NewFile =~ s/.minimal//;
                } else {
                    $NewFile .= ".reconstituted";
                }
                print MetaFileOut "$Mode = $NewFile $Comment\n";
                system("ReconstituteGrid $File $dt");
                print("Reconstituted File: $File \n");
            } else {
                print MetaFileOut "$line\n";
            }
            if($line =~ m/^[^\#]*\[/) { last; }
        }
    } else {
        print MetaFileOut "$line\n";
    }
}

close(MetaFileIn);
close(MetaFileOut);
