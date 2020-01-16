#! /usr/bin/perl

$NumArgs = $#ARGV + 1;
$TarBall = $ARGV[0];
$dt = "";
if($NumArgs>1) { $dt = $ARGV[1]; }

$OriginalMetadata = $TarBall;
$OriginalMetadata =~ s/.minimal.tgz//;
system("mkdir ${OriginalMetadata}Reconstructed");
chdir("${OriginalMetadata}Reconstructed")
    or die "\nCouldn't create the directory ${OriginalMetadata}Reconstructed\n\n";
system("tar -xvzf ../${TarBall}");
system("ln -sf ../ReconstituteGrid");

$MetadataOut = $OriginalMetadata;
$MetadataIn  = ${MetadataOut} . ".minimal";

# Handle minor disagrements between the metafile
# name and the archive name
@files = <*>;
foreach $file (@files) {
    if ($file =~ m/$MetadataIn/ix) {
        $MetadataIn = $file;
    }
}

open(MetaFileIn,  "<${MetadataIn}");
open(MetaFileOut, ">${MetadataOut}");

while(<MetaFileIn>) {
    my($line) = $_;
    chomp($line);
    if($line =~ m/^[^\#]*\[ht-ampphi-data\]/) { # Look for the [ht-ampphi-data] section
        print MetaFileOut "[ht-data]\n";
        while(<MetaFileIn>) {
            my($line) = $_;
            chomp($line);
            if($line =~ m/(.*)= *(.*).minimal */) {
                print MetaFileOut "$1 = $2\n";
                print("ReconstituteGrid $2.minimal $dt\n");
                system("ReconstituteGrid $2.minimal $dt") == 0
                    or die "\nCouldn't execute ReconstituteGrid.\nIt might not be on the path.\nYou might want to link it into this directory.\n\n";
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
