#! /usr/bin/perl

# Use this script with a metadata file as
# 
#   MinimizeFullWaveform.pl MetadataFileName
# 
# The metadata file is searched for an [ht-data] section.
# The data file names found in that section are converted
# to a minimal time grid, and output with file names
# appended by '.minimal'.  A new metadata file named
# MetadataFileName.minimal is produced, referring to
# these minimal files.  Finally, the whole set is tarred
# into an archive named MetadataFileName.minimal.tgz.
# The new metadata file and tarball can be submitted to
# the repository.
# 
# Note that the 'MinimizeGrid' routine in this directory
# needs to be in the path when this file is run.  You
# might want to link it from the data directory.
# 
# Also note that there are two optional arguments to the
# 'MinimizeGrid' function.  These can be passed as
# trailing arguments to this script in the same way,
# and will be used on each of the minimized data files.



$NumArgs = $#ARGV + 1;
$MetadataIn = $ARGV[0];
$AmpTol = "";
$PhiTol = "";
if($NumArgs>1) { $AmpTol = $ARGV[1]; }
if($NumArgs>2) { $PhiTol = $ARGV[2]; }

$MetadataOut = ${MetadataIn} . ".minimal";
$TarCommand = "tar -czvf ${MetadataIn}.minimal.tgz " . $MetadataOut;

open(MetaFileIn,  "<${MetadataIn}");
open(MetaFileOut, ">${MetadataOut}");

while(<MetaFileIn>) {
    my($line) = $_;
    chomp($line);
    if($line =~ m/^[^\#]*\[ht-data\]/) { # Look for the [ht-data] section
        print MetaFileOut "[ht-ampphi-data]\n";
        while(<MetaFileIn>) {
            my($line) = $_;
            chomp($line);
            if($line =~ m/(.*)= *(\S*)/) {
                print MetaFileOut "$1 = $2.minimal\n";
                print("MinimizeGrid $2 $AmpTol $PhiTol\n");
                system("MinimizeGrid $2 $AmpTol $PhiTol\n") == 0
                    or die "\nCouldn't execute MinimizeGrid.\nIt might not be on the path.\nYou might want to link it into this directory.\n\n";
                $TarCommand .= " $2.minimal";
            } else {
                print MetaFileOut "$line\n";
            }
            if($line =~ m/\[/) { last; }
        }
    } else {
        print MetaFileOut "$line\n";
    }
}

close(MetaFileIn);
close(MetaFileOut);

system($TarCommand);
