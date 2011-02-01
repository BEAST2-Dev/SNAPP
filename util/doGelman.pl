#!/usr/bin/perl

$nFiles = $#ARGV+1;

if ($nFiles <= 1) {
    print "Usage: $0 [-percent XX] <file1.log> <file2.log> ... <filen.log>\n";
    print "where fileX.log a log file produced by snap/beast2\n";
    print "-percent XX : XX percent of the data (excluding burn in) will be used\n";
    print "Prints Gelman statistics for all columns in the log files\n";
    print "This assumes all log files have the same nr of columns and same nr of samples\n";
    exit;
}



$nSamplePercentage = 100;
$nFileStart = 0;
if ($ARGV[0] eq '-percent') {
    $nSamplePercentage = $ARGV[1];
    $nFileStart = 2;
    $nFiles -= 2;
}
doGelman();



sub doGelman() {
    $sHeader = "";
    for ($iFile = 0; $iFile < $nFiles; $iFile++) {
        $sFile = $ARGV[$nFileStart + $iFile];
        print "Processing $sFile\n";
        open(FIN,$sFile) or die "Cannot open file >$sFile< for reading";
        # count samples
        $nSamples = 0;
        $nData = 0;
        while ($s=<FIN>) {
            if (($sHeader eq "") && ($s !~ /#/) && ($s =~ /[a-zA-Z]/)) {
                $sHeader = $s;
            }
            if ($s =~ /^([0-9]+)\t/) {
                $nSamples = $1;
                $nData++;
            }
        }
        close FIN;
        $nBurnIn = ($nSamples/10)*($nSamplePercentage/100);
        $nMaxData = $nSamples*$nSamplePercentage/100;


        # process file
        open(FIN,$sFile) or die "Cannot open file >$sFile< for reading";
        # count samples
        $nData = 0;
        $nCols = 0;
        while ($s=<FIN>) {
            if ($s =~ /^([0-9]+)\t/) {
                $nSamples = $1;
                if (($nBurnIn <= $nSamples) && ($nSamples <= $nMaxData)) {
                    @s = split('\t',$s);
                    $nCols = $#s-1;
                    for($i=1;$i<=$#s;$i++) {
                        $fMean[$iFile][$i] = $fMean[$iFile][$i] + $s[$i];
                        $fVar[$iFile][$i]  = $fVar[$iFile][$i]  + $s[$i] * $s[$i];
                    }
                    $nData++;
                }
            }
        }
        close FIN;
    }

print "Data points=$nData\n";

    for ($iFile = 0; $iFile < $nFiles; $iFile++) {
        for($i=1;$i<=$nCols;$i++) {
            $fMean[$iFile][$i] = $fMean[$iFile][$i] / $nData;
            $fVar[$iFile][$i]  = $fVar[$iFile][$i]  / $nData - $fMean[$iFile][$i] * $fMean[$iFile][$i];
        }
    }

    #x = sum(M[i],i=1..m)
    for($i=1;$i<=$nCols;$i++) {
        for ($iFile = 0; $iFile < $nFiles; $iFile++) {
            $x[$i] = $x[$i] + $fMean[$iFile][$i];
        }
        $x[$i] = $x[$i]/$nFiles;
    }

    #B =  sum((M[i] - x)^2, i=1..m) / (m-1)
    for($i=1;$i<=$nCols;$i++) {
        for ($iFile = 0; $iFile < $nFiles; $iFile++) {
            $B[$i] = $B[$i] + ($fMean[$iFile][$i] - $x[$i])*($fMean[$iFile][$i] - $x[$i]) / ($nFiles-1);
        }
    }

    #W = sum(V[i],i=1..m)/m
    for($i=1;$i<=$nCols;$i++) {
        for ($iFile = 0; $iFile < $nFiles; $iFile++) {
            $W[$i] = $W[$i] + $fVar[$iFile][$i] / $nFiles;
        }
    }


    @sHeader = split("\t",$sHeader);
    print "\n(mean,variance) summary\nvariable";
    for ($iFile = 0; $iFile < $nFiles; $iFile++) {
        print " $ARGV[$nFileStart + $iFile]";
    }
    print "\n";
    for($i=1;$i<=$nCols;$i++) {
        print "$sHeader[$i] ";
        for ($iFile = 0; $iFile < $nFiles; $iFile++) {
            print  "$fMean[$iFile][$i],$fVar[$iFile][$i] ";
        }
        print "\n";
    }


    #R = ((n-1)/n*W + B ) / W
    for($i=1;$i<=$nCols;$i++) {
        if ($W[$i]>0) {
            $R[$i] = ((($nData-1.0)/$nData) * $W[$i] + $B[$i]) / $W[$i];
        }
    }

    print "\nGelman statistics\n";
    for($i=1;$i<=$nCols;$i++) {
        print "$sHeader[$i] ";
        print sqrt($R[$i]);
        print "  B=".$B[$i] * $nData;
        print "  W=".$W[$i];
        print "\n";
    }
    print "\n";
}
