#!/usr/bin/perl


$nFiles = $#ARGV+1;

if ($nFiles <= 1) {
    print "Usage: $0 <file1.log> <file2.log> ... <filen.log>\n";
    print "where fileX.log a log file produced by snap/beast2\n";
    print "Prints Gelman statistics for all columns in the log files\n";
    print "This assumes all log files have the same nr of columns and same nr of samples\n";
    exit;
}

$sHeader = "";
for ($iFile = 0; $iFile < $nFiles; $iFile++) {
    $sFile = $ARGV[$iFile];
    print "Processing $sFile\n";
    open(FIN,$sFile) or die "Cannot open file >$sFile< for reading";
    # count samples
    $nSamples = 0;
    while ($s=<FIN>) {
        if (($sHeader eq "") && ($s !~ /#/) && ($s =~ /[a-zA-Z]/)) {
            $sHeader = $s;
        }
        if ($s =~ /^([0-9]+)\t/) {
            $nSamples = $1;
        }
    }
    close FIN;
    $nBurnIn = $nSamples/10;


    # process file
    open(FIN,$sFile) or die "Cannot open file >$sFile< for reading";
    # count samples
    $nData = 0;
    $nCols = 0;
    while ($s=<FIN>) {
        if ($s =~ /^([0-9]+)\t/) {
            $nSamples = $1;
            if ($nBurnIn <= $nSamples) {
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


for($i=1;$i<=$nCols;$i++) {
    for ($iFile = 0; $iFile < $nFiles; $iFile++) {
        print  "$x[$i],$B[$i],$W[$i] ";
    }
    print "\n";
}


#R = ((n-1)/n*W + B ) / W
for($i=1;$i<=$nCols;$i++) {
    if ($W[$i]>0) {
        $R[$i] = ((($nData-1.0)/$nData) * $W[$i] + $B[$i]) / $W[$i];
    }
}

print $sHeader;
@sHeader = split("\t",$sHeader);
for($i=1;$i<=$nCols;$i++) {
    print "$sHeader[$i] ";
    print sqrt($R[$i]);
    print "\n";
}
print "\n";

