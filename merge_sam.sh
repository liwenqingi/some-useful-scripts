header=0.sam
files="0.sam 1.sam 2.sam"
output=out.sam
 
(grep ^@ $header; for f in $files; do grep -v ^@ $f; done) > $output
