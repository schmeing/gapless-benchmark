#!/bin/bash
cd $1
# Get maximum size from all plots
maxsize=$(cat <(pdfinfo figure1.pdf) <(pdfinfo figure2.pdf) <(pdfinfo figure3.pdf) <(pdfinfo figure4.pdf) | grep 'Page size:' | awk '{x=x<$3?$3:x; y=y<$5?$5:y}END{print x,y}')

# Extend all plots to maximum size (with whitespace, no scaling: All plots are expected to be at the same scalling and the size difference caused by trimming)
for fi in figure[1-4].pdf; do
    cropvals=$(cat <(pdfinfo $fi) | grep 'Page size:' | awk -v maxval="$maxsize" '{split(maxval,mv," ");printf "%i %i %i %i", (mv[1]-$3)/2, (mv[2]-$5)/2, (mv[1]-$3+1)/2, (mv[2]-$5+1)/2}');  pdfcrop --margins "$cropvals" $fi c$fi 1>/dev/null
done

# Merge plots to have final width of Landscape A4: 841.89 x 595.276 pts
pdfjam cfigure[1-4].pdf --delta '10px 20px' --papersize '{842px,5000px}' --nup 2x2 --outfile final.pdf 2>/dev/null
pdfcrop --margins '0 20 0 0' final.pdf cfinal.pdf 1>/dev/null

# Add letters (extend all sides to get page size pdf with letter at correct position: left, top, right, bottom)
echo "a" | enscript -B -f LMSans10-Bold20 -o- 2>/dev/null | ps2pdf - > a.pdf 
pdfcrop a.pdf ca.pdf 1>/dev/null #a(11 x 12 pts)
# $3: page width, $5: page height, l[3]: letter width, l[5]: letter height, cor: correction for letters (like a and c) that do need 4 pts extra at the top, pw: panel width, ph: panel height, px: x position, py, y position
cropvals=$(cat <(pdfinfo cfinal.pdf) <(pdfinfo ca.pdf) | grep 'Page size:' | awk '{cor=4; getline line; split(line,l," "); pw=($3+10)/2; ph=$5/2; px=0; py=cor; printf "%i %i %i %i", px, py, $3-px-l[3], $5-py-l[5] }'); pdfcrop --margins "$cropvals" a.pdf ca.pdf 1>/dev/null
pdftk ca.pdf background cfinal.pdf output scfinal.pdf

echo "b" | enscript -B -f LMSans10-Bold20 -o- 2>/dev/null | ps2pdf - > b.pdf
pdfcrop b.pdf cb.pdf 1>/dev/null #b(11 x 16 pts)
cropvals=$(cat <(pdfinfo cfinal.pdf) <(pdfinfo cb.pdf) | grep 'Page size:' | awk '{cor=0; getline line; split(line,l," "); pw=($3+10)/2; ph=$5/2; px=pw; py=cor; printf "%i %i %i %i", px, py, $3-px-l[3], $5-py-l[5] }'); pdfcrop --margins "$cropvals" b.pdf cb.pdf 1>/dev/null
pdftk cb.pdf background scfinal.pdf output scfinal2.pdf && mv -f scfinal2.pdf scfinal.pdf

echo "c" | enscript -B -f LMSans10-Bold20 -o- 2>/dev/null | ps2pdf - > c.pdf
pdfcrop c.pdf cc.pdf 1>/dev/null #c(11 x 12 pts)
cropvals=$(cat <(pdfinfo cfinal.pdf) <(pdfinfo cc.pdf) | grep 'Page size:' | awk '{cor=4; getline line; split(line,l," "); pw=($3+10)/2; ph=$5/2; px=0; py=ph+cor; printf "%i %i %i %i", px, py, $3-px-l[3], $5-py-l[5] }'); pdfcrop --margins "$cropvals" c.pdf cc.pdf 1>/dev/null
pdftk cc.pdf background scfinal.pdf output scfinal2.pdf && mv -f scfinal2.pdf scfinal.pdf

echo "d" | enscript -B -f LMSans10-Bold20 -o- 2>/dev/null | ps2pdf - > d.pdf
pdfcrop d.pdf cd.pdf 1>/dev/null #d(11 x 16 pts)
cropvals=$(cat <(pdfinfo cfinal.pdf) <(pdfinfo cd.pdf) | grep 'Page size:' | awk '{cor=0; getline line; split(line,l," "); pw=($3+10)/2; ph=$5/2; px=pw; py=ph+cor; printf "%i %i %i %i", px, py, $3-px-l[3], $5-py-l[5] }'); pdfcrop --margins "$cropvals" d.pdf cd.pdf 1>/dev/null
pdftk cd.pdf background scfinal.pdf output scfinal2.pdf && mv -f scfinal2.pdf scfinal.pdf
