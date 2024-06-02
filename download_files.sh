DOWNLOAD_SIZE=.200MB
CUT_SIZE=128MiB
mkdir data
curl -o data/sources.gz https://pizzachili.dcc.uchile.cl//texts/code/sources$DOWNLOAD_SIZE.gz
curl -o data/english.gz https://pizzachili.dcc.uchile.cl//texts/nlang/english$DOWNLOAD_SIZE.gz
curl -o data/dna.gz https://pizzachili.dcc.uchile.cl//texts/dna/dna$DOWNLOAD_SIZE.gz

(
cd data
gzip -d sources.gz
gzip -d english.gz
gzip -d dna.gz
head -c $CUT_SIZE dna > dna.tmp && mv dna.tmp dna
head -c $CUT_SIZE english > english.tmp && mv english.tmp english
head -c $CUT_SIZE sources > sources.tmp && mv sources.tmp sources
)