WDIR=/home/rcf-40/jqu/staging/lift
cd $WDIR
proj=dosSantos-Mouse-2015

## copy data from smithlab to hpc
rsync -ravml --safe-links --include="*/results_mm9/*.meth"  --exclude-from 'exclude-list.txt' jqu@smithlab.usc.edu:/labdata/methylome/public/${proj}  ./

## Location of tools and reference files
SOURCE=mm9
TARGET=mm10
INDEX=~jqu/panfs/cache/mouse/mm9/mm9-mm10.index
CHAIN=~jqu/panfs/cache/liftOverchains/mm9ToMm10.over.chain.gz
TARGETSIZES=~jqu/panfs/cache/mouse/mm10/mm10.chrom.sizes
METHPIPE=~jqu/panfs/tools/git/methpipe/bin/
SCRIPT_DIR="/home/rcf-40/jqu/cmb-01/tools/script/WGBS/"
TOOL_DIR="/home/rcf-40/jqu/cmb-01/tools/GenomeBrowser/"
export PATH=${PATH}:${METHPIPE}:${TOOL_DIR}

## lift methylome
echo $proj
# ===================
# lift methcounts
# ===================
for file in `find  ${WDIR}/${proj} -type f -name "*.meth" | grep results_${SOURCE}`; do
  SRCDIR=$(dirname $file)
  TRGTDIR=${SRCDIR/$SOURCE/$TARGET};
  mkdir $TRGTDIR
  cd $TRGTDIR
  touch lifted_${SOURCE}_${TARGET}
  for SAMPLE in `ls $SRCDIR/*.meth`;  do
    NAME=$(basename $SAMPLE)
    NAME=${NAME/.meth/}
    # ===================
    # liftover methcounts
    # ===================
    # make sure conform to new format
    ${SCRIPT_DIR}/to_methcount ${SAMPLE} ${NAME}.meth.${SOURCE}
    # fast liftover
    fast-liftover -p -v -i ${INDEX} -f ${NAME}.meth.${SOURCE} -t ${NAME}.meth.lift2${TARGET}
    # sort
    LC_ALL=C sort -k1,1 -k2,2g -k3,3 -o ${NAME}.meth.sorted ${NAME}.meth.lift2${TARGET}
    # only keep 1-1 mapping
    lift-filter -v -u -o ${NAME}.meth ${NAME}.meth.sorted && \
    rm ${NAME}.meth.sorted ${NAME}.meth.lift2${TARGET} ${NAME}.meth.${SOURCE}
    # ===================
    # liftover allelic score
    # ===================
    if [ -f ../results_${SOURCE}/${NAME}.allelic ]; then
      echo "lift allelic"
      awk '{OFS="\t"; split($4,a,":"); 
      print $1,$2, "+", a[1], $5,a[2]}' < ../results_${SOURCE}/${NAME}.allelic  \
      > ${NAME}.allelic.${SOURCE}
      fast-liftover -p -v -i ${INDEX} -f ${NAME}.allelic.${SOURCE} -t ${NAME}.allelic.lift2${TARGET}
      LC_ALL=C sort -k1,1 -k2,2g -k3,3 -o ${NAME}.allelic.sorted ${NAME}.allelic.lift2${TARGET}
      lift-filter -v -u -o ${NAME}.allelic ${NAME}.allelic.sorted && \
      rm ${NAME}.allelic.sorted ${NAME}.allelic.lift2${TARGET} ${NAME}.allelic.${SOURCE}
      awk '{OFS="\t"; print $1,$2,$2+1, $4":"$6, $5}' < ${NAME}.allelic  > ${NAME}.allelic.tmp && \
      mv ${NAME}.allelic.tmp ${NAME}.allelic 
    fi
    # ===================
    # Call HMR PMD levels
    # ===================
    hmr ${NAME}.meth -v -p ${NAME}.hmrparams -o ${NAME}.hmr
    pmd ${NAME}.meth -v -p ${NAME}.pmdparams -o ${NAME}.pmd
    levels ${NAME}.meth -o ${NAME}.levels
    # ===================
    # liftover AMR
    # ===================
    if [ -f ../results_${SOURCE}/${NAME}.amr ]; then
       echo "lift AMR"
       liftOver ../results_${SOURCE}/${NAME}.amr ${CHAIN} ${NAME}.amr unlifted && rm unlifted
       LC_ALL=C sort -k1,1 -k2,2g -k3,3g  ${NAME}.amr -o ${NAME}.amr.sorted
       bedtools merge -i ${NAME}.amr.sorted > ${NAME}.amr && rm ${NAME}.amr.sorted
    else
       echo "No results_${SOURCE}/${NAME}.amr to liftover"
    fi
    # ===================
    # copy bsrate
    # ===================
    if [ -f ../results_${SOURCE}/${NAME}.bsrate ]; then
      cp ../results_${SOURCE}/${NAME}.bsrate ./
    fi
    # ===================
    # Make tracks
    # ===================
    mkdir  ../tracks_${TARGET}
    ${SCRIPT_DIR}/meth_to_readwig.sh ${NAME}.meth | \
    wigToBigWig /dev/stdin ${TARGETSIZES} ../tracks_${TARGET}/${NAME}.read.bw
    ${SCRIPT_DIR}/meth_to_methwig.sh ${NAME}.meth | \
    wigToBigWig /dev/stdin ${TARGETSIZES} ../tracks_${TARGET}/${NAME}.meth.bw

    if [ -f ${NAME}.allelic ]; then
      cut -f 1-3,5 ${NAME}.allelic | \
      wigToBigWig /dev/stdin ${TARGETSIZES} ../tracks_${TARGET}/${NAME}.allelic.bw
    fi

    for suffix in `echo amr pmd hmr`; do
      if [ -f ${NAME}.${suffix} ]; then
        echo "making track for" $suffix
        cut -f 1-3 ${NAME}.${suffix} > ${NAME}.${suffix}.tmp && \
        bedToBigBed ${NAME}.${suffix}.tmp ${TARGETSIZES} ../tracks_${TARGET}/${NAME}.${suffix}.bb && \
        rm ${NAME}.${suffix}.tmp
      fi
    done
  done
  
  # create symbolic link for higher level results and tracks
  cd  ${WDIR}/${proj} 
  for sampledir in `ls -d ${PWD}/*_*/`; do 
    cd $sampledir
    if [ -d results_mm10 ]; then
      echo "exist"
    else
      echo "not exist"
      DIR=`find . -type d -name "results_mm10" | awk '{print length($0), $0}' | sort -k1,1g | head -1 | awk '{print $2}'`
      mkdir results_mm10
      cd results_mm10
      for i in `ls "."${DIR}/*.*`; do 
        f=$(basename $i)
        f=${f/_R[1-9]/}
        f=${f/_L[1-9]/}
        ln -s $i $f
      done 
      cd ../
      DIR=`find . -type d -name "tracks_mm10"| awk '{print length($0), $0}' | sort -k1,1g | head -1 | awk '{print $2}'`
      mkdir tracks_mm10
      cd tracks_mm10
      for i in `ls "."${DIR}/*.*`; do 
        f=$(basename $i)
        f=${f/_R[1-9]/}
        f=${f/_L[1-9]/}
        ln -s $i $f
      done 
    fi
  done
done

