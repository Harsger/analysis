filename=$1
specifier=$2
destination=$3
hadd -k ${specifier}_${filename}_results.root ${filename}*${specifier}*
if [ "$?" != "0" ]
then
	echo " ERROR : hadd failed "
	exit 1
fi
mv ${specifier}_${filename}_results.root $destination
if [ "$filename" != "" ]
then
	rm ${filename}*${specifier}*
fi
