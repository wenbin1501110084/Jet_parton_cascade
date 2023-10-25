This is the framework to couple the jet parton shower and parton cascade

## Install the LHAPDF ##
version="6.5.4"
echo "Building LHAPDF ${version} ... "
tar -xf utilities/LHAPDF-${version}.tar.gz
LIBPATH=`echo $PWD/LHAPDF_Lib`
(
    cd LHAPDF-${version}/
    ./configure --prefix=`echo $LIBPATH` --disable-python
    make -j${num_of_cores}
    make install
)
rm -fr LHAPDF-${version}

echo "downloading pdfsets ... "
wget http://lhapdfsets.web.cern.ch/lhapdfsets/current/NNPDF31_nnlo_as_0118.tar.gz -O- | tar xz -C $PWD/LHAPDF_Lib/share/LHAPDF

## Then compile pythia8 with LHAPDF ##
./configure --prefix=/wsu/home/he/he92/he9215/soft_ware/Pythia8_package/pythia8310_install --with-lhapdf6=/wsu/home/he/he92/he9215/soft_ware/LHAPDF_package/LHAPDF_install

Then please remember to copy the Makefile and Makefile.inc in the example folder to the pythia parton folder

