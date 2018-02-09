# cosmolike_light

CosmoLike installation (on JPL HPC, similar and easier for your laptop since you have roo rights):

- install fftw3
	a) goto http://www.fftw.org/download.html
	b) download v3.3.4, unpack tar file, copy directory to zodiac into your home directory
	c) goto directory , configure --prefix=/home/teifler --with-pic , make, make install
- install gsl 2.1
	a) goto http://www.gnu.org/software/gsl/
	b) download, unpack tar file, copy directory to zodiac into your home directory 
	c) goto directory , configure --prefix=/home/teifler , make, make install
- install anaconda
	https://docs.anaconda.com/anaconda/install/linux

- add home/teifler/anaconda2/bin to bash or csh (bash is done automatically) for csh add 
setenv LD_LIBRARY_PATH /home/teifler/lib

- install yaml: conda install -c anaconda yaml
- install numpy: conda install -c anaconda numpy
- install scipy: conda install -c anaconda scipy
- install emcee: conda install -c astropy emcee
- install matplotlib: conda install -c conda-forge matplotlib
