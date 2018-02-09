# cosmolike_light

CosmoLike installation (on JPL HPC, similar and easier for your laptop since you have root rights):

- install fftw3
	a) goto http://www.fftw.org/download.html
	b) download v3.3.4, unpack tar file, copy directory to zodiac into your home directory
	c) goto directory , configure --prefix=/home/YOURUSRNAME --with-pic , make, make install
	
- install gsl 2.1
	a) goto http://www.gnu.org/software/gsl/
	b) download, unpack tar file, copy directory to zodiac into your home directory 
	c) goto directory, configure --prefix=/home/YOURUSRNAME , make, make install
	
- install anaconda
	https://docs.anaconda.com/anaconda/install/linux

- add home/teifler/anaconda2/bin to bash or csh (bash is done automatically) for csh add 
setenv LD_LIBRARY_PATH /home/YOURUSRNAME/lib

- install yaml: type "conda install -c anaconda yaml"
- install numpy: type "conda install -c anaconda numpy"
- install scipy: type "conda install -c anaconda scipy"
- install emcee: type "conda install -c astropy emcee"
