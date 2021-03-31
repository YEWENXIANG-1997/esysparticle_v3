# installation gide
Auther: sawano
Date: 2021.03.31

## general instalation gide
https://answers.launchpad.net/esys-particle/+faq/719
This link has instructions for installation. However, it doesn't work in most cases.

## How to get esysparticle
```git clone "this repository"```


## prerequierd
- gcc
- python
- automake
- autoconf
- libtool
- mpi (openmpi or mpich)
- boost (>=1.47.0)

### gcc
The important thing about the compiler is that mpi, boost and python are all built with the same compiler.
The recommended compiler is gcc. intel-based compilers may be able to build esysparticle, but I have failed to build it using icpcc and icc. It's worth a try.

In an environment with multiple versions of gcc, problems are likely to occur with library linking. In particular, links related to glibc and lstdc++ can be a problem. In most cases, the build fails because the path to the directory containing the library of the desired version of gcc is not followed, and a different version of the library is called. In such a case, create "???.conf" in "/etc/ld.so.conf.d/", write the directory that contains the gcc library you want to use, and do ldconfig to solve the problem. I think you can also just add the path to LD_LIBRARY_PATH.

###  python
It is easiest to use python, which is installed by default. If you want to use a different python, please note that you need to explicitly specify the python to use when building the boost python library.

###  automake
If you are using CentOS8, you can install it with yum as follows.
```yum -y install automake```
On CentOS7, esysparticle fails to build because libtool is old. Therefore, you need to build automake from source. For example, please refer to the following URL.

https://stackoverflow.com/questions/33278928/how-to-overcome-aclocal-1-15-is-missing-on-your-system-warning

If you build a new automake from source, maybe you don't need to yum autoconf? (I am not sure. If you get a warning, put them in.)

###  autoconf
```yum -y install autoconf```

###  libtool
On CentOS8, 
```yum -y install libtool```

On CentOS7, the build of esysparticle fails due to an older version of libtool. Please build it from the source.

https://noknow.info/it/os/install_libtool_from_source?lang=ja


### mpi
Please install openmpi or mpich. The one that has worked for me so far is openmpi. mpich should also work, but I haven't tried it.
The machines in Sato lab has ifort installed, but when ifort is installed, the intel mpi is also included. And by default, when mpirun/mpiexec is called, the Intel mpi is called. Of course, the Intel mpi is built with the Intel compiler, so as mentioned above, if boost and python are not built with inte, the esysparticle build will fail.
When installing openmpi or mpich, it is recommended that you do not install it in a directory that has a path through it by default, to avoid conflicts with intel's mpi. The best way for now is to create a file that exports the path to mpi with a name like .bashrc_mpi and activate it when building esysparticle or before running it.


### boost
The boost that you put in with yum will usually fail to build the esysparticle, so build it from source.
The important thing is to use gcc to build boost, specifying the python to be used. If you want to install using the default python, you can use boostrap to specify the prefix and then run . /b2 -prefix=xyz install. If you want to specify a different python, specify the path of the python you want to use in booststrap as in bsshell.sh of the root of boost in hpctech4.


## install for esysparticle
- First, run ./autogen.sh. At this point, automake will check for the required libraries and warn you if they are not present. If they are missing, install them with yum. Note that you don't need bzr. If you get a warning "copying.lib' not found in /usr/share/libtool/libltdl'", this means that you have an older version of libtool. Download and build the new libtool source.
- Next, run ./configure to specify the prefix, path to the boost, etc. Refer to instshell.sh in hpctech4. You can also run ./configure -help to see the options. You may get aclocal errors or warnings about missing boost libraries, and you may need to tweak configure directly, but this is pretty environment-dependent, so I can't give you clear instructions. Good luck.
- Finally, use make -j ## of proc & make install to build and install.