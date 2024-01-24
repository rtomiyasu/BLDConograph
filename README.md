# Instructions for BLDConograph
[BLDConograph version 1.0](https://github.com/rtomiyasu/BLDConograph/tree/main/BLDConograph1_0_03_win) is a standalone program to execute the Bravais lattice determination by the method introduced in the following reference (download page in SourceForge.jp; Example of output files).

This code is also used in the powder auto-indexing software Conograph. (Conograph GUI version).

## NEWS
### 2016/8/6
- The output format for base-centered monoclinic cells was corrected.
### 2014/2/16
- The following bugs were fixed in the newest version. (Note that the corresponding part of the above paper is correct.)
One of the three matrices used for face- and body-centrings was incorrect.
In the output file, the value between <Distance> and </Distance> was replaced by its square root. (The former should be called the "norm".)
- The program did not compute the distances (and the norms) properly, when the input cell and the primitive cell of the output convensional cell are not reduced with regard to the same basis.

## How to use the BLDConograph program.
1. BLDConograph requires the following two input files (examples can be found in the sample folder):
    - A *.inp.xml file that includes information about input parameters ([Example](https://github.com/rtomiyasu/BLDConograph/blob/main/BLDConograph1_0_03_win/sample/sample1(Tetragonal(I)_3.84%2C3.84%2C20.1%2C90%2C90%2C90)/HERMES_Sr327_250K.inp.xml))
    - A cntl.inp.xml file that includes the names of the *.inp.xml file and the output file ([Example](https://github.com/rtomiyasu/BLDConograph/blob/main/BLDConograph1_0_03_win/sample/sample1(Tetragonal(I)_3.84%2C3.84%2C20.1%2C90%2C90%2C90)/cntl.inp.xml))
1. Copy one of the folders from the sample folder. Modify the contents of the two xml-files and the name of the *.inp.xml file if necessary. If you change the name of the *.inp.xml file, then it will be necessary to modify the contents of the cntl.inp.xml file accordingly.
1. Open a command prompt or terminal window in your operating system. Change the current folder to the folder that contains the modified cntl.inp.xml file.
1. Enter the absolute path to the BLDConograph.exe file on the command line and execute BLDConograph.