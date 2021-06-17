## Project description
This is the IHEP-LQCD collaboration's repository for using [USQCD](https://www.usqcd.org/) [chroma](http://jeffersonlab.github.io/chroma/) library with customized addons.

## Compilation
Requirment:
- For GPU: Chroma (with QDP-JIT, QUDA)
- For CPU: Chroma (with QDP++)
Build & Install:
- setup chroma's environemnt, e.g. `module load lqcd/chroma/double/latest-openmpi` on IHEP's GPU cluster.
 ```
	mkdir build install && pushd build 
	cmake -DCMAKE_INSTALL_PREFIX=../install ..
	make 
	make install
 ```

## Develop
Create a new branch (e.g. feature, bugfix ...) to develop and test the code and then merge into the **devel** branch.
Seek for code review if possiable and merge into **master** branch finally.

## Code format
use `./clang-format-all.sh` to format the code before commit, e.g. `./clang-format-all.sh .` to format all the code under root directory
or specify a directory to format.
