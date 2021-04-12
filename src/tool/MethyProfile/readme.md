## MethyProfile: 

----------

### What is MethyProfile:

MethyProfile is focus on generating methylation profile based on methylation information (BED file) and genome information (GFF3 file), which is upstream of the methylation analysis process.

-----------

### Command Line:

usage:

```bash
methyprofile [options] input.bed input.gff3 (output.methyprofile.txt)
```
options:
```
   -P, --promoter=n    - Analysing promoters methylation information with n bp. The default is 2000bp.
   -l, --single-list   - Gathering single gene information for a gene list behind.
   -h, --help          - Show this message.
```

-----------

### Building from source:

To build MethyProfile manually, please make sure cmake on your UNIX system (MacOS or LINUX) computer.

Run these command below:

```
mkdir build
cd build
cmake ..
```

Then you can execute `./methyprofile` in current directory.

Besides, to enable plug-in saving, execute command below:

```
cd etc
./build
cp -R ./etc ./build/  # Copy ./etc to the dirctory of methyprofile 
```

-----------------

### More:

You can change macros in `config.h` and rebuild the project to meet the actual needs. For more information, please see code comment in `config.h`.
