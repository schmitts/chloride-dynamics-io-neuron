::--------------------------------------------------------------------------------------------
:: Script to compile dll from .mod files quickly
:: V1.0
:: Author: Christopher Brian Currin
:: Set up:
::	1) Place script in folder with .mod files
::	2) Delete all files from C:\nrn\MOD (create folder if needed) [First run only]
::	3) run script to compile .mod files
::	3.1) changes in .mod files in working directory are reflected
::	3.2) push enter after compilation to place new .dll in working directory
::	3.2.1) ensure any nrn instances using the to-be-overwritten .dll are closed
::	4) any mod files to be removed from .dll should be deleted from C:\nrn\MOD directly
::--------------------------------------------------------------------------------------------

DEL /F /S /Q /A "c:\nrn\MOD\*"

XCOPY *.mod c:\nrn\MOD /y
:: /m only updated files copied
:: /e all subdirectories too
:: /y confirm all
set "var=%cd%"
cd /d c:\nrn\MOD
:: The “/d” parameter is used to change the current drive to a specific folder from another disk volume.
c:\nrn/mingw/bin/sh c:\nrn/lib/mknrndll.sh /c\nrn
XCOPY nrnmech.dll "%var%" /m /y