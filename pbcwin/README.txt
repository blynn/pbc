THE PBC LIBRARY - VISUAL C++ EDITION

* NOTE: This project has prerequisites. Read the "COMPILING" section before
        attempting to build the solution.

This is the Visual C++ version of the PBC library. It is capable of compiling
the core library in shared or dynamic, and 32- or 64-bit versions. The Visual
C++ solution is also able to compile the various binaries included with PBC,
but only shared library configurations are provided for these binaries. These
projects were created using Visual C++ 2013.


DIRECTORIES
------------------------------------------------------------------------------
The following directories are used solely by the Visual C++ edition:

   pbcwin          - Root directory for all Visual C++ data
   pbcwin/projects - Contains projects for the library and all binaries
   pbcwin/bin      - Binaries for examples, tests, and the calculator
   pbcwin/dll      - Dynamic libraries
   pbcwin/lib      - Shared libraries
   pbcwin/obj      - All intermediate build files


COMPILING
------------------------------------------------------------------------------
There are several requirements for compiling PBC in Visual C++:

   MPIR  - Windows port of GMP; required for the core library
   Bison - Required for compiling the "PBC" calculator application
   Flex  - Required for compiling the "PBC" calculator application

To build anything, you will need to install MPIR. If you wish to build the PBC
calculator application, you will also need to install Bison and Flex. See
below for details.


INSTALLING MPIR
------------------------------------------------------------------------------
To install MPIR, visit http://mpir.org/ and download the latest sources. While
you may extract them anywhere, the PBC project files expect them to be in a
directory named "mpir" at the same level as your PBC directory.

For example, if you have placed PBC in:

   C:\Users\Username\Documents\pbc

then you will need to extract the MPIR sources in:

   C:\Users\Username\Documents\mpir

If you would prefer to place MPIR in a different location, then you will need
to modify the VC++ include directories in each project file manually.

Once you have extracted MPIR, you must compile it. If you are using Visual
C++ 2013, open "mpir\build.vc12\mpir.sln". The MPIR solution contains many
projects for different processor architectures (the code is highly optimized).
The PBC library will use the appropriate MPIR projects (DLL vs lib),
configurations (Debug vs Release), and platform (Win32 vs x64) based on your
settings in the IDE. It is recommended that you compile all eight versions for
your specific processor architecture.

For example, if you are using an Intel CPU with the Sandy Bridge architecture,
you would perform the following steps:

   1) Compile "dll_mpir_sandybridge" in 32-bit debug mode
   2) Compile "lib_mpir_sandybridge" in 32-bit debug mode
   3) Compile "dll_mpir_sandybridge" in 32-bit release mode
   4) Compile "lib_mpir_sandybridge" in 32-bit release mode
   5) Compile "dll_mpir_sandybridge" in 64-bit debug mode
   6) Compile "lib_mpir_sandybridge" in 64-bit debug mode
   7) Compile "dll_mpir_sandybridge" in 64-bit release mode
   8) Compile "lib_mpir_sandybridge" in 64-bit release mode

Note that the MPIR solution only supports building one project at a time.
   
If you are using a different processor architecture, compile the appropriate
projects. If your processor has no optimized builds available, you can instead
compile the "dll_mpir_gc" and "lib_mpir_gc" projects, although you will suffer
a performance penalty.


INSTALLING BISON AND FLEX
------------------------------------------------------------------------------
The PBC calculator application requires the generation of a parser and lexer.
The Bison and Flex tools are used to accomplish this. While these applications
are not required to compile the core PBC library for use in your programs,
you must install them to compile the calculator.

Install Bison for Windows from:
http://gnuwin32.sourceforge.net/packages/bison.htm

Install Flex for Windows from:
http://gnuwin32.sourceforge.net/packages/flex.htm

It is recommended that you install these applications using the setup files.
After they have been installed, you will need to add the binaries to your PATH
environment variable. On recent versions of Windows, open:

   Control Panel > System > Advanced system settings > Advanced
   > Environment Variables

You can now edit the PATH variable for the system to add the Bison and Flex
binaries. For example, if your path is currently set to "MyPath" and you
installed Bison and Flex to "C:\Program Files (x86)\GnuWin32", set PATH to:

   MyPath;C:\Program Files (x86)\GnuWin32\bin

You should now be able to compile the "pbc" project without issue. Note that
this project is only set to compile in the "Debug" and "Release"
configurations.


USING PBC
------------------------------------------------------------------------------
To use the PBC library in your own applications, you must link to it in the
standard manner. Include files are located in the "include" directory, and
libraries are found in "pbcwin\lib\$(Platform)\$(Configuration)". DLLs are
placed in "pbcwin\dll\$(Platform)\$(Configuration)".

For convenience, we reiterate the process for Visual C++ 2013 applications
wishing to statically link against PBC:

In the properties dialog for your new project, open:

   Configuration Properties > VC++ Directories

Open the "Include Directories" field and add the "include" directory for PBC.
For example, if you installed PBC to:

   C:\Users\Username\Documents\pbc

and you installed MPIR to:

   C:\Users\Username\Documents\mpir

and your current include directories value was "$(VC_IncludePath)", set it to:

   C:\Users\Username\Documents\pbc\include;
   C:\Users\Username\Documents\mpir\lib\$(Platform)\$(Configuration);
   $(VC_IncludePath)

(Note that this should appear as one line in the property pages, but it is
split in this document to restrict line length. You can select "Edit..." from
the dropdown menu in the field to add individuals paths using a GUI, if you
would prefer to do so)

Likewise, you need to add the library path to the "Library Directories" field.
For example, if your library path was "$(VC_LibraryPath_x64)", set it to:

   C:\Users\Username\Documents\pbc\lib\$(Platform)\$(Configuration);
   C:\Users\Username\Documents\mpir\lib\$(Platform)\$(Configuration);
   $(VC_LibraryPath_x64)

Note that this will only work if you maintain the default "Debug" and
"Release" configurations for your project. If you change them, manually set
the "$(Configuration)" part of the path to point to the correct directory.

If your project is an executable, then open:

   Configuration Properties > Linker > Input

If your project is a library, then open this instead:

   Configuration Properties > Librarian > general

In either case, modify the "Additional Dependencies" field to include
"pbclib.lib" as a dependency. You should now be able to link against PBC.