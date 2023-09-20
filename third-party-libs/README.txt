# Third-Party Libraries


This director contains a set of libraries that are required by
FCLib. Rather than selectively distributing files from these
libraries, we have opted to include the original tar files in this
release and use patches to make modifications that FCLib requires. You
should be able to build these libraries simply by doing a "make" in
the root directory of this bundle. That make will unpack the
libraries, patch them, run configure, and do a build. All build files
are installed in the "bundle" directory.

Once you have built these TPLs, edit fclib's `Makefile.include` and
point NETCDF_HOME at the bundle directory (use a full path).

Please see the release information inside the tar files for
information about the individual libraries.



