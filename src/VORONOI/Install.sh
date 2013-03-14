# Install/unInstall package files in LAMMPS
# edit 2 Makefile.package files to include/exclude package info

if (test $1 = 1) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*voronoi[^ \t]* //g' ../Makefile.package
    sed -i -e 's|^PKG_SYSINC =[ \t]*|&$(voronoi_SYSINC) |' ../Makefile.package
    sed -i -e 's|^PKG_SYSLIB =[ \t]*|&$(voronoi_SYSLIB) |' ../Makefile.package
    sed -i -e 's|^PKG_SYSPATH =[ \t]*|&$(voronoi_SYSPATH) |' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    sed -i -e '/^include.*VORONOI.*$/d' ../Makefile.package.settings
    # multiline form needed for BSD sed on Macs
    sed -i -e '4 i \
include ..\/VORONOI\/Makefile.lammps
' ../Makefile.package.settings
  fi

  cp -p compute_laguerre_atom.h ..
  cp -p compute_voronoi_atom.h ..
  cp -p fix_membrane_gran.h ..
  cp -p compute_laguerre_atom.cpp ..
  cp -p compute_voronoi_atom.cpp ..
  cp -p fix_membrane_gran.cpp ..

elif (test $1 = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*voronoi[^ \t]* //g' ../Makefile.package
  fi

  if (test -e ../Makefile.package.settings) then
    sed -i -e '/^include.*VORONOI.*$/d' ../Makefile.package.settings
  fi

  rm -f ../compute_laguerre_atom.h
  rm -f ../compute_voronoi_atom.h
  rm -f ../fix_membrane_gran.h
  rm -f ../compute_laguerre_atom.cpp
  rm -f ../compute_voronoi_atom.cpp
  rm -f ../fix_membrane_gran.cpp
fi
